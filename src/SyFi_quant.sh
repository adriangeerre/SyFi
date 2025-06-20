#!/bin/bash

# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Gijs Selten, Florian Lamouche and Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk / g.selten@uu.nl
# Created Date  : 01/10/2024
# version       : '1.0'
# ---------------------------------------------------------------------------
# Pipeline (bash) to perform fingerprint identification from microbiome data.
# ---------------------------------------------------------------------------

function logo() {
	# Colors
	bold=$(tput bold)
	normal=$(tput sgr0)

	# Arrays
	array_keep=("None" "Salmon output" "All")
	array_verbose=("Quiet" "Sample" "All")

	# Variables
	if [[ $1 != "help" ]]; then
		readfolder=$1
		fingerprintfolder=$2
		readtype=$3
		minscorefraction=$4
		threads=$5
		keepfiles=${array_keep[${6}]}
		verbose=${array_keep[${7}]}

		# Logo
		echo ""
		echo "${bold}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${normal}"
		echo "${bold}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${normal}"
		echo "${bold}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${normal}"
		echo "${bold}|____|/___/    /____/   |___|¯¯    |________|${normal}"
		echo ""
		echo "${bold}----------------------------------------------${normal}"
		echo "${bold}Summary:${normal}"
		echo "    Read folder: ${readfolder}"
		echo "    Fingerprint folder: ${fingerprintfolder}"
		echo "    Read type: ${readtype}"
		echo "    Minimum pseudoalignment identity score: ${minscorefraction}"
		echo "    Threads: ${threads}"
		echo "    Keep files: ${keepfiles}"
		echo "    Verbose: ${verbose}"
		echo ""
		printf "${bold}Start:${normal}\n"
		date "+    date: %d/%m/%Y"
		date "+    time: %H:%M:%S"
		echo ""
		echo "${bold}Progress:${normal}"
		checkProgress ${readfolder}
		echo "${bold}----------------------------------------------${normal}"
		echo ""
	else
		# Logo
		echo ""
		echo "${bold}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${normal}"
		echo "${bold}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${normal}"
		echo "${bold}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${normal}"
		echo "${bold}|____|/___/    /____/   |___|¯¯    |________|${normal}"
		echo ""
	fi
}

function usage()
{
	# Colors
	bold=$(tput bold)
	normal=$(tput sgr0)

	# Logo
	logo help

	# Usage
	echo "Usage: SyFi.sh quant -i <READ_FOLDER> -f <FINGERPRINT_FOLDER> -t <THREADS>"
	printf "\n"
	echo "${bold}REQUIRED:${normal}"
	echo "# Input"
	echo "  -i  | --read_folder                      Folder containing the sample reads. The software assumes that the folder contains sub-folders for each sample containing the fastq.gz files of either single end reads or paired end reads."
	echo "  -f  | --fingerprint_folder               SyFi main module output folder that contains the Fingerprints. If the user wants to exclude bacterial members, a personalized folder can be created with the SynCom members of interest using the defined folder structure."
	printf "\n"
	echo "${bold}OPTIONAL:${normal}"
	echo "# Read type"
	echo "  -r  | --read_type                        Paired or single end reads [paired or single] (default: paired)."
	printf "\n"
	echo "# Pseudoalignment percentage identity score:"
	echo "  -m  | --minscorefraction                 Percentage identity score that Salmon uses for pseudoaligning metagenomic reads to the SyFi-generated fingerprints (default: 0.95)."
	printf "\n"
	echo "# Output options:"
	echo "  -k  | --keep_files                       Keep temporary files [0: Minimum, 1: Salmon output, or 2: All] (default: 0)."
	echo "  -v  | --verbose                          Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2)."
	printf "\n"
}

### Variables

# Colors
red=$(tput setaf 1)
yellow=$(tput setaf 220)
blue=$(tput setaf 27)
green=$(tput setaf 10)
normal=$(tput sgr0)

# Default variables
THREADS=1
READ_TYPE='paired'
COPY_NUMBER_NORMALIZATION='SyFi'
MINSCOREFRACTION=0.95
KEEPF=0
VERBOSE=2

### Parameters

# display whether the obligatory files are called
if [[ $# -lt 2  && $1 != "-h" && $1 != "--help" && $1 != "--citation" && $1 != "--folder_structure" && $1 != "--marker_copy_number_format" ]]; then
	printf "\n${red}ERROR:${normal} You must provide at least the read and fingerprint folders.\n"
	usage
	exit 2
fi

#Get parameters
while [[ "$1" > 0 ]]; do
	case $1 in
		-i | --read_folder)
			shift
			READ_FOLDER=$1
			shift
			;;
		-f | --fingerprint_folder)
			shift
			FINGERPRINT_FOLDER=$1
			shift
			;;
		-r | --read_type)
			shift
			READ_TYPE=$1
			shift
			;;
		-m | --minscorefraction)
			shift
			MINSCOREFRACTION=$1
			shift
			;;
		-t | --threads)
			shift
			THREADS=$1
			shift
			;;
		-k | --keep_files)
			shift
			KEEPF=$1
			shift
			;;
		-v | --verbose)
			shift
			VERBOSE=$1
			shift
			;;
		-h | --help)
			usage
			exit
			;;
		--citation)
			citation
			exit
			;;
		--folder_structure)
			folder_structure
			exit
			;;
		--marker_copy_number_format)
			marker_copy_number_format
			exit
			;;
		*)
			echo $1
			echo "ERROR: Missing parameter."
			exit
			;;
	esac
done

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
	printf "\n${red}Execution halted by user.${normal}\n"
	exit
}

### Checks

# Check: Software dependencies
software=(blastn bwa-mem2 samtools gzip spades.py seqtk seqkit kallisto Rscript gatk3 plot-bamstats bcftools salmon)
for pckg in ${software[@]}
do 
	type ${pckg} 2> /dev/null 1>&2 
	if [ $? != 0 ]
	then
		printf "\n${red}ERROR:${normal} ${pckg} missing. Install or activate SyFi conda environment.\n\n"
		exit
fi
done

# Check: read folder
if [[ ! -d ${READ_FOLDER} ]]; then
	printf "\n${red}ERROR:${normal} Folder ${READ_FOLDER} missing.\n\n"
	exit
fi

# Check: fingerprint folder
if [[ ! -d ${FINGERPRINT_FOLDER} ]]; then
	printf "\n${red}ERROR:${normal} Folder ${FINGERPRINT_FOLDER} missing.\n\n"
	exit
fi

# Check: Keep file
if [[ ! "$KEEPF" =~ ^[0-9]+$ || ${KEEPF} -gt 2 || ${KEEPF} -lt 0 ]]; then
	printf "\n${red}ERROR:${normal} Keep temporary files should be a number between 0 and 2 [0: None, 1: Salmon output, or 2: All].\n\n"
	exit
fi

# Check: Force value
#if [[ ! "$FORCE" =~ ^[0-9]+$ || ${FORCE} -gt 3 || ${FORCE} -lt 0 ]]; then
#  printf "\n${red}ERROR:${normal} Force value should be a number between 0 and 3 [0: None, 1: All, 2: Skipped or 3: Failed].\n\n"
#  exit
#fi

# Check: progress table
function checkProgress() {
	# Colors
	r=$(tput setaf 1)
	y=$(tput setaf 220)
	g=$(tput setaf 10)
	n=$(tput sgr0)

	# Variables
	READ_FOLDER=$1

	# Check file
	if [ -e progress_part_2.txt ]; then
		total=$(ls $READ_FOLDER | wc -l)
		success=$(grep "Success" progress_part_2.txt | wc -l)
		skipped=$(grep "Skipped" progress_part_2.txt | wc -l)
		failed=$(grep "Failed" progress_part_2.txt | wc -l)
		printf "\tSuccess: ${g} ${success} ${n}\n\tSkipped: ${y} ${skipped} ${n}\n\tFailed: ${r}  ${failed} ${n}\n\tTotal:    ${total}\n"
	else
		total=$(ls $READ_FOLDER | wc -l)
		printf "\tSuccess: ${g} 0 ${n}\n\tSkipped: ${y} 0 ${n}\n\tFailed: ${r}  0 ${n}\n\tTotal:    ${total}\n"
	fi

}

### Functions

## FUNCTION: pseudoalignment
function pseudoalignment() {
	# Variables
	read_folder=$1
	sample=$2
	fingerprints=$3
	readtype=$4
	minscorefraction=$5
	threads=$6

	if [ ${VERBOSE} -eq 2 ]; then printf "\n"; fi
	if [ ${VERBOSE} -ge 1 ]; then echo ${sample}; fi
	if [ ${VERBOSE} -eq 2 ]; then printf "Performing pseudoalignment of reads to fingerprints"; fi

	mkdir -p 80-Pseudoalignment/${sample}

	if [ ${readtype} = "paired" ]; then
		 salmon quant -i ${fingerprints} -l IU -1 ${read_folder}/${sample}/${sample}_R1.fastq.gz -2 ${read_folder}/${sample}/${sample}_R2.fastq.gz --validateMappings -o 80-Pseudoalignment/${sample} --minScoreFraction ${minscorefraction} --writeMappings 80-Pseudoalignment/${sample}/salmon_alignments.sam --writeUnmappedNames -p ${threads} &> /dev/null

	elif [ ${readtype} = "single" ]; then
		 salmon quant -i ${fingerprints} -l U -1 ${read_folder}/${sample}/${sample}.fastq.gz --validateMappings -o 80-Pseudoalignment/${sample} --minScoreFraction ${minscorefraction} --writeMappings 80-Pseudoalignment/${sample}/salmon_alignments.sam -p ${threads} &> /dev/null
	fi

	if [ ${VERBOSE} -eq 2 ]; then printf "; Done\n"; fi

}

## FUNCTION: normalizing the count table for marker copy numbers
normalize_isolate_counts() {
	local count_table="$1"
	local copy_number_file="$2"
	local output_file="$3"

	# Use awk to perform the normalization
	awk '
	BEGIN {
		FS = "\t"
		OFS = "\t"
	}

	# Read copy numbers into an array
	NR == FNR {
		if (FNR > 1) { # Skip header line
				copy_numbers[$1] = $2
		}
		next
	}

	# Process the count table
	NR == 1 {
		header = $0
		print header
		next
	}

	{
	asv = $1
	if (asv in copy_numbers) {
		copy_number = copy_numbers[asv]
		printf "%s", asv
		for (i = 2; i <= NF; i++) {
			if (copy_number == 0) {
				printf "\tNaN"
			} else {
				printf "\t%.2f", $i / copy_number
			}
		}
		print ""
	} else {
		printf "%s", asv
		for (i = 2; i <= NF; i++) {
			printf "\tNaN"
		}
		print ""
	}
	}
	' "$copy_number_file" "$count_table" > "$output_file"
}

### Execution

#---------------------------- #
# I. Preparation of the files #
#---------------------------- #

# Call logo
logo ${READ_FOLDER} ${FINGERPRINT_FOLDER} ${READ_TYPE} ${MINSCOREFRACTION} ${THREADS} ${KEEPF} ${VERBOSE}

#Make a new directory
mkdir -p 80-Pseudoalignment 90-Output

#----------------------------- #
# Fingerprint index generation #
#----------------------------- #

if [ ${VERBOSE} -eq 2 ]; then printf "\nConcatenating fingerprints\n"; fi

#Concatenate all the fingerprints from the fingerprint folder
for subf in $(ls ${FINGERPRINT_FOLDER}); do

	## CHECK: Log file
	if [[ -f "01-Logs/quant/log_${subf}.txt" ]]; then
		# Remove old log file
		rm -rf 01-Logs/quant/log_${subf}.txt
		touch 01-Logs/quant/log_${subf}.txt
	else
		# Touch Log File
		mkdir -p 01-Logs/quant
		touch 01-Logs/quant/log_${subf}.txt
	fi

	if [ -f "${fingerprintfolder}/${subf}/${subf}_all_haplotypes.fasta" ]; then
		cat ${fingerprintfolder}/${subf}/${subf}_all_haplotypes.fasta >> 90-Output/SyFi_Fingerprints.fasta
		printf "\n" >> 90-Output/SyFi_Fingerprints.fasta
	else
		printf "Warning: ${fingerprintfolder}/${subf}/${subf}_all_haplotypes.fasta is missing.\n" >>  01-Logs/quant/log_${subf}.txt
	fi
done

#Remove empty lines, place headers to isolate names
sed -i '/^$/d' 90-Output/SyFi_Fingerprints.fasta
sed -i 's/_all_haplotypes//' 90-Output/SyFi_Fingerprints.fasta
sed -i 's/_ALL_HAPLOTYPES//' 90-Output/SyFi_Fingerprints.fasta

#Build Salmon index
salmon index -t 90-Output/SyFi_Fingerprints.fasta -i Fingerprint_index -k 31 &>> 01-Logs/quant/log_salmon.txt
Fingerprint_index="Fingerprint_index"

#------------------------ #
# Copy Number preparation #
#------------------------ #

# SyFi copy number normalized
cat Summary.tsv | cut -f1,11 | grep -v "#" | awk '{ $2 = int($2 + 0.5); print }' | sed '1i\isolate\tcopy_number' | sed 's/ /\t/g' >> 90-Output/copy_number.tsv

#-----------------------------------------------------------#
# Pseudoalignment of samples including confidence intervals #
#---------------------------------------------------------- #

#Loop through the samples in the read folder and convert the pseudolignment sam output to confidence values 
for subf in $(ls ${READ_FOLDER}) ; do
	pseudoalignment ${READ_FOLDER} ${subf} ${Fingerprint_index} ${READ_TYPE} ${MINSCOREFRACTION} ${THREADS}

	awk '
  	  $0 !~ /^@/ {
    	    read_id = $1
	    tx = $3
	    match($0, /AS:i:([0-9]+)/, m)
	    if (m[1] != "") {
	      score = m[1]
	      print read_id, tx, score
	    }
	  }
	' 80-Pseudoalignment/${sample}/salmon_alignments.sam > 80-Pseudoalignment/${sample}/alignments.tsv

	sort -k1,1 80-Pseudoalignment/${sample}/alignments.tsv > 80-Pseudoalignment/${sample}/sorted_alignments.tsv

	SCORE_EXP="1.0"

	awk -v score_exp=${SCORE_EXP} '
	  {
	    if ($1 != prev_read && NR > 1) {
	      # Process all stored lines for previous read
	      for (i = 1; i <= n; i++) {
	        score = scores[i]
	        p = exp(-score_exp * (best - score))
	        tx_sum[txs[i]] += p
	        tx_count[txs[i]] += 1
	      }
	      n = 0
	    }
	    prev_read = $1
	    txs[++n] = $2
	    scores[n] = $3
	    if (n == 1 || $3 > best) {
	      best = $3
	    }
	  }
	  END {
	    # Final read group
	    for (i = 1; i <= n; i++) {
	      score = scores[i]
	      p = exp(-score_exp * (best - score))
	      tx_sum[txs[i]] += p
	      tx_count[txs[i]] += 1
	    }

	    # Output average P(aln) per transcript
	    printf "%-20s %-10s %-10s %-10s\n", "Transcript", "Count", "Sum_P(aln)", "Avg_P(aln)"
	    for (tx in tx_sum) {
	      avg = tx_sum[tx] / tx_count[tx]
	      printf "%-20s %-10d %-10.4f %-10.4f\n", tx, tx_count[tx], tx_sum[tx], avg
	    }
	  }
	' 80-Pseudoalignment/${sample}/sorted_alignments.tsv >> 80-Pseudoalignment/${sample}/confidence_output.tsv

	rm 80-Pseudoalignment/${sample}/alignments.tsv
	rm 80-Pseudoalignment/${sample}/sorted_alignments.tsv
	rm 80-Pseudoalignment/${sample}/salmon_alignments.sam
done

#------------------------------------------------- #
# Assemble the Salmon output into an isolate table #
#------------------------------------------------- #

if [ ${VERBOSE} -eq 2 ]; then printf "\nCollect pseudoalignment results\n"; fi

#Retrieve the pseudoaligment hits from Salmon
for subf in $(ls ${READ_FOLDER}); do
	cut -f 1,5 80-Pseudoalignment/${subf}/quant.sf > 80-Pseudoalignment/${subf}/output.txt
done

# Concatenate the single Salmon outputs into a isolate count table
# Create a variable with the list of files to process
files=$(ls 80-Pseudoalignment/*/output.txt)

# Run paste.awk on these files
awk -f src/paste.awk ${files} > 90-Output/raw_output_table.txt

#Remove header from Salmon output 
sed -i '1d' 90-Output/raw_output_table.txt

#Add first line with the correct sample column headers
samples_header="sample_id"

# Append the file names from the READ_FOLDER, separated by tabs
samples_header+="$(ls ${READ_FOLDER} | xargs -I{} printf "\t%s" {})"

#Insert it into the microbiome table
sed -i "1i $samples_header" 90-Output/raw_output_table.txt

#--------------------------------- #
# Normalize for marker copy number #
#--------------------------------- #

if [ ${VERBOSE} -eq 2 ]; then printf "Normalize pseudoalignment results\n"; fi

# SyFi normalized
normalize_isolate_counts "90-Output/raw_output_table.txt" "90-Output/copy_number.tsv" "90-Output/norm_output_table.txt"
sed -i '1d' 90-Output/norm_output_table.txt
sed -i "1i $samples_header" 90-Output/norm_output_table.txt

#-------------------------- #
# Collect Confidence values #
#-------------------------- #

if [ ${VERBOSE} -eq 2 ]; then printf "\nCollect confidence values of pseudoalignments\n"; fi

#Retrieve the confidence values of pseudoaligment hits from Salmon
for subf in $(ls ${READ_FOLDER}); do
        sed 's/  */\t/g' 80-Pseudoalignment/${subf}/confidence_output.tsv | cut -f 1,4 > 80-Pseudoalignment/${subf}/intermediate_conf_output.txt
done

# Collate the Salmon confidence values of individual pseudoalignment runs
# Create a variable with the list of files to process
files=$(ls 80-Pseudoalignment/*/intermediate_conf_output.txt)

# Run paste.awk on these files
awk -f src/paste.awk ${files} > 90-Output/confidence_output_table.txt

#Remove header from Salmon output
sed -i '1d' 90-Output/confidence_output_table.txt

#Add first line with the correct sample column headers
samples_header="sample_id"

# Append the file names from the READ_FOLDER, separated by tabs
samples_header+="$(ls ${READ_FOLDER} | xargs -I{} printf "\t%s" {})"

#Insert it into the microbiome table
sed -i "1i $samples_header" 90-Output/confidence_output_table.txt

#----------------------- #
# Collect unmapped reads #
#----------------------- #

printf "Sample\tNo_mapped_reads\tNo_unmapped_reads\tproportion_mapped\n" > 90-Output/unmapped_reads.txt

for subf in $(ls ${READ_FOLDER}); do
    mapped=$(cut -f5 80-Pseudoalignment/${subf}/quant.sf | grep -v "NumReads" | awk '{sum += $1} END {print sum}')
    unmapped=$(wc -l < 80-Pseudoalignment/${subf}/aux_info/unmapped_names.txt)
    total=$(awk -v m=$mapped -v u=$unmapped 'BEGIN {print m + u}')
    percentage=$(awk -v m=$mapped -v t=$total 'BEGIN {if (t > 0) printf "%.6f", m / t; else print "NaN"}')
    printf "%s\t%d\t%d\t%s\n" "$subf" "$mapped" "$unmapped" "$percentage" >> 90-Output/unmapped_reads.txt
done

#--------- #
# Clean up #
#--------- #

if [ ${KEEPF} -eq 0 ]; then
	rm -r Fingerprint_index &> /dev/null
	rm -r 80-Pseudoalignment/*/aux_info 80-Pseudoalignment/*/cmd_info.json  80-Pseudoalignment/*/lib_format_counts.json  80-Pseudoalignment/*/libParams  80-Pseudoalignment/*/logs &> /dev/null
fi

if [ ${KEEPF} -eq 1 ]; then
	rm -r Fingerprint_index &> /dev/null
fi

# END
if [ ${VERBOSE} -eq 2 ]; then printf "\nSyFi mic-drop\n"; fi
