#!/bin/bash

# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Gijs Selten, Florian Lamouche and Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk / g.selten@uu.nl
# Created Date  : 
# version       : '1.0'
# ---------------------------------------------------------------------------
# Pipeline (bash) to retrieve amplicon fingerprints from gene fingerprints using in silico primers
# ---------------------------------------------------------------------------

function logo() {
	# Colors
	bold=$(tput bold)
	normal=$(tput sgr0)

	# Arrays
	array_keep=("None" "BAMs" "All")
	array_verbose=("Quiet" "Sample" "All")
	array_force=("None" "All" "Skipped" "Failed")

	# Variables
	if [[ $1 != "help" ]]; then
		folder=$1
		fpprimer=$2
		rpprimer=$3
		minlength=$4
		maxlength=$5
		keepfiles=${array_keep[${6}]}
		verbose=${array_keep[${7}]}
		force=${array_force[${8}]}

		# Logo
		echo ""
		echo "${bold}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${normal}"
		echo "${bold}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${normal}"
		echo "${bold}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${normal}"
		echo "${bold}|____|/___/    /____/   |___|¯¯    |________|${normal}"
		echo ""
		echo "${bold}----------------------------------------------${normal}"
		echo "${bold}Summary:${normal}"
		echo "    Folder: ${folder}"
		echo "    Forward primer: ${fpprimer}"
		echo "    Reverse primer: ${rpprimer}"
		echo "    Minimum length: ${minlength}"
		echo "    Maximum length: ${maxlength}"
		echo "    Keep files: ${keepfiles}"
		echo "    Verbose: ${verbose}"
		echo "    Force: ${force}"
		echo ""
		printf "${bold}Start:${normal}\n"
		date "+    date: %d/%m/%Y"
		date "+    time: %H:%M:%S"
		echo ""
		echo "${bold}Progress:${normal}"
		checkProgress ${folder}
		echo "${bold}----------------------------------------------${normal}"
		echo ""
	elif [[ $8 != 0 ]]; then
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
	echo "Usage: ./$0 -i <INPUT_FOLDER> -f <FORWARD_PRIMER> -r <REVERSE_PRIMER> -n <MINIMUM_LENGTH> -x <MAXIMUM_LENGTH>"
	printf "\n"
	echo "${bold}REQUIRED:${normal}"
	echo "# Input"
	echo "  -i   | --input_folder      The 70-Fingerprint folder containing SyFi output executed on full length target"
	echo "  -fp  | --forward_primer    Forward primer for extracting amplicon of interest"
	echo "  -rp  | --reverse_primer    Reverse primer for extracting amplicon of interest"
	echo "  -n   | --minimum_length    Minimum length of the amplicon of interest"
	echo "  -x   | --maximum_length    Maximum length of the amplicon of interest" 
	printf "\n"
	echo "${bold}OPTIONAL:${normal}"
	echo "# Output options:"
	echo "  -k  | --keep_files       Keep temporary files [0: Minimum, 1: BAM's, or 2: All] (default: 0)."
	echo "  -v  | --verbose          Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2)."
	echo "  -f  | --force            Force re-computation of computed samples [0: None, 1: All, 2: Skipped, or 3: Failed] (default: 0)."
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
FORCE=0
KEEPF=0
VERBOSE=2

### Parameters

# display usage if 
if [[ $# -lt 5  && $1 != "-h" && $1 != "--help" && $1 != "--citation" ]]; then
	printf "\n${red}ERROR:${normal} You must provided at least the SyFi output folder, the primer sequences, and the minimum and maximum length of the amplicon.\n"
	usage
	exit 2
fi

#Get parameters
while [[ "$1" > 0 ]]; do
	case $1 in
		-i | --input_folder)
			shift
			INPUT_FOLDER=$1
			shift
			;;
		-fp | --forward_primer)
			shift
			FORWARD_PRIMER=$1
			shift
			;;
		-rp | --reverse_primer)
			shift
			REVERSE_PRIMER=$1
			shift
			;;
		-n | --minimum_length)
			shift
			MINIMUM_LENGTH=$1
			shift
			;;
		-x | --maximum_length)
			shift
			MAXIMUM_LENGTH=$1
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
		-f | --force)
			shift
			FORCE=$1
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
	printf "\n${red}Execution halted by user.${normal}\n" >> 01-Logs/log_${subf}.txt
	exit
}

### Checks

# Check: Software dependencies
software=(qiime)
for pckg in ${software[@]}
do 
	type ${pckg} 2> /dev/null 1>&2 
	if [ $? != 0 ]
	then
		printf "\n${red}ERROR:${normal} ${pckg} missing. Install or activate SyFi conda environment.\n\n"
		exit
fi
done

# Check: Input folder
if [[ ! -d ${INPUT_FOLDER} ]]; then
	printf "\n${red}ERROR:${normal} Folder ${INPUT_FOLDER} missing.\n\n"
	exit
fi

#Check forward primer
if [[ -z ${FORWARD_PRIMER} ]]; then
	printf "\n${red}ERROR:${normal} Forward primer missing.\n\n"
	exit 
fi

#Check reverse primer
if [[ -z ${REVERSE_PRIMER} ]]; then
	printf "\n${red}ERROR:${normal} Reverse primer missing.\n\n"
	exit
fi

#Check minimum length of amplicon sequence
if [[ -z ${MINIMUM_LENGTH} ]]; then
	printf "\n${red}ERROR:${normal} Minimum length of the amplicon sequence is missing.\n\n"
	exit
fi

#Check reverse primer
if [[ -z ${MAXIMUM_LENGTH} ]]; then
	printf "\n${red}ERROR:${normal} Maximum length of the amplicon sequence is missing.\n\n"
	exit
fi

# Check: Keep file
if [[ ! "$KEEPF" =~ ^[0-9]+$ || ${KEEPF} -gt 2 || ${KEEPF} -lt 0 ]]; then
	printf "\n${red}ERROR:${normal} Keep temporary files should be a number between 0 and 2 [0: None, 1: BAM's, or 2: All].\n\n"
	exit
fi

# Check: Force value
if [[ ! "$FORCE" =~ ^[0-9]+$ || ${FORCE} -gt 3 || ${FORCE} -lt 0 ]]; then
	printf "\n${red}ERROR:${normal} Force value should be a number between 0 and 3 [0: None, 1: All, 2: Skipped or 3: Failed].\n\n"
	exit
fi

# Check: progress table
function checkProgress() {
	# Colors
	r=$(tput setaf 1)
	y=$(tput setaf 220)
	g=$(tput setaf 10)
	n=$(tput sgr0)

	# Variables
	INPUT_FOLDER=$1
}

### Functions

# Function: Retrieving amplicon from gene marker
function amplicon() {
	#THIS SHOULD BE REPLACED BY SYFI ALREADY HAVING QIIME2 INSTALLED
	#eval "$(conda shell.bash hook)" && conda activate qiime2-2022.8

	# Variables
	mode=$1 # unique/multiple
	subf=$2
	fpprimer=$3
	rpprimer=$4
	minlength=$5
	maxlength=$6

	# Create folder
	mkdir -p 71-Amplicon/${subf}

	if [ $mode == "unique" ];then
	 # Log
	 if [ ${VERBOSE} -eq 2 ]; then printf "Retrieving amplicon from single haplotype\n"; fi
	 printf "\n\n### Amplicon ###\n\n" >> 01-Logs/log_${subf}.txt

	 # Set input sequence for amplicon retrieval
	 GENPATH=70-Fingerprints/${subf}
	 QZAPATH=71-Amplicon/${subf}
	 OUTPATH=71-Amplicon/${subf}

	 for file in "$GENPATH"/*.fasta
	 do
		sed -i 's/[a-z]/\U&/g' $file
		FNAME=${file%%.fasta*}
		FNAME=$( basename $FNAME )
		qiime tools import \
			 --input-path $file \
			 --output-path $QZAPATH/${FNAME}.qza \
			 --type 'FeatureData[Sequence]'
		qiime feature-classifier extract-reads \
			 --i-sequences $QZAPATH/${FNAME}.qza \
			 --p-f-primer ${fpprimer} \
			 --p-r-primer ${rpprimer} \
			 --p-min-length ${minlength} \
			 --p-max-length ${maxlength} \
			 --o-reads $QZAPATH/${FNAME}_amp.qza
		qiime tools export \
			 --input-path $QZAPATH/${FNAME}_amp.qza \
			 --output-path $OUTPATH/${FNAME}_amp.fna

		#Clean up the folder
		mv $OUTPATH/${FNAME}_amp.fna/* $OUTPATH/${FNAME}.fasta
		rm -r $OUTPATH/${FNAME}_amp.fna/ $QZAPATH/${FNAME}_amp.qza $QZAPATH/${FNAME}.qza

		done
	#This should also be removed later
		#conda deactivate

	elif [ $mode == "multiple" ];then
		# Log
		if [ ${VERBOSE} -eq 2 ]; then printf "Retrieving amplicons from multiple haplotypes\n"; fi
		printf "\n\n### Amplicon ###\n\n" >> 01-Logs/log_${subf}.txt

		# Set input sequence for amplicon retrieval
	 GENPATH=70-Fingerprints/${subf}
	 QZAPATH=71-Amplicon/${subf}
	 OUTPATH=71-Amplicon/${subf}

	 for file in "$GENPATH"/seq_h*.fasta
	 do
		sed -i 's/[a-z]/\U&/g' $file
		FNAME=${file%%.fasta*}
		FNAME=$( basename $FNAME )
		qiime tools import \
			 --input-path $file \
			 --output-path $QZAPATH/${FNAME}.qza \
			 --type 'FeatureData[Sequence]'
		qiime feature-classifier extract-reads \
			 --i-sequences $QZAPATH/${FNAME}.qza \
			 --p-f-primer ${fpprimer} \
			 --p-r-primer ${rpprimer} \
			 --p-min-length ${minlength} \
			 --p-max-length ${maxlength} \
			 --o-reads $QZAPATH/${FNAME}_amp.qza
		qiime tools export \
			 --input-path $QZAPATH/${FNAME}_amp.qza \
			 --output-path $OUTPATH/${FNAME}_amp.fna

		#Clean up the folder
		mv $OUTPATH/${FNAME}_amp.fna/* $OUTPATH/${FNAME}.fasta
		rm -r $OUTPATH/${FNAME}_amp.fna/ $QZAPATH/${FNAME}_amp.qza $QZAPATH/${FNAME}.qza

		done
	#This should also be removed later
		#conda deactivate
	fi

	# Keep haplotypes filtered in integration
	keep=($(tail -n +2 60-Integration/${subf}/integration.tsv | cut -f 1))
	for i in $(ls 71-Amplicon/${subf} | sed 's/.fasta//g'); do
		if [[ ! "${keep[*]}" =~ "${i}" && ${mode} == "multiple" ]]; then
			rm -f 71-Amplicon/${subf}/${i}.fasta
		fi
	done

 # Concatenate haplotypes (NF-1 to obtain "per_haplotype")
		cat 60-Integration/${subf}/integration.tsv | awk '{print $1 "/" $(NF-1)}' | tail -n +2 > 71-Amplicon/${subf}/tmp.tsv
		for h in $(cat 71-Amplicon/${subf}/tmp.tsv); do
			seq=$(echo ${h} | cut -d "/" -f 1)
			num=$(echo ${h} | cut -d "/" -f 2)
			for n in $(seq ${num}); do
				# End when haplotype ratio > 1
				if [[ ${n} == ${num} && ${n} != 1 ]]; then
					grep -v "^>" 71-Amplicon/${subf}/${seq}.fasta >> 71-Amplicon/${subf}/${subf}_all_haplotypes.tmp.fasta
				# End when haplotype ratio == 1
				elif [[ ${n} == ${num} && ${n} == 1 ]]; then
					grep -v "^>" 71-Amplicon/${subf}/${seq}.fasta >> 71-Amplicon/${subf}/${subf}_all_haplotypes.tmp.fasta
				# Addition of haplotypes before ending
				else
					grep -v "^>" 71-Amplicon/${subf}/${seq}.fasta >> 71-Amplicon/${subf}/${subf}_all_haplotypes.tmp.fasta
					echo "NNNNNNNNNN" >> 71-Amplicon/${subf}/${subf}_all_haplotypes.tmp.fasta
				fi
			done
		done
		printf ">${subf}_all_haplotypes\n" > 71-Amplicon/${subf}/${subf}_all_haplotypes.fasta
		cat 71-Amplicon/${subf}/${subf}_all_haplotypes.tmp.fasta | tr -d "\n" >> 71-Amplicon/${subf}/${subf}_all_haplotypes.fasta

 
}


function CleanFiles() {
	# Variables
	SAMPLE=$1
	KEEPF=$2

	# Array
	keep=("./01-Logs/log_${SAMPLE}.txt" "./10-Blast/${SAMPLE}.tsv" "./11-Sequences/${SAMPLE}/${SAMPLE}.fasta" "./20-Alignment/${SAMPLE}/${SAMPLE}.fasta" "./20-Alignment/${SAMPLE}/${SAMPLE}_R1.fastq.gz" "./20-Alignment/${SAMPLE}/${SAMPLE}_R2.fastq.gz" "./30-VariantCalling/${SAMPLE}/variants/${SAMPLE}.vcf.gz" "./40-Phasing/${SAMPLE}/${SAMPLE}_phased.vcf.gz" "./50-Haplotypes/${SAMPLE}/clean_${SAMPLE}_haplotypes.fasta" "./60-Integration/${SAMPLE}/abundance.tsv" "./60-Integration/${SAMPLE}/copy_number.tsv" "./60-Integration/${SAMPLE}/integration.tsv" "./70-Fingerprints/${SAMPLE}/${SAMPLE}_all_haplotypes.fasta" "./71-Amplicon/${SAMPLE}/${SAMPLE}_all_haplotypes.fasta")

	# Remove files
	if [ ${KEEPF} -eq 0 ]; then
		# Find and remove
		for file in $(find . -type f -name "*${SAMPLE}*"); do
			if [[ ! "${keep[*]}" =~ "${file}" ]]; then
				if [[ $(echo ${file} | grep "71-Amplicon" | wc -l) != 0 ]]; then
					continue
				elif [[ $(echo ${file} | grep -e "seq_h" -e "assembly_h" -e ".fasta" | grep -e "71-Amplicon" | wc -l) != 0 ]]; then
					continue
				else
					rm -rf ${file}
				fi
			fi
		done

		# Remove empty directories
		find . -type d -empty -delete
		# Remove non-captured files
		rm -rf 71-Amplicon/${SAMPLE}/*tmp*

	elif [ ${KEEPF} -eq 1 ]; then
		# Find and remove
		for file in $(find . -type f -name "*${SAMPLE}*"); do
			if [[ ! "${keep[*]}" =~ "${file}" ]]; then
				if [[ $(echo ${file} | grep "71-Amplicon" | wc -l) != 0 ]]; then
					continue
				elif [[ $(echo ${file} | grep -e "seq_h" -e "assembly_h" -e ".fasta" | grep -e "71-Amplicon" | wc -l) != 0 ]]; then
					continue
				elif [[ $(echo ${file} | grep ".sort.bam" | wc -l) != 0 ]]; then
					continue
				else
					rm -rf ${file}
				fi
			fi
		done    

		# Remove empty directories
		find . -type d -empty -delete
		# Remove non-captured files
		rm -rf 71-Amplicon/${SAMPLE}/*tmp*

		#Make sure permissions are correct
		chmod uog+rwx 71-Amplicon/${SAMPLE}
		chmod uog+rwx 71-Amplicon/${SAMPLE}/*
	fi
}

### Execution

#----------------------------------------------------------------------------- #
# I. Retrieving amplicon sequences from haplotypes and build a new fingerprint #
#----------------------------------------------------------------------------- #

# Call logo
logo ${INPUT_FOLDER} ${FORWARD_PRIMER} ${REVERSE_PRIMER} ${MINIMUM_LENGTH} ${MAXIMUM_LENGTH} ${KEEPF} ${VERBOSE} ${FORCE}

mkdir -p 71-Amplicon

for subf in $(ls ${INPUT_FOLDER}); do

	# Continue with old log file
	touch 01-Logs/log_${subf}.txt

	if [ ${VERBOSE} -ge 1 ]; then date "+start time: %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt; fi
	if [ ${VERBOSE} -ge 1 ]; then printf "${blue}Sample:${normal} $subf\n" | tee -a 01-Logs/log_${subf}.txt; fi

	#Check whether there is only one haplotype or not
	if [[ -f 30-VariantCalling/${subf}/variants/${subf}.vcf.gz && $(zcat 30-VariantCalling/${subf}/variants/${subf}.vcf.gz | grep -v "#" | wc -l) == 0 && $(grep "^>" 20-Alignment/${subf}/${subf}.fasta | wc -l) == 1 ]]; then

	# One haplotype
	amplicon 'unique' ${subf} ${FORWARD_PRIMER} ${REVERSE_PRIMER} ${MINIMUM_LENGTH} ${MAXIMUM_LENGTH}

	else
	# Multiple haplotypes
	amplicon 'multiple' ${subf} ${FORWARD_PRIMER} ${REVERSE_PRIMER} ${MINIMUM_LENGTH} ${MAXIMUM_LENGTH}
 fi

 if [ -f 71-Amplicon/${subf}/${subf}_all_haplotypes.fasta ]; then
	 printf "${subf}\tAmplicon retrieved\n" >> progress_amp.txt
 fi

 CleanFiles ${subf} ${KEEPF}

 date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
 echo ""

done

# Final format

if [ ${VERBOSE} -ge 1 ]; then
	printf "Finish:\n"
	date "+    date: %d/%m/%Y" 
	date "+    time: %H:%M:%S"
	printf "\n"
fi
