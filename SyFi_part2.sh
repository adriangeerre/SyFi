#!/bin/bash

# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Gijs Selten, Florian Lamouche and Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk / g.selten@uu.nl
# Created Date  : 06/04/2022 - 07/11/2022
# version       : '1.0'
# ---------------------------------------------------------------------------
# Pipeline (bash) to perform fingerprint identification from microbiome data.
# ---------------------------------------------------------------------------

function logo() {
  # Colors
  w=$(tput setaf 231)
  n=$(tput sgr0)

  # Arrays
  array_keep=("None" "Salmon output" "All")
  array_verbose=("Quiet" "Sample" "All")

  # Variables
  if [[ $1 != "help" ]]; then
    readfolder=$1
    fingerprintfolder=$2
    readtype=$3
    copynumbernormalization=$4
    copynumber=$5
    minscorefraction=$6
    threads=$7
    keepfiles=${array_keep[${8}]}
    verbose=${array_keep[${9}]}

    # Logo
    echo ""
    echo "${w}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${n}"
    echo "${w}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${n}"
    echo "${w}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${n}"
    echo "${w}|____|/___/    /____/   |___|¯¯    |________|${n}"
    echo ""
    echo "----------------------------------------------"
    echo "Summary:"
    echo "    Read folder: ${readfolder}"
    echo "    Fingerprint folder: ${fingerprintfolder}"
    echo "    Read type: ${readtype}"
    echo "    Copy number normalization: ${copynumbernormalization}"
    echo "    Copy number file: ${copynumber}"
    echo "    Minimum pseudoalignment identity score: ${minscorefraction}"
    echo "    Threads: ${threads}"
    echo "    Keep files: ${keepfiles}"
    echo "    Verbose: ${verbose}"
    echo ""
    printf "Start:\n"
    date "+    date: %d/%m/%Y"
    date "+    time: %H:%M:%S"
    echo ""
    echo "Progress:"
    checkProgress ${readfolder}
    echo "----------------------------------------------"
    echo ""
  else
    # Logo
    echo ""
    echo "${w}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${n}"
    echo "${w}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${n}"
    echo "${w}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${n}"
    echo "${w}|____|/___/    /____/   |___|¯¯    |________|${n}"
    echo ""
  fi
}

function usage()
{
  # Colors
  w=$(tput setaf 231)
  n=$(tput sgr0)

  # Logo
  logo help

  # Usage
  echo "${w}Usage: ./$0 -i <READ_FOLDER> -f <FINGERPRINT_FOLDER> -t <THREADS>${n}"
  printf "\n"
  echo "${w}REQUIRED:${n}"
  echo "# Input"
  echo "  -i  | --read_folder                      Folder containing the sample reads. The software assumes that the folder contains sub-folders for each sample containing the fastq.gz files of either single end reads or paired end reads. For more details, execute <pipeline --folder_structure>."
  echo "  -f  | --fingerprint_folder               SyFi part 1 output folder that contains the Fingerprints. If the user wants to exclude bacterial members, a personalized folder can be created with the SynCom members of interest using the same folder structure as the SyFi part 1 output."
  printf "\n"
  echo "${w}OPTIONAL:${n}"
  echo "# Read type"
  echo "  -r  | --read_type                        Paired or single end reads [paired or single] (default: paired)."
  printf "\n"
  echo "# Copy number normalization:"
  echo "  -n  | --copy_number_normalization        Whether the microbiome or isolate count table needs to be normalized for marker copy number, and whether the marker copy number should be taken from the SyFi output or from another source [SyFi, other, or no] (default: SyFi)."
  echo "  -c  | --copy_number_file                 If user specified 'other' for the marker copy number normalization, here the file should be supplied for the marker copy number. The format of this personalized file should be the name of the SynCom member followed by the marker copy number. For details, execute <pipeline --marker_copy_number_format>."
  printf "\n"
  echo "# Pseudoalignment percentage identity score:"
  echo "  -m  | --minscorefraction                 Percentage identity score that Salmon uses for pseudoaligning metagenomic reads to the SyFi-generated fingerprints (default: 0.95)."
  printf "\n"
  echo "# Output options:"
  echo "  -k  | --keep_files                       Keep temporary files [0: Minimum, 1: Salmon output, or 2: All] (default: 0)."
  echo "  -v  | --verbose                          Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2)."
  printf "\n"
  echo "# Display:"
  echo "  -h  | --help                             Display help."
  echo "  --citation                               Display citation."
  echo "  --folder_structure                       Display required folder structure."
  echo "  --marker_copy_number_format              Display required marker copy number file format"
  printf "\n"
}

function folder_structure()
{
  printf "\n"
  echo "The pipeline assumes that the sample reads are organized in sub-folders inside of the input folder (-i | --read_folder). Each sub-folder should contain the reads either in single end or paired end format (with extension .fastq.gz)."
  echo "For example:"
  echo ""
  echo "read_folder/
          └── sample_1
              ├── sample_1_R1.fastq.gz
              └── sample_1_R2.fastq.gz
          └── samplen_2
              ├── sample_2_R1.fastq.gz
              └── sample_2_R2.fastq.gz"
  printf "\n"
}

function marker_copy_number_format()
{
  printf "\n"
  echo "The pipeline takes the marker copy number from the SyFi part 1 output (Summary.tsv) as default (-n SyFi). Personalized marker copy numbers can be supplied as a separate tsv in the following format: "
  echo ""
  echo "isolate	copy_number"
  echo "strain_1	2"
  echo "strain_2	4"
  printf "\n"
}

function citation()
{
  logo help
  printf "\n"
  echo "Thank you for using SyFi. Please, cite:"
  echo ""
  echo "<Include citation>"
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

# Define software path
#SYFI_BASE=$(dirname $(whereis ${0} | cut -d " " -f 2))

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
    -n | --copy_number_normalization)
      shift
      COPY_NUMBER_NORMALIZATION=$1
      shift
      ;;
    -c | --copy_number_file)
      shift
      COPY_NUMBER_FILE=$1
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
  printf "\n${red}Execution halted by user.${normal}\n" >> 10-Pseudoaligment_logs/log_${subf}.txt
  exit
}

### Checks

# Check: Software dependencies
#software=(blastn bwa-mem2 samtools gzip spades.py seqtk seqkit kallisto Rscript gatk plot-bamstats bcftools)
#for pckg in ${software[@]}
#do 
#  type ${pckg} 2> /dev/null 1>&2 
#  if [ $? != 0 ]
 # then
 #   printf "\n${red}ERROR:${normal} ${pckg} missing. Install or activate SyFi conda environment.\n\n"
 #   exit
#fi
#done

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

  printf "\n"
  echo "Performing pseudoalignment of reads to fingerprints"

  mkdir -p 80-Pseudoalignment/${sample}

  if [ ${readtype} = "paired" ]; then
     ~/salmon/salmon-latest_linux_x86_64/bin/salmon quant -i ${fingerprints} -l IU -1 ${read_folder}/${sample}/${sample}_R1.fastq.gz -2 ${read_folder}/${sample}/${sample}_R2.fastq.gz --validateMappings -o 80-Pseudoalignment/${sample} --minScoreFraction ${minscorefraction} -p ${threads} 1> /dev/null

  elif [ ${readtype} = "single" ]; then
     ~/salmon/salmon-latest_linux_x86_64/bin/salmon quant -i ${fingerprints} -l U -1 ${read_folder}/${sample}/${sample}.fastq.gz --validateMappings -o 80-Pseudoalignment/${sample} --minScoreFraction ${minscorefraction} -p ${threads} /dev/null
  fi

  printf "\n"
  echo "Finished pseudoalignment of reads of ${sample}"

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
              printf "\t%.2f", $i / copy_number
          }
          print ""
      } else {
          print "ASV not found in copy numbers:", asv
      }
  }
  ' "$copy_number_file" "$count_table" > "$output_file"
}

### Execution

#---------------------------- #
# I. Preparation of the files #
#---------------------------- #

# Call logo
logo ${READ_FOLDER} ${FINGERPRINT_FOLDER} ${READ_TYPE} ${COPY_NUMBER_NORMALIZATION} ${COPY_NUMBER_FILE} ${MINSCOREFRACTION} ${THREADS} ${KEEPF} ${VERBOSE}

#Make a new directory
mkdir -p 80-Pseudoalignment 

#----------------------------- #
# Fingerprint index generation #
#----------------------------- #

#Concatenate all the fingerprints from the fingerprint folder
printf "\nConcatenating fingerprints\n"
for subf in $(ls ${FINGERPRINT_FOLDER}); do
  cat ${fingerprintfolder}/${subf}/${subf}_all_haplotypes.fasta >> SyFi_Fingerprints.fasta ;
  printf "\n" >> SyFi_Fingerprints.fasta  ;
done

#Remove empty lines, place headers to isolate names
sed -i '/^$/d' SyFi_Fingerprints.fasta
sed -i 's/_all_haplotypes//' SyFi_Fingerprints.fasta
sed -i 's/_ALL_HAPLOTYPES//' SyFi_Fingerprints.fasta


#Build Salmon index
~/salmon/salmon-latest_linux_x86_64/bin/salmon index -t SyFi_Fingerprints.fasta -i Fingerprint_index -k 31 1> /dev/null
Fingerprint_index="Fingerprint_index"

#------------------------ #
# Copy Number preparation #
#------------------------ #

#If statement that default copy number file is taken if not specified otherwise. New file created for copy number normalization

if [ ${copynumbernormalization} == "SyFi" ]; then
  cat Summary.tsv | cut -f1,11 | grep -v "#" | awk '{ $2 = int($2 + 0.5); print }' | sed '1i\isolate\tcopy_number' | sed 's/ /\t/g'  >> copy_number.tsv
  elif [ ${copynumbernormalization} == "other" ] ; then
    mv ${copynumber} copy_number.tsv
  elif [ ${copynumbernormalization} == "no" ] ; then
    ls ${FINGERPRINT_FOLDER} | sed 's/$/\t1/g' |  sed '1i\isolate\tcopy_number' >> copy_number.tsv
  fi

#--------------------------- #
# Pseudoalignment of samples #
#--------------------------- #

#Loop through the samples in the read folder

for subf in $(ls ${READ_FOLDER}) ; do
  pseudoalignment ${READ_FOLDER} ${subf} ${Fingerprint_index} ${READ_TYPE} ${MINSCOREFRACTION} ${THREADS}
done

#------------------------------------------------- #
# Assemble the Salmon output into an isolate table #
#------------------------------------------------- #

#Retrieve the pseudoaligment hits from Salmon
for subf in $(ls ${READ_FOLDER}); do
  cut -f 1,5 80-Pseudoalignment/${subf}/quant.sf >> 80-Pseudoalignment/${subf}/output.txt
done

# Concatenate the single Salmon outputs into a isolate count table
# Create a variable with the list of files to process
files=$(ls 80-Pseudoalignment/*/output.txt)

# Run paste.awk on these files
awk -f src/paste.awk ${files} > raw_output_table.txt

#Remove header from Salmon output 
sed -i '1d' raw_output_table.txt

#Add first line with the correct sample column headers
samples_header="sample_id"

# Append the file names from the READ_FOLDER, separated by tabs
samples_header+="$(ls ${READ_FOLDER} | xargs -I{} printf "\t%s" {})"

#Insert it into the microbiome table
sed -i "1i $samples_header" raw_output_table.txt

#--------------------------------- #
# Normalize for marker copy number #
#--------------------------------- #

normalize_isolate_counts "raw_output_table.txt" "copy_number.tsv" "norm_output_table.txt"
sed -i '1d' norm_output_table.txt
sed -i "1i $samples_header" norm_output_table.txt

#--------- #
# Clean up #
#--------- #

mkdir -p 90-Output
mv SyFi_Fingerprints.fasta 90-Output
mv raw_output_table.txt 90-Output
mv norm_output_table.txt 90-Output
rm -r Fingerprint_index
rm -r 80-Pseudoalignment/*/aux_info  80-Pseudoalignment/*/cmd_info.json  80-Pseudoalignment/*/lib_format_counts.json  80-Pseudoalignment/*/libParams  80-Pseudoalignment/*/logs
rm 0


printf "\nSyFi mic-drop"
