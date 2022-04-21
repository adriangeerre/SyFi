#!/bin/bash

# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Gijs Selten, Florian Lamouche and Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 06/04/2022
# version       : '1.0'
# ---------------------------------------------------------------------------
# Pipeline (bash) to perform fingerprint identification from microbiome data.
# ---------------------------------------------------------------------------

function usage()
{
    echo "Usage: $0 -i <INPUT_FOLDER> -s <SEARCH_TARGET> -p <PREFIX> -t <THREADS>"
    printf "\n"
    echo "  -i  | --input_folder    Folder containing input genomes and reads. The software assumes that the folder contains sub-folders for each strain. For more detail, execute <pipeline --folder_structure> (REQUIRED)."
    echo "  -s  | --search_target   Genomic region of interest in fasta format, e.g., 16S (REQUIRED)."
    echo "  -p  | --prefix          Prefix for output files (default: project)."
    echo "  -t  | --threads         Number of threads (default: 1)."
    echo "  -h  | --help            Display help."
    echo "  -c  | --citation        Display citation."
    printf "\n"
}

function folder_structure()
{
  printf "\n"
  echo "The pipeline assumes that the genomes and reads are organized in sub-folders inside of the input folder (-i | --input_folder). Each sub-folder should contain the genome (.fasta) and the reads (.fastq.gz)."
  echo "For example:"
  echo ""
  echo "input_folder/
          └── sub-folder_1
              ├── strain_1_R1.fastq.gz
              ├── strain_1_R2.fastq.gz
              └── strain_1_genome.fasta
          └── sub-folder_2
              ├── strain_2_R1.fastq.gz
              ├── strain_2_R2.fastq.gz
              └── strain_2_genome.fasta"
  printf "\n"
}

function citation()
{
  printf "\n"
  echo "Thank you for using <CHANGE_NAME>. Please, cite:"
  echo ""
  echo "<Include citation>"
  printf "\n"
}

# Default variables
PREFIX='project'
THREADS=1

# display usage if 
if [[ $# -lt 4  && $1 != "-h" && $1 != "--help" && $1 != "-c" && $1 != "--citation" && $1 != "--folder_structure" ]]; then
  echo "ERROR: You must provided at least the assembly and alignment files."
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
    -s | --search_target)
      shift
      SEARCH_TARGET=$1
      shift
      ;;
    -p | --prefix)
      shift
      PREFIX=$1
      shift
      ;;
    -t | --threads)
      shift
      THREADS=$1
      shift
      ;;  
    -h | --help)
      usage
      exit
      ;;
    -c | --citation)
      citation
      exit
      ;;
    --folder_structure)
      folder_structure
      ;;
    *)
      echo "ERROR: Missing parameter."
      exit
      ;;
  esac
done


# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
  printf "\nExecution halted by user.\n"
  exit
}

### Checks

# Input folder
if [[ ! -d ${INPUT_FOLDER} ]]; then
  echo "ERROR: Folder ${INPUT_FOLDER} missing."
  exit
fi

### Execution

## Mapping (BLAST)
mkdir -p 10-Blast 11-Sequences
for subf in $(ls ${INPUT_FOLDER}); do
  mkdir -p 11-Sequences/${subf}
  blastn -query ${INPUT_FOLDER}/${subf}/*.fasta -subject ${target} -strand both -outfmt "6 std qseq" > 10-Blast/${subf}.tsv
  printf "> ${subf}\n$(cat 10-Blast/CHA0_modified.tsv | head -n 1 | cut -f 13 | sed 's/-//g')" > 11-Sequences/${subf}/${subf}.fasta
done 

# List of strains
# ls 11-Sequences/*.fasta | cut -d "." -f 1 > list.txt

## Alignment, SamToBam & BamToFastq (BWA + Samtools)
mkdir -p 20-Alignment
for subf in $(ls ${INPUT_FOLDER}); do
  mkdir -p 20-Alignment/${subf}
  bwa-mem2 index 11-Sequences/${subf}/${subf}.fasta
  bwa-mem2 mem 11-Sequences/${subf}/${subf}.fasta ${INPUT_FOLDER}/${subf}/${subf}_R1.fastq.gz ${INPUT_FOLDER}/${subf}/${subf}_R2.fastq.gz -t ${THREADS} > 20-Alignment/${subf}/${subf}.sam
  samtools view -bS 20-Alignment/${subf}/${subf}.sam -@ ${THREADS} > 20-Alignment/${subf}/${subf}.bam
  samtools collate 20-Alignment/${subf}/${subf}.bam 20-Alignment/${subf}/${subf}.collate
  samtools fastq -1 20-Alignment/${subf}/${subf}_16S_R1.fastq -2 20-Alignment/${subf}/${subf}_16S_R2.fastq -s 20-Alignment/${subf}/${subf}_leftover.fastq 20-Alignment/${subf}/${subf}.collate.bam
  rm 20-Alignment/${subf}/${subf}.sam #20-Alignment/${subf}/${subf}.bam
  gzip 20-Alignment/${subf}/${subf}_*.fastq
done

## Variant Calling (variant_calling.py)

# Nested folders
# while read line ; do mkdir reads_"$line" ; mkdir reads_"$line"/"$line" ; mv 16S_"$line".fasta reads_"$line"/"$line"/. ; done<list.txt
mkdir -p 30-VariantCalling
python3 variant_calling.py -r 16S_"$line".fasta -s reads_"$line" -o output_"$line" -n project_"$line" -f .fastq.gz -p 2 -c 8

