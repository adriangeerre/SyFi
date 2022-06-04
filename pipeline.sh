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


### Functions

## FUNCTION: Variant Calling
function variantCalling() {
  # Variables
  bam_file=$1
  reference_fasta=$2
  reference_dict=$3
  input_bam=$4
  projectname=$5
  output_dir=$6
  threads=$7

  # Intro messages
  printf "\n"
  echo "Performing Variant Calling"

  # Kmer plots (avoid for the moment)

  # Mark duplicates (Mapping is already done)
  printf "\n"
  echo "Mark BAM duplicates"
  markDuplicates ${bam_file} ${output_dir} ${threads}

  # Haplotype caller
  printf "\n"
  echo "Haplotype calling"
  haplotypeCaller ${reference_fasta} ${reference_dict} ${output_dir} ${input_bam}

  # Joint genotyping
  printf "\n"
  echo "Join genotyping"
  jointGenotype ${reference_fasta} ${output_dir} ${projectname} ${threads}

}

# FUNCTION: Mark Duplicates
function markDuplicates() {
  # Arguments
  bam_file=$1
  output_dir=$2
  threads=$3

  # Variables
  sample=$(basename ${bam_file} | sed 's/.bam//g')
  output_fn="${output_dir}/mapped_filtered/${sample}_stats.txt)"

  # Execution
  gatk MarkDuplicates -I ${bam_file} -O ${output_dir}/${sample}.filtered.bam -M ${output_dir}/${sample}.filtered.bam-metrics.txt -AS --REMOVE_DUPLICATES true --VERBOSITY ERROR --CREATE_INDEX true --TMP_DIR ${output_dir}

  # Index BAM
  # samtools index -@ ${threads} -b ${filtered_bam} 

  # Plot BAM
  samtools stats -@ ${threads} -d ${bam_file} > ${output_fn}
  plot-bamstats -p ${output_dir}/mapped_filtered/ ${output_fn}

}

# FUNCTION: Haplotype Caller
function haplotypeCaller() {
  # Arguments
  reference_fasta=$1
  reference_dict=$2
  output_dir=$3
  input_bam=$4

  # Variables (-t genomic -m joint -a False -i False)
  sample=$(basename ${input_bam} | cut -d "." -f 1)
  erc_mode="GVCF"
  min_thr=30
  output_mode="EMIT_VARIANTS_ONLY"
  ploidy=2
  
  # Create dictionary
  gatk CreateSequenceDictionary -R {reference_fasta} -O {reference_dict}

  # Add read group to BAM
  gatk AddOrReplaceReadGroups -I ${input_bam} -O ${output_dir}/${sample}.filtered.readgroup.bam -LB lib1 -PL ILLUMINA -PU unit1 -SM ${sample}
  samtools index -b ${output_dir}/${sample}.filtered.readgroup.bam ${output_dir}/${sample}.filtered.readgroup.bai -@ ${threads}

  # Execution
  gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=${output_dir}/genotyped/" HaplotypeCaller -ERC ${erc_mode} --verbosity ERROR -VS LENIENT --native-pair-hmm-threads ${threads} -ploidy ${ploidy} -stand-call-conf ${min_thr} -I ${output_dir}/${sample}.filtered.readgroup.bam -O ${output_dir}/genotyped/${sample}.g.vcf.gz -R ${reference_fasta} --output-mode ${output_mode}
}

# FUNCTION: Joint Genotyping
function jointGenotype() {
  # Arguments
  reference_fasta=$1
  output_dir=$2
  projectname=$3
  threads=$4

  # Variables
  input_gvcf="${output_dir}/genotyped/${sample}.g.vcf.gz"
  intervals_fn="${output_dir}/reference/intervals.list"
  workspace_dir="${output_dir}/variants/db_workspace"
  merge_intervals="--merge-input-intervals"
  database=${workspace_dir}
  output_vcf="${output_dir}/variants/${projectname}.vcf.gz"
  comp_fn="${output_dir}/variants/${projectname}.vchk"
  plot_dir="${output_dir}/variants/${projectname}/"

  # Create intervals
  cat ${output_dir}/${sample}.fasta.fai | awk '{print $1":1-"$2}' > ${intervals_fn}

  # Execution
  # tabix -p vcf {input_gvcf} # Already done by gatk

  gatk --java-options "-Xmx32g -Xms8g -Djava.io.tmpdir={output_dir}/variants/" GenomicsDBImport -V ${input_gvcf} -L ${intervals_fn} --tmp-dir ${output_dir}/variants/ --genomicsdb-workspace-path ${workspace_dir} --batch-size 70 --seconds-between-progress-updates 120 --reader-threads ${threads} ${merge_intervals} # Run genotype (Create database)

  gatk --java-options "-Xmx32g -Djava.io.tmpdir={output_dir}/variants/" GenotypeGVCFs -V gendb://${database} -R ${reference_fasta} -O ${output_vcf} --tmp-dir ${output_dir}/variants/ -L ${intervals_fn} -G StandardAnnotation --seconds-between-progress-updates 120 # Run genotype
  bcftools query -l ${output_vcf}

  bcftools stats -F ${reference_fasta} -s- ${output_vcf} > ${comp_fn}
  plot-vcfstats -p ${plot_dir} -s ${comp_fn}

}


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
  # Alignment
  bwa-mem2 index 11-Sequences/${subf}/${subf}.fasta
  bwa-mem2 mem 11-Sequences/${subf}/${subf}.fasta ${INPUT_FOLDER}/${subf}/${subf}_R1.fastq.gz ${INPUT_FOLDER}/${subf}/${subf}_R2.fastq.gz -t ${THREADS} > 20-Alignment/${subf}/${subf}.sam
  # Sam to BAM
  samtools view -bS 20-Alignment/${subf}/${subf}.sam -@ ${THREADS} > 20-Alignment/${subf}/${subf}.bam
  # Sort BAM (Coordinate) for Variant Call
     20-Alignment/${subf}/${subf}.sort.sam -O bam 20-Alignment/${subf}/${subf}.sam -@ ${THREADS}
  # Obtain Fastq's
  samtools collate 20-Alignment/${subf}/${subf}.bam 20-Alignment/${subf}/${subf}.collate
  samtools fastq -1 20-Alignment/${subf}/${subf}_R1.fastq -2 20-Alignment/${subf}/${subf}_R2.fastq -s 20-Alignment/${subf}/${subf}_leftover.fastq 20-Alignment/${subf}/${subf}.collate.bam
  # Clean folder
  rm 20-Alignment/${subf}/${subf}.sam #20-Alignment/${subf}/${subf}.bam
  gzip 20-Alignment/${subf}/${subf}_*.fastq
done

## Variant Calling (variant_calling.py)

# Nested folders
# while read line ; do mkdir reads_"$line" ; mkdir reads_"$line"/"$line" ; mv 16S_"$line".fasta reads_"$line"/"$line"/. ; done<list.txt
mkdir -p 30-VariantCalling
for subf in $(ls ${INPUT_FOLDER}); do
  mkdir -p 30-VariantCalling/${subf}

  variantCalling ...

  #python3 00-scripts/variant_calling/variant_calling.py -r 11-Sequences/${subf}/${subf}.fasta -s 20-Alignment -o 30-VariantCalling/${subf} -n project_${subf} -f .fastq.gz -p 2 -c 8
done




# Phasing
mkdir -p 40-Phasing
bash genome_phase.sh -s ${subf} -t 2 -r 11-Sequences/${subf}/${subf}.fasta -v 30-VariantCalling/${subf}/variants/*.vcf.gz -o 40-Phasing/${subf}
