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

  # Kmer plots (Avoided at the moment)

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
  gatk MarkDuplicates -I ${bam_file} -O ${output_dir}/mapped_filtered/${sample}.filtered.bam -M ${output_dir}/mapped_filtered/${sample}.filtered.bam-metrics.txt -AS --REMOVE_DUPLICATES true --VERBOSITY ERROR --CREATE_INDEX true --TMP_DIR ${output_dir}/mapped_filtered

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
  gatk AddOrReplaceReadGroups -I ${input_bam} -O ${output_dir}/genotyped/${sample}.filtered.readgroup.bam -LB lib1 -PL ILLUMINA -PU unit1 -SM ${sample}
  samtools index -b ${output_dir}/genotyped/${sample}.filtered.readgroup.bam ${output_dir}/genotyped/${sample}.filtered.readgroup.bai -@ ${threads}

  # Execution
  gatk --java-options "-Xmx${mem}g -Djava.io.tmpdir=${output_dir}/genotyped/" HaplotypeCaller -ERC ${erc_mode} --verbosity ERROR -VS LENIENT --native-pair-hmm-threads ${threads} -ploidy ${ploidy} -stand-call-conf ${min_thr} -I ${output_dir}/genotyped/${sample}.filtered.readgroup.bam -O ${output_dir}/genotyped/${sample}.g.vcf.gz -R ${reference_fasta} --output-mode ${output_mode}
}

# FUNCTION: Joint Genotyping
function jointGenotype() {
  # Arguments
  reference_fasta=$1
  output_dir=$2
  projectname=$3
  threads=$4

  # Variables
  sample=$(echo ${projectname} | sed 's/project_//g')
  input_gvcf="${output_dir}/genotyped/${sample}.g.vcf.gz"
  intervals_fn="${output_dir}/reference/intervals.list"
  workspace_dir="${output_dir}/variants/db_workspace"
  merge_intervals="--merge-input-intervals"
  database=${workspace_dir}
  output_vcf="${output_dir}/variants/${projectname}.vcf.gz"
  comp_fn="${output_dir}/variants/${projectname}.vchk"
  plot_dir="${output_dir}/variants/${projectname}/"

  # Create intervals
  cat 11-Sequences/${sample}/${sample}.fasta.fai | awk '{print $1":1-"$2}' > ${intervals_fn}

  # Execution
  # tabix -p vcf {input_gvcf} # Already done by gatk

  gatk --java-options "-Xmx32g -Xms8g -Djava.io.tmpdir={output_dir}/variants/" GenomicsDBImport -V ${input_gvcf} -L ${intervals_fn} --tmp-dir ${output_dir}/variants/ --genomicsdb-workspace-path ${workspace_dir} --batch-size 70 --seconds-between-progress-updates 120 --reader-threads ${threads} ${merge_intervals} # Run genotype (Create database)

  gatk --java-options "-Xmx32g -Djava.io.tmpdir={output_dir}/variants/" GenotypeGVCFs -V gendb://${database} -R ${reference_fasta} -O ${output_vcf} --tmp-dir ${output_dir}/variants/ -L ${intervals_fn} -G StandardAnnotation --seconds-between-progress-updates 120 # Run genotype
  bcftools query -l ${output_vcf}

  bcftools stats -F ${reference_fasta} -s- ${output_vcf} > ${comp_fn}
  plot-vcfstats -p ${plot_dir} -s ${comp_fn}

}


### Execution

#-------------- #
# I. Haplotypes #
#-------------- #

## Mapping (BLAST)
mkdir -p 10-Blast 11-Sequences
for subf in $(ls ${INPUT_FOLDER}); do
  mkdir -p 11-Sequences/${subf}
  blastn -query ${INPUT_FOLDER}/${subf}/*.fasta -subject ${target} -strand both -outfmt "6 std qseq" > 10-Blast/${subf}.tsv
  printf ">${subf}\n$(cat 10-Blast/${subf}.tsv | head -n 1 | cut -f 13 | sed 's/-//g')" > 11-Sequences/${subf}/${subf}.fasta
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
  # Obtain BAM of mapped reads (properly pair)
  samtools view -q 30 -f 0x2 20-Alignment/${subf}/${subf}.bam > 20-Alignment/${subf}/${subf}.mapped.bam 
  # Obtain Fastq's
  samtools collate 20-Alignment/${subf}/${subf}.mapped.bam 20-Alignment/${subf}/${subf}.collate
  samtools fastq -1 20-Alignment/${subf}/${subf}_R1.fastq -2 20-Alignment/${subf}/${subf}_R2.fastq -s 20-Alignment/${subf}/${subf}_leftover.fastq 20-Alignment/${subf}/${subf}.collate.bam
  # Clean folder
  rm 20-Alignment/${subf}/${subf}.sam #20-Alignment/${subf}/${subf}.bam
  gzip 20-Alignment/${subf}/${subf}_*.fastq
done

## Variant Calling & Phasing

# Folders
mkdir -p 30-VariantCalling 40-Phasing
for subf in $(ls ${INPUT_FOLDER}); do
  # Create folders
  mkdir -p 30-VariantCalling/${subf}
  mkdir -p 30-VariantCalling/${subf}/mapped_filtered 30-VariantCalling/${subf}/genotyped 30-VariantCalling/${subf}/reference 30-VariantCalling/${subf}/variants

  # Variant call
  variantCalling 20-Alignment/${subf}/${subf}.sort.bam 11-Sequences/${subf}/${subf}.fasta 11-Sequences/${subf}/${subf}.dict 30-VariantCalling/${subf}/mapped_filtered/CHA0_modified.filtered.bam project_${subf} 30-VariantCalling/${subf} ${THREADS}

  # Phasing
  bash 00-scripts/genome_phase.sh -s ${subf} -t ${THREADS} -r 11-Sequences/${subf}/${subf}.fasta -v 30-VariantCalling/${subf}/variants/${subf}.vcf.gz -o 40-Phasing/${subf}

  # Concatenate haplotypes and rename headers
  mkdir -p 50-Haplotypes/${subf}
  hnum=$(cat 40-Phasing/${subf}/${subf}_assembly_h*.fasta | sed 's/^> />/g' | grep "^>" | wc -l)

  for n in $(seq ${hnum})
  do
    if [ ${n} -eq "1" ]; then
      cat 40-Phasing/${subf}/${subf}_assembly_h*.fasta | sed 's/^> />/g' | sed -z "s|CHA0_modified|CHA0_modified_h${n}|${n}" > tmp
      mv tmp 50-Haplotypes/${subf}/${subf}_haplotypes.fasta
    else
      cat 50-Haplotypes/${subf}/${subf}_haplotypes.fasta | sed -z "s|CHA0_modified|CHA0_modified_h${n}|${n}" > tmp
      mv tmp 50-Haplotypes/${subf}/${subf}_haplotypes.fasta
    fi
  done

  # Seqkit duplicate removal
  seqkit rmdup -s 50-Haplotypes/${subf}/${subf}_haplotypes.fasta > 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta

  # Kallisto (Only applied to strains with more that one haplotype)
  if [ $(grep "^>" 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta | wc -l) -gt "2" ]
  then
    mkdir -p 60-Kallisto
    kallisto index -i 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta
    kallisto quant -i 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx -o 60-Kallisto/${subf} 00-Data/${subf}/${subf}_R1.fastq.gz 00-Data/${subf}/${subf}_R2.fastq.gz
  fi
  
  # Filter haplotypes
  Rscript filterHaplotypes.R 60-Kallisto/${subf}/abundance.tsv
done

# Merge Kallisto output
cat 60-Kallisto/*/abundance.tsv > 60-Kallisto/kallisto_output.tsv

# --------------- #
# II. Copy Number #
# --------------- #

for subf in $(ls ${INPUT_FOLDER}); do
  # Variables
  lgen=$(cat 00-Data/${subf}/${subf}.fasta | grep -v "^>" | tr -d "\n" | wc -c)
  braw=$(zcat 00-Data/${subf}/${subf}_R[12].fastq.gz | paste - - - - | cut -f 2 | tr -d "\n" | wc -c)
  l16S=$(cat 11-Sequences/${subf}/${subf}.fasta | grep -v "^>" | tr -d "\n" | wc -c)
  b16S=$(zcat 20-Alignment/${subf}/${subf}_R[12].fastq.gz | paste - - - - | cut -f 2 | tr -d "\n" | wc -c)

  # Compute ratio (round <1 to 1)
  cnum=$(echo "(${b16S}/${l16S}) / (${braw}/${lgen})" | bc -l)
  if (( $(echo "${cnum} < 0.5" | bc -l) )); then
    cnum="1"
  elif (( $(echo "${cnum} > 25.5" | bc -l) )); then
    cnum="25"
  fi

  # File
  printf "${subf}\t${lgen}\t${braw}\t${l16S}\t${b16S}\t${cnum}\n" > copy_number.tsv
done

# ---------------- #
# III. Integration #
# ---------------- #

#

# ---------------- #
# IV. Fingerprints #
# ---------------- #
