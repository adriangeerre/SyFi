#!/bin/bash

# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Gijs Selten, Florian Lamouche and Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
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
  array_keep=("None" "BAMs" "All")
  array_verbose=("Quiet" "Sample" "All")
  array_force=("None" "All" "Skipped" "Failed")

  # Variables
  if [[ $1 != "help" && $9 != 0 ]]; then
    folder=$1
    target=$2
    lendev=$3
    cutoff=$4
    threads=$5
    minmem=$6
    maxmem=$7
    keepfiles=${array_keep[${8}]}
    verbose=${array_keep[${9}]}
    force=${array_force[${10}]}

    # Logo
    echo ""
    echo "${w}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${n}"
    echo "${w}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${n}"
    echo "${w}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${n}"
    echo "${w}|____|/___/    /____/   |___|¯¯    |________|${n}"
    echo ""
    echo "----------------------------------------------"
    echo "Summary:"
    echo "    Folder: ${folder}"
    echo "    Target: ${target}"
    echo "    Deviation length: ${lendev}"
    echo "    Threads: ${threads}"
    echo "    Minimum memory: ${minmem} GB"
    echo "    Maximum memory: ${maxmem} GB"
    echo "    Keep files: ${keepfiles}"
    echo "    Verbose: ${verbose}"
    echo "    Force: ${force}"
    echo ""
    printf "Start:\n"
    date "+    date: %d/%m/%Y"
    date "+    time: %H:%M:%S"
    echo ""
    echo "Progress:"
    checkProgress ${folder}
    echo "----------------------------------------------"
    echo ""
  elif [[ $9 != 0 ]]; then
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
  echo "${w}Usage: ./$0 -i <INPUT_FOLDER> -s <SEARCH_TARGET> -t <THREADS>${n}"
  printf "\n"
  echo "${w}REQUIRED:${n}"
  echo "# Input"
  echo "  -i  | --input_folder     Folder containing input genomes and reads. The software assumes that the folder contains sub-folders for each strain. For more details, execute <pipeline --folder_structure>."
  echo "  -s  | --search_target    Genomic region of interest in fasta format, e.g., 16S."
  printf "\n"
  echo "${w}OPTIONAL:${n}"
  echo "# Haplotype deviation:"
  echo "  -l  | --len_deviation    Total base-pairs for the haplotypes to deviate from the target length upstream and downstream (defaut: 100 bp)."
  echo "  -c | --cutoff            Maximum ratio deviation between haplotypes per sample. This parameter defined how much can an haplotype deviate from the minimum haplotype ratio (default: 25)."
  printf "\n"
  echo "# Input extension:"
  echo "  --fasta-extension        Reference file extension (default: fasta)."
  echo "  --fastq-extension        Illumina reads file extension (default: fastq.gz)."
  printf "\n"
  echo "# Computation:"
  echo "  -t  | --threads          Number of threads (default: 1)."
  echo "  -mn | --min_memory       Minimum memory required in GB (default: 4GB)."
  echo "  -mx | --max_memory       Maximum memory required in GB (default: 8GB)."
  printf "\n"
  echo "# Output options:"
  echo "  -k  | --keep_files       Keep temporary files [0: None, 1: BAM's, or 2: All] (default: 0)."
  echo "  -v  | --verbose          Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2)."
  echo "  -f  | --force            Force re-computation of computed samples [0: None, 1: All, 2: Skipped, or 3: Failed] (default: 0)."
  printf "\n"
  echo "# Display:"
  echo "  -h  | --help             Display help."
  echo "  --citation               Display citation."
  echo "  --folder_structure       Display required folder structure."
  printf "\n"
}

function folder_structure()
{
  printf "\n"
  echo "The pipeline assumes that the genomes and reads are organized in sub-folders inside of the input folder (-i | --input_folder). Each sub-folder should contain the genome (.fasta) and the reads (.fastq.gz)."
  echo "For example:"
  echo ""
  echo "input_folder/
          └── strain_1
              ├── strain_1_R1.fastq.gz
              ├── strain_1_R2.fastq.gz
              └── strain_1.fasta
          └── strain_2
              ├── strain_2_R1.fastq.gz
              ├── strain_2_R2.fastq.gz
              └── strain_2.fasta"
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
FAEXT='fasta'
FQEXT='fastq.gz'
MIN_MEM=4
MAX_MEM=8
BPDEV=300
CUTOFF=25
FORCE=0
KEEPF=0
VERBOSE=2

# Define software path
SYFI_BASE=$(/usr/bin/dirname $(/usr/bin/realpath ${0}))
GENOME_PHASE=$(echo ${SYFI_BASE}/src/genome_phase.sh)
FILTER_HAPLOTYPES=$(echo ${SYFI_BASE}/src/filterHaplotypes.R)
INTEGRATION=$(echo ${SYFI_BASE}/src/Integration.R)

### Parameters

# display usage if 
if [[ $# -lt 4  && $1 != "-h" && $1 != "--help" && $1 != "-c" && $1 != "--citation" && $1 != "--folder_structure" ]]; then
  printf "\n${red}ERROR:${normal} You must provided at least the assembly and alignment files.\n"
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
    -l | --len_deviation)
      shift
      BPDEV=$1
      shift
      ;;
    -c | --cutoff)
      shift
      CUTOFF=$1
      shift
      ;;
    --fasta_extension)
      shift
      FAEXT=$1
      shift
      ;;
    --fastq_extension)
      shift
      FQEXT=$1
      shift
      ;;
    -t | --threads)
      shift
      THREADS=$1
      shift
      ;;
    -mn | --min_mem)
      shift
      MIN_MEM=$1
      shift
      ;;
    -mx | --max_mem)
      shift
      MAX_MEM=$1
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
    --folder_structure)
      folder_structure
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
software=(blastn bwa-mem2 samtools gzip spades.py seqtk seqkit kallisto Rscript gatk plot-bamstats bcftools plot-vcfstats)
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

# Check: Length deviation above zero
  tl=$(grep -v "^>" ${SEARCH_TARGET} | wc -c)
  min_tl=$((tl-${BPDEV}))
  if [ ${min_tl} -le 0 ]; then
    printf "\n${red}ERROR:${normal} length deviation is equal or larger than the target length (target length: ${tl}).\n\n"
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

  # Check file
  if [ -e progress.txt ]; then
    total=$(ls $INPUT_FOLDER | wc -l)
    success=$(grep "Success" progress.txt | wc -l)
    skipped=$(grep "Skipped" progress.txt | wc -l)
    failed=$(grep "Failed" progress.txt | wc -l)
    printf "\tSuccess: ${g} ${success} ${n}\n\tSkipped: ${y} ${skipped} ${n}\n\tFailed: ${r}  ${failed} ${n}\n\tTotal:    ${total}\n"
  else
    total=$(ls $INPUT_FOLDER | wc -l)
    printf "\tSuccess: ${g} 0 ${n}\n\tSkipped: ${y} 0 ${n}\n\tFailed: ${r}  0 ${n}\n\tTotal:    ${total}\n"
  fi

}

### Functions

## FUNCTION: Variant Calling
function variantCalling() {
  # Variables
  bam_file=$1
  reference_fasta=$2
  reference_dict=$3
  input_bam=$4
  output_dir=$5
  threads=$6
  min_mem=$7
  max_mem=$8

  # Intro messages
  printf "\n"
  echo "Performing Variant Calling"

  # Mark duplicates (Mapping is already done)
  printf "\n"
  echo "Mark BAM duplicates"
  markDuplicates ${bam_file} ${output_dir} ${threads}

  # Haplotype caller
  printf "\n"
  echo "Haplotype calling"
  haplotypeCaller ${reference_fasta} ${reference_dict} ${output_dir} ${input_bam} ${max_mem}

  # Joint genotyping
  printf "\n"
  echo "Join genotyping"
  jointGenotype ${reference_fasta} ${output_dir} ${threads} ${max_mem} ${min_mem}

}

# FUNCTION: Mark Duplicates
function markDuplicates() {
  # Arguments
  bam_file=$1
  output_dir=$2
  threads=$3

  # Variables
  sample=$(basename ${bam_file} | sed 's/.rebuild.sort.bam//g')
  output_fn="${output_dir}/mapped_filtered/${sample}_stats.txt"

  # Execution
  gatk MarkDuplicates -I ${bam_file} -O ${output_dir}/mapped_filtered/${sample}.filtered.bam -M ${output_dir}/mapped_filtered/${sample}.filtered.bam-metrics.txt -AS --REMOVE_DUPLICATES true --VERBOSITY ERROR --CREATE_INDEX true --TMP_DIR ${output_dir}/mapped_filtered

  # Plot BAM
  samtools stats -d ${bam_file} > ${output_fn}
  plot-bamstats -p ${output_dir}/mapped_filtered/ ${output_fn}
}

# FUNCTION: Haplotype Caller
function haplotypeCaller() {
  # Arguments
  reference_fasta=$1
  reference_dict=$2
  output_dir=$3
  input_bam=$4
  max_mem=$5

  # Variables (-t genomic -m joint -a False -i False)
  sample=$(basename ${input_bam} | cut -d "." -f 1)
  erc_mode="GVCF"
  min_thr=30
  output_mode="EMIT_VARIANTS_ONLY"
  ploidy=2
  
  # Create dictionary
  gatk CreateSequenceDictionary -R ${reference_fasta} -O ${reference_dict}

  # Add read group to BAM
  gatk AddOrReplaceReadGroups -I ${input_bam} -O ${output_dir}/genotyped/${sample}.filtered.readgroup.bam -LB lib1 -PL ILLUMINA -PU unit1 -SM ${sample}
  samtools index -b ${output_dir}/genotyped/${sample}.filtered.readgroup.bam ${output_dir}/genotyped/${sample}.filtered.readgroup.bai -@ ${threads}

  # Index fasta (without -o for samtools version 1.7 or below)
  samtools faidx ${reference_fasta} #-o ${reference_fasta}.fai

  # Execution
  gatk --java-options "-Xmx${max_mem}g -Djava.io.tmpdir=${output_dir}/genotyped/" HaplotypeCaller -ERC ${erc_mode} --verbosity ERROR -VS LENIENT --native-pair-hmm-threads ${threads} -ploidy ${ploidy} -stand-call-conf ${min_thr} -I ${output_dir}/genotyped/${sample}.filtered.readgroup.bam -O ${output_dir}/genotyped/${sample}.g.vcf.gz -R ${reference_fasta} --output-mode ${output_mode}
}

# FUNCTION: Joint Genotyping
function jointGenotype() {
  # Arguments
  reference_fasta=$1
  output_dir=$2
  threads=$3
  max_mem=$4
  min_mem=$5

  # Variables
  sample=$(basename ${reference_fasta} | sed 's/.fasta//g')
  input_gvcf="${output_dir}/genotyped/${sample}.g.vcf.gz"
  intervals_fn="${output_dir}/reference/intervals.list"
  workspace_dir="${output_dir}/variants/db_workspace"
  merge_intervals="--merge-input-intervals"
  database=${workspace_dir}
  output_vcf="${output_dir}/variants/${sample}.vcf.gz"
  comp_fn="${output_dir}/variants/${sample}.vchk"
  plot_dir="${output_dir}/variants/${sample}/"

  # Create intervals
  cat 20-Alignment/${sample}/${sample}.fasta.fai | awk '{print $1":1-"$2}' > ${intervals_fn}

  # Execution
  # tabix -p vcf {input_gvcf} # Already done by gatk

  gatk --java-options "-Xmx${max_mem}g -Xms${min_mem}g -Djava.io.tmpdir={output_dir}/variants/" GenomicsDBImport -V ${input_gvcf} -L ${intervals_fn} --tmp-dir ${output_dir}/variants/ --genomicsdb-workspace-path ${workspace_dir} --batch-size 70 --seconds-between-progress-updates 120 --reader-threads ${threads} ${merge_intervals} # Run genotype (Create database)

  gatk --java-options "-Xmx${max_mem}g -Djava.io.tmpdir={output_dir}/variants/" GenotypeGVCFs -V gendb://${database} -R ${reference_fasta} -O ${output_vcf} --tmp-dir ${output_dir}/variants/ -L ${intervals_fn} -G StandardAnnotation --seconds-between-progress-updates 120 # Run genotype
  # bcftools query -l ${output_vcf}

  bcftools stats -F ${reference_fasta} -s- ${output_vcf} > ${comp_fn}
  plot-vcfstats -p ${plot_dir} -s ${comp_fn}
}

# FUNCTION: Copy Number
function copyNumber() {
  # Arguments
  INPUT_FOLDER=$1
  subf=$2

  if [ ${VERBOSE} -eq 2 ]; then printf "Copy number; "; fi

  # Log
  printf "\n\n### Target Copy Number ###\n\n" >> 01-Logs/log_${subf}.txt

  # Header
  printf "Strain\tGenome_length\tGenome_nbases\tTarget_length\tTarget_nbases\tCopy_number\n" > 60-Integration/${subf}/copy_number.tsv

  # Variables
  # Get length of assembly
  lgen=$(cat ${INPUT_FOLDER}/${subf}/${subf}.fasta | grep -v "^>" | tr -d "\n" | wc -c)
  # Get number of bases in assembly reads (Speed up?)
  braw=$(zcat ${INPUT_FOLDER}/${subf}/${subf}_R[12].fastq.gz | paste - - - - | cut -f 2 | tr -d "\n" | wc -c)

  # Get length of longest recovered target
  l16S=$(cut -f 14 20-Alignment/${subf}/flanking/max.tsv)

  # Number of bases in selected reads
  start=$(cut -f 7 20-Alignment/${subf}/flanking/max.tsv)
  end=$(cut -f 8 20-Alignment/${subf}/flanking/max.tsv)
  b16S=$(bedtools genomecov -d -ibam 20-Alignment/${subf}/${subf}.rebuild.sort.bam | awk -v start=${start} '($2 > start)' | awk -v end=${end} '($2 < end)' | awk '{sum+=$3;} END{print sum;}')

  # Compute ratio (round <1 to 1)
  cnum=$(echo "(${b16S}/${l16S}) / (${braw}/${lgen})" | bc -l)
  if (( $(echo "${cnum} < 0.5" | bc -l) )); then
    cnum="1"
  elif (( $(echo "${cnum} > 25.5" | bc -l) )); then
    cnum="25"
  fi

  # File
  printf "${subf}\t${lgen}\t${braw}\t${l16S}\t${b16S}\t${cnum}\n" >> 60-Integration/${subf}/copy_number.tsv
}

# FUNCTION: Fingerprint
function fingerPrint() {
  # Variables
  mode=$1 # unique/multiple
  subf=$2

  # Create folder
  mkdir -p 70-Fingerprints/${subf}

  if [ $mode == "unique" ];then
    # Log
    if [ ${VERBOSE} -eq 2 ]; then printf "Fingerprint [unique]\n"; fi
    printf "\n\n### Fingerprint ###\n\n" >> 01-Logs/log_${subf}.txt

    # Copy recovered target
    cp 20-Alignment/${subf}/${subf}.fasta 70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta

    # Rename fasta header
    cat 70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta | sed "s/^>NODE.*/>${subf}_all_haplotypes/g" >> 70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta

  elif [ $mode == "multiple" ];then
    # Log
    if [ ${VERBOSE} -eq 2 ]; then printf "Fingerprint\n"; fi
    printf "\n\n### Fingerprint ###\n\n" >> 01-Logs/log_${subf}.txt

    # Separate fasta into multiples fasta
    seqkit split -i 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta -O 70-Fingerprints/${subf} &>> 01-Logs/log_${subf}.txt

    # Rename haplotypes files
    for i in $(ls -d 70-Fingerprints/${subf}/*); do mv ${i} $(echo ${i} | sed "s/clean_${subf}_haplotypes.part_//g"); done &>> 01-Logs/log_${subf}.txt

    # Keep haplotypes filtered in integration
    keep=($(tail -n +2 60-Integration/${subf}/integration.tsv | cut -f 1))
    for i in $(ls 70-Fingerprints/${subf} | sed 's/.fasta//g'); do
      if [[ ! "${keep[*]}" =~ "${i}" ]]
      then
        rm -f 70-Fingerprints/${subf}/${i}.fasta
      fi
    done

    # Concatenate haplotypes (NF-1 to obtain "per_haplotype")
    cat 60-Integration/${subf}/integration.tsv | awk '{print $1 "/" $(NF-1)}' | tail -n +2 > 60-Integration/${subf}/tmp.tsv
    for h in $(cat 60-Integration/${subf}/tmp.tsv); do
      seq=$(echo ${h} | cut -d "/" -f 1)
      num=$(echo ${h} | cut -d "/" -f 2)
      for n in $(seq ${num}); do
        if [[ ${n} == ${num} && ${n} != 1 ]]; then
          grep -v "^>" 70-Fingerprints/${subf}/${seq}.fasta >> 70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
        else
          grep -v "^>" 70-Fingerprints/${subf}/${seq}.fasta >> 70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
          echo "NNNNNNNNNN" >> 70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
        fi
      done 
    done
    printf ">${subf}_all_haplotypes\n" > 70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta
    cat 70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta | tr -d "\n" >> 70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta
  fi
}

function CleanFiles() {
  # Variables
  mode=$1
  subf=$2
  KEEPF=$3

  # Clean folder
  # MODE: Final
  if [ ${mode} == "Final" ]; then
    if [ ${KEEPF} -eq 0 ]; then
      # BAMs
      rm -rf 20-Alignment/${subf}/${subf}.bam 20-Alignment/${subf}/${subf}.*.bam
      rm -rf 30-VariantCalling/${subf}/mapped_filtered/*.ba*
      rm -rf 40-Phasing/${subf}/mapped
    fi
    if [ ${KEEPF} -lt 2 ]; then
      # 11-Sequences
      rm -rf 11-Sequences/${subf}/${subf}.fasta.*
      # 20-Alignment
      rm -rf 20-Alignment/${subf}/${subf}.sam 20-Alignment/${subf}/${subf}.dict 20-Alignment/${subf}/${subf}.fasta.* 20-Alignment/${subf}/spades 20-Alignment/${subf}/${subf}.rebuild.sam
      # 30-VariantCalling
      rm -rf 30-VariantCalling/${subf}/genotyped 30-VariantCalling/${subf}/mapped_filtered 30-VariantCalling/${subf}/reference  30-VariantCalling/${subf}/variants/tmp_* 30-VariantCalling/${subf}/variants/db_workspace 30-VariantCalling/${subf}/variants/${subf}
      # 40-Phasing
      rm -rf 40-Phasing/${subf}/${subf}_reads.txt 
      # 50-Haplotypes
      rm -rf 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx
      # 60-Integration
      rm -rf 60-Integration/${subf}/abundance.h5 60-Integration/${subf}/run_info.json 60-Integration/${subf}/tmp.tsv
      # 70-Fingerprints
      rm -rf 70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
    fi
  # MODE: ReadRecov
  elif [ ${mode} == "ReadRecov" ]; then
    if [ ${KEEPF} -eq 0 ]; then
      # BAMs
      rm -rf 20-Alignment/${subf}/${subf}.bam 20-Alignment/${subf}/${subf}.*.bam
    fi
    if [ ${KEEPF} -lt 2 ]; then
      # 11-Sequences
      rm -rf 11-Sequences/${subf}/${subf}.fasta.*
      # 20-Alignment
      rm -rf 20-Alignment/${subf}/${subf}.sam 20-Alignment/${subf}/${subf}.dict 20-Alignment/${subf}/${subf}.fasta.* 20-Alignment/${subf}/spades 20-Alignment/${subf}/${subf}.rebuild.sam
    fi
  # MODE: Phasing
  elif [ ${mode} == "Phasing" ]; then
    if [ ${KEEPF} -eq 0 ]; then
      # BAMs
      rm -rf 30-VariantCalling/${subf}/mapped_filtered/*.ba*
      rm -rf 40-Phasing/${subf}/mapped
    fi
    if [ ${KEEPF} -lt 2 ]; then
      # 11-Sequences
      rm -rf 11-Sequences/${subf}/${subf}.fasta.*
      # 20-Alignment
      rm -rf 20-Alignment/${subf}/${subf}.sam 20-Alignment/${subf}/${subf}.dict 20-Alignment/${subf}/${subf}.fasta.* 20-Alignment/${subf}/spades 20-Alignment/${subf}/${subf}.rebuild.sam
      # 30-VariantCalling
      rm -rf 20-Alignment/${subf}/${subf}.sam 20-Alignment/${subf}/${subf}.dict 20-Alignment/${subf}/${subf}.fasta.* 20-Alignment/${subf}/spades 20-Alignment/${subf}/${subf}.rebuild.sam
      # 40-Phasing
      rm -rf 40-Phasing/${subf}/${subf}_reads.txt
    fi
  # MODE: Integration
  elif [ ${mode} == "Integration" ]; then
    if [ ${KEEPF} -eq 0 ]; then
      # BAMs
      rm -rf 30-VariantCalling/${subf}/mapped_filtered/*.ba*
      rm -rf 40-Phasing/${subf}/mapped
    fi
    if [ ${KEEPF} -lt 2 ]; then
      # 11-Sequences
      rm -rf 11-Sequences/${subf}/${subf}.fasta.*
      # 20-Alignment
      rm -rf 20-Alignment/${subf}/${subf}.sam 20-Alignment/${subf}/${subf}.dict 20-Alignment/${subf}/${subf}.fasta.* 20-Alignment/${subf}/spades 20-Alignment/${subf}/${subf}.rebuild.sam
      # 30-VariantCalling
      rm -rf 20-Alignment/${subf}/${subf}.sam 20-Alignment/${subf}/${subf}.dict 20-Alignment/${subf}/${subf}.fasta.* 20-Alignment/${subf}/spades 20-Alignment/${subf}/${subf}.rebuild.sam
      # 40-Phasing
      rm -rf 40-Phasing/${subf}/${subf}_reads.txt
      # 50-Haplotypes
      rm -rf 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx
      # 60-Integration
      rm -rf 60-Integration/${subf}/abundance.h5 60-Integration/${subf}/run_info.json 60-Integration/${subf}/tmp.tsv
    fi
  fi
}

### Execution

#-------------- #
# I. Haplotypes #
#-------------- #

# Call logo
logo ${INPUT_FOLDER} ${SEARCH_TARGET} ${BPDEV} ${CUTOFF} ${THREADS} ${MIN_MEM} ${MAX_MEM} ${KEEPF} ${VERBOSE} ${FORCE}

for subf in $(ls ${INPUT_FOLDER}); do
 
  ## CHECK: Redirect workflow given status and force
  if [ -f progress.txt ]; then
    status=$(grep -w ${subf} progress.txt | cut -f 2)
    if [[ ${status} == "Success" ]] && [[ ${FORCE} != 1 ]]; then
      continue
    elif [[ ${status} == "Skipped" ]] && [[ ${FORCE} != 2 ]]; then
      continue
    elif [[ ${status} == "Failed" ]] && [[ ${FORCE} != 3 ]]; then
      continue
    elif [[ ${status} == "" ]] && [[ ${FORCE} != 0 ]]; then
      continue
    fi

    # CHECK: Avoid re-computation
    if [[ ${FORCE} != 0 ]]; then
      # Remove line in progress.txt if exists already
      if [[ $(grep -w ${subf} progress.txt | wc -l) -eq 1 ]]; then
        grep -wv ${subf} progress.txt > tmp
        mv tmp progress.txt
      fi

      rm 10-Blast/${subf}.tsv
      # Remove results for sample
      for fld in 11-Sequences 20-Alignment 30-VariantCalling 40-Phasing 50-Haplotypes 60-Integration 70-Fingerprints; do
        rm -rf ${fld}/${subf}
      done
    fi
  fi

  ## CHECK: Log file
  if [[ -f 01-Logs/log_${subf}.txt ]]; then
    # Remove old log file
    rm -rf 01-Logs/log_${subf}.txt
    touch 01-Logs/log_${subf}.txt
  else
    # Touch Log File
    mkdir -p 01-Logs
    touch 01-Logs/log_${subf}.txt
  fi

  ## CHECK: Input in folder
  # Fasta input
  if [[ ! -f ${INPUT_FOLDER}/${subf}/${subf}.${FAEXT} ]]; then
    printf "\n${red}ERROR:${normal} File ${INPUT_FOLDER}/${subf}/${subf}.${FAEXT} missing, please use the \"-e/--extension\" argument if the extension is not \"fasta\"." | tee -a 01-Logs/log_${subf}.txt
    if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} computation will be skipped.\n"; fi
    printf "${subf}\tSkipped\tWrong reference extension\n" >> progress.txt
    continue
  fi

  # Fastq input
  if [[ ! -f ${INPUT_FOLDER}/${subf}/${subf}_R1.${FQEXT} || ! -f ${INPUT_FOLDER}/${subf}/${subf}_R2.${FQEXT} ]]; then
    printf "\n${red}ERROR:${normal} One or both illumina reads ${INPUT_FOLDER}/${subf}/${subf}_R[1/2].${FQEXT} are missing, please use the \"-e/--extension\" argument if the extension is not \"fastq-gz\"." | tee -a 01-Logs/log_${subf}.txt
    if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} computation will be skipped.\n"; fi
    printf "${subf}\tSkipped\tWrong reads extension\n" >> progress.txt
    continue
  fi

  if [ ${VERBOSE} -ge 1 ]; then date "+start time: %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt; fi
  if [ ${VERBOSE} -ge 1 ]; then printf "${blue}Sample:${normal} $subf\n" | tee -a 01-Logs/log_${subf}.txt; fi

  ## ------------------------------------
  ## Target recovery from contigs (BLAST)
  ## ------------------------------------

  if [ ${VERBOSE} -eq 2 ]; then printf "Performing: Mapping; "; fi

  # Create folder
  mkdir -p 10-Blast 11-Sequences/${subf}
  # Blastn
  blastn -query ${INPUT_FOLDER}/${subf}/${subf}.${FAEXT} -subject ${SEARCH_TARGET} -strand both -outfmt "6 std qseq" > 10-Blast/${subf}.tsv
  printf ">${subf}\n$(cat 10-Blast/${subf}.tsv | sort -n -k4 | tail -n 1 | cut -f 13 | sed 's/-//g')" > 11-Sequences/${subf}/${subf}.fasta

  # CHECK: Absent target match (Perhaps, remove small hits!)
  if [ $(cat 11-Sequences/${subf}/${subf}.fasta | wc -l) -lt 1 ]; then
    if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} No target was found for ${subf}. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt; fi
    printf "${subf}\tFailed\tNo target recovered from reference\n" >> progress.txt
    date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
    echo ""
    continue
  fi

  ## --------------------------------------------------
  ## Recover pair-end reads for target (BWA & Samtools)
  ## --------------------------------------------------

  # Create folder
  mkdir -p 20-Alignment/${subf}

  if [ ${VERBOSE} -eq 2 ]; then printf "Alignment; "; fi

  # Alignment
  printf "### BWA mapping ###\n\n" > 01-Logs/log_${subf}.txt
  bwa-mem2 index 11-Sequences/${subf}/${subf}.fasta &>> 01-Logs/log_${subf}.txt
  bwa-mem2 mem 11-Sequences/${subf}/${subf}.fasta ${INPUT_FOLDER}/${subf}/${subf}_R1.${FQEXT} ${INPUT_FOLDER}/${subf}/${subf}_R2.${FQEXT} -t ${THREADS} 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.sam

  # Sam to BAM
  samtools view -b 20-Alignment/${subf}/${subf}.sam -@ ${THREADS} 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.bam
  # Sort BAM (Coordinate) for Variant Call
  samtools sort -o 20-Alignment/${subf}/${subf}.sort.bam -O bam 20-Alignment/${subf}/${subf}.bam -@ ${THREADS} 2>> 01-Logs/log_${subf}.txt
  # Obtain BAM of mapped reads (properly pair)
  samtools view -b -q 30 -f 0x2 20-Alignment/${subf}/${subf}.bam 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.mapped.bam

  if [ ${VERBOSE} -eq 2 ]; then printf "Reads recovery; "; fi

  # Obtain Fastq's
  printf "\n\n### Reads recovery ###\n\n" >> 01-Logs/log_${subf}.txt
  samtools collate 20-Alignment/${subf}/${subf}.mapped.bam 20-Alignment/${subf}/${subf}.collate 2>> 01-Logs/log_${subf}.txt
  samtools fastq -1 20-Alignment/${subf}/${subf}_R1.fastq -2 20-Alignment/${subf}/${subf}_R2.fastq -s 20-Alignment/${subf}/${subf}_leftover.fastq 20-Alignment/${subf}/${subf}.collate.bam 2>> 01-Logs/log_${subf}.txt
  gzip -f 20-Alignment/${subf}/${subf}_*.fastq

  ## --------------------------------------------------
  ## Rebuild target from reads (SPAdes, BWA & Samtools)
  ## --------------------------------------------------

  # CHECK: Absent R1/R2 for de novo assembly
  if [[ $(zcat 20-Alignment/${subf}/${subf}_R1.fastq.gz | wc -l) -eq 0 || $( zcat 20-Alignment/${subf}/${subf}_R2.fastq.gz | wc -l) -eq 0 ]]; then
    if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} No target reads were recovered for ${subf}. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt; fi
    printf "${subf}\tFailed\tNo reads recovered for target\n" >> progress.txt
    # Clean files/folders
    CleanFiles ${subf} ${KEEPF} "ReadRecov"
    # Date
    date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
    echo ""
    continue
  fi

  if [ ${VERBOSE} -eq 2 ]; then printf "Contig re-build; "; fi

  # SPAdes (Unambigous nucleotides assembly)
  printf "\n\n### Contig re-build ###\n\n" >> 01-Logs/log_${subf}.txt
	spades.py -1 20-Alignment/${subf}/${subf}_R1.fastq.gz -2 20-Alignment/${subf}/${subf}_R2.fastq.gz -t ${THREADS} -o 20-Alignment/${subf}/spades -m ${MAX_MEM} &>> 01-Logs/log_${subf}.txt

  # Size select SPAdes recovered target
  seqtk seq -L ${min_tl} 20-Alignment/${subf}/spades/contigs.fasta > 20-Alignment/${subf}/spades/contigs.seqtk.fasta
  if [ $(grep "^>" 20-Alignment/${subf}/spades/contigs.seqtk.fasta | wc -l) -eq 0 ]; then
    printf "\n${red}ERROR:${normal} No target was recovered for ${subf} or the target recovered was too small. In the second case, make \"-l/--len_deviation\" larger. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt
    printf "${subf}\tSkipped\tRecovered target length below minimum length\n" >> progress.txt
    date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
    echo ""
    continue
  fi
  
  # Define target without flanking regions
  mkdir -p 20-Alignment/${subf}/flanking
  blastn -subject ${SEARCH_TARGET} -query 20-Alignment/${subf}/spades/contigs.seqtk.fasta -outfmt "6 std sseq" > 20-Alignment/${subf}/flanking/${subf}.target.tsv

  # Size select blast hits
  while read -r line;do
    total=$(echo ${line} | cut -d " " -f 4)
    mist=$(echo ${line} | cut -d " " -f 6)
    len=$((total-mist))
    if [ $len -ge $min_tl ]; then
      printf "${line}\t${len}\n"
    fi
  done < 20-Alignment/${subf}/flanking/${subf}.target.tsv > 20-Alignment/${subf}/flanking/${subf}.target.sizeclean.tsv
  
  # Check if hit was not recovered
  if [ $(cat 20-Alignment/${subf}/flanking/${subf}.target.sizeclean.tsv | wc -l) == 0 ]; then
    printf "\n${red}ERROR:${normal} No target was recovered for ${subf} or the target recovered was too small. In the second case, make \"-l/--len_deviation\" larger. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt
    printf "${subf}\tSkipped\tRecovered target length below minimum length\n" >> progress.txt
    date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
    echo ""
    continue
  fi

  # Select maximum (in all cases)
  cat 20-Alignment/${subf}/flanking/${subf}.target.sizeclean.tsv | sort -n -k14 | tail -n 1 > 20-Alignment/${subf}/flanking/max.tsv
  cut -f 1 20-Alignment/${subf}/flanking/max.tsv > 20-Alignment/${subf}/flanking/max.header.txt

  # Select recovered sequence given header of maximum
  seqtk subseq 20-Alignment/${subf}/spades/contigs.seqtk.fasta 20-Alignment/${subf}/flanking/max.header.txt > 20-Alignment/${subf}/${subf}.fasta

  # Align
  printf "\n\n### Contig re-alignment (BWA) ###\n\n" >> 01-Logs/log_${subf}.txt
  bwa-mem2 index 20-Alignment/${subf}/${subf}.fasta &>> 01-Logs/log_${subf}.txt
  bwa-mem2 mem 20-Alignment/${subf}/${subf}.fasta 20-Alignment/${subf}/${subf}_R1.fastq.gz 20-Alignment/${subf}/${subf}_R2.fastq.gz -t ${THREADS} 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.rebuild.sam
	
  samtools view -b 20-Alignment/${subf}/${subf}.rebuild.sam -@ ${THREADS} 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.rebuild.bam
  samtools sort -o 20-Alignment/${subf}/${subf}.rebuild.sort.bam -O bam 20-Alignment/${subf}/${subf}.rebuild.bam -@ ${THREADS} 2>> 01-Logs/log_${subf}.txt
  #samtools view -b -q 30 -f 0x2 20-Alignment/${subf}/${subf}.rebuild.bam 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.rebuild.mapped.bam 

  ## -------------------------------------------
  ## Variant Calling & Phasing (GATK & BCFTools)
  ## -------------------------------------------

  # CHECK: Absent rebuilt BAM
  if [[ ! -f 20-Alignment/${subf}/${subf}.rebuild.sort.bam || ! -f 20-Alignment/${subf}/${subf}.fasta ]]; then
    if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}}WARNING:${normal} No target reads were recovered for ${subf}. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt; fi
    printf "${subf}\tFailed\tNo reads recovered for target\n" >> progress.txt
    # Clean files/folders
    CleanFiles ${subf} ${KEEPF} "ReadRecov"
    # Date
    date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
    echo ""
    continue
  fi

  # Create folders
  mkdir -p 30-VariantCalling/${subf}
  mkdir -p 30-VariantCalling/${subf}/mapped_filtered 30-VariantCalling/${subf}/genotyped 30-VariantCalling/${subf}/reference 30-VariantCalling/${subf}/variants

  if [ ${VERBOSE} -eq 2 ]; then printf "Variant calling; "; fi

  # Variant call
  printf "\n\n### Variant calling ###\n\n" >> 01-Logs/log_${subf}.txt
  variantCalling 20-Alignment/${subf}/${subf}.rebuild.sort.bam 20-Alignment/${subf}/${subf}.fasta 20-Alignment/${subf}/${subf}.dict 30-VariantCalling/${subf}/mapped_filtered/${subf}.filtered.bam 30-VariantCalling/${subf} ${THREADS} ${MIN_MEM} ${MAX_MEM} &>> 01-Logs/log_${subf}.txt

  # CHECK: No variants recovered
  if [[ -f 30-VariantCalling/${subf}/variants/${subf}.vcf.gz && $(zcat 30-VariantCalling/${subf}/variants/${subf}.vcf.gz | grep -v "#" | wc -l) != 0 ]]; then
    # MULTIPLE HAPLOTYPES

    # Create folder
    mkdir -p 40-Phasing/${subf}

    if [ ${VERBOSE} -eq 2 ]; then printf "Phasing; "; fi

    # Phasing
    printf "\n\n### Phasing ###\n\n" >> 01-Logs/log_${subf}.txt
    bash ${GENOME_PHASE} -s ${subf} -t ${THREADS} -r 20-Alignment/${subf}/${subf}.fasta -v 30-VariantCalling/${subf}/variants/${subf}.vcf.gz -i 20-Alignment -o 40-Phasing/${subf} &>> 01-Logs/log_${subf}.txt

    # CHECK: Absent haplotypes
    if [[ ! -f 40-Phasing/${subf}/${subf}_assembly_h1.fasta || ! -f 40-Phasing/${subf}/${subf}_assembly_h2.fasta ]]; then
      if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} No haplotypes were recovered for ${subf}. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt; fi
      printf "${subf}\tFailed\tNo haplotypes were found\n" >> progress.txt
      # Clean files/folders
      CleanFiles ${subf} ${KEEPF} "Phasing"
      # Date
      date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
      echo ""
      continue
    fi

    # Haplotype length correction: remove short matches
    mkdir -p 50-Haplotypes/${subf}
    seqtk seq -L ${min_tl} <(cat 40-Phasing/${subf}/${subf}_assembly_h*.fasta) > 50-Haplotypes/${subf}/${subf}_haplotypes.fasta 

    # Rename headers
    hnum=$(cat 50-Haplotypes/${subf}/${subf}_haplotypes.fasta | sed 's/^> />/g' | grep "^>" | wc -l)

    for n in $(seq ${hnum})
    do
      if [ ${n} -eq "1" ]; then
        cat 50-Haplotypes/${subf}/${subf}_haplotypes.fasta | sed 's/_[0-9]_length_.*//g' | sed -z "s|NODE|seq_h${n}|${n}" > tmp
        mv tmp 50-Haplotypes/${subf}/${subf}_haplotypes.fasta
      else
        cat 50-Haplotypes/${subf}/${subf}_haplotypes.fasta | sed -z "s|NODE|seq_h${n}|1" > tmp
        mv tmp 50-Haplotypes/${subf}/${subf}_haplotypes.fasta
      fi
    done

    # Haplotype duplicate removal (SeqKit)
    seqkit rmdup -s 50-Haplotypes/${subf}/${subf}_haplotypes.fasta 2>> 01-Logs/log_${subf}.txt > 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta

    ## ------------------------------------
    ## Haplotype abundance ratio (Kallisto)
    ## ------------------------------------

    # CHECK: Absent clean haplotypes
    if [ ! -f 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta ]; then
      if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} No haplotypes were recovered for ${subf}. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt; fi
      printf "${subf}\tSkipped\No haplotypes were recovered\n" > progress.txt
      date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
      echo ""
      continue
    fi

    # Kallisto (Only applied to strains with more that one haplotype)
    if [ $(grep "^>" 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta | wc -l) -gt "1" ]
    then
      # Create folder
      mkdir -p 60-Integration

      if [ ${VERBOSE} -eq 2 ]; then printf "Abundance ratio; "; fi

      # Kallisto
      printf "\n\n### Kallisto ###\n\n" >> 01-Logs/log_${subf}.txt
      kallisto index -i 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta &>> 01-Logs/log_${subf}.txt
      kallisto quant -i 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx -o 60-Integration/${subf} ${INPUT_FOLDER}/${subf}/${subf}_R1.fastq.gz ${INPUT_FOLDER}/${subf}/${subf}_R2.fastq.gz -t ${THREADS} &>> 01-Logs/log_${subf}.txt

      # Filter haplotypes
      cp 60-Integration/${subf}/abundance.tsv 60-Integration/${subf}/abundance.orig.tsv
      Rscript ${FILTER_HAPLOTYPES} -i 60-Integration/${subf}/abundance.tsv -c ${CUTOFF} &>> 01-Logs/log_${subf}.txt

      # --------------- #
      # II. Copy Number #
      # --------------- #

      # Copy number
      copyNumber ${INPUT_FOLDER} ${subf}

      # ---------------- #
      # III. Integration #
      # ---------------- #

      # CHECK: Absent kallisto output
      if [ ! -f 60-Integration/${subf}/abundance.tsv ]; then
        if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} Missing kallisto output for ${subf}. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt; fi
        printf "${subf}\tFailed\tKallisto could not determined the haplotype abundances\n" >> progress.txt
        # Clean files/folders
        CleanFiles ${subf} ${KEEPF} "Integration"
        # Date
        date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
        echo ""
        continue
      fi

      if [ ${VERBOSE} -eq 2 ]; then printf "Integration; "; fi
      # Log
      printf "\n\n### Integration ###\n\n" >> 01-Logs/log_${subf}.txt

      # Integrate step I and II
      Rscript ${INTEGRATION} -r 60-Integration/${subf}/abundance.tsv -c 60-Integration/${subf}/copy_number.tsv -i ${subf} -m "multiple" &>> 01-Logs/log_${subf}.txt

      # ---------------- #
      # IV. Fingerprints #
      # ---------------- #

      # CHECK: Absent integration output
      if [ ! -f 60-Integration/${subf}/integration.tsv ]; then
        if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} Missing integration output for ${subf}. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt; fi
        printf "${subf}\tFailed\tIntegration could not be performed\n" >> progress.txt
        # Clean files/folders
        CleanFiles ${subf} ${KEEPF} "Integration"
        # Date
        date "+end time: %d/%m/%Y - %H:%M:%S" | tee -a 01-Logs/log_${subf}.txt
        echo ""
        continue
      fi

      # Fingerprint
      fingerPrint 'multiple' ${subf}
    fi
  elif [[ -f 30-VariantCalling/${subf}/variants/${subf}.vcf.gz && $(zcat 30-VariantCalling/${subf}/variants/${subf}.vcf.gz | grep -v "#" | wc -l) == 0 ]]; then
    # ONE HAPLOTYPE

    if [ ${VERBOSE} -eq 2 ]; then printf "Abundance ratio [Avoid]; "; fi

    # --------------- #
    # II. Copy Number #
    # --------------- #

    # Create folder
    mkdir -p 60-Integration/${subf}

    # Copy number
    copyNumber ${INPUT_FOLDER} ${subf}

    # ---------------- #
    # III. Integration #
    # ---------------- #

    if [ ${VERBOSE} -eq 2 ]; then printf "Integration [unique]; "; fi
    # Log
    printf "\n\n### Integration ###\n\n" >> 01-Logs/log_${subf}.txt

    # Create folder
    mkdir -p 60-Integration/${subf}

    # Integrate only step II 
    Rscript 00-Scripts/Integration.R -r "None" -c 60-Integration/${subf}/copy_number.tsv -i ${subf} -m 'unique' &>> 01-Logs/log_${subf}.txt

    # ---------------- #
    # IV. Fingerprints #
    # ---------------- #

    # Fingerprint
    fingerPrint 'unique' ${subf}

  else
    # Absent variant file
    if [ ${VERBOSE} -eq 2 ]; then printf "\n${yellow}WARNING:${normal} VCF file missing for ${subf}. Computation will be skipped.\n" | tee -a 01-Logs/log_${subf}.txt; fi
    continue
  fi

  if [ -f 70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta ]; then
    printf "${subf}\tSuccess\tComputation finished\n" >> progress.txt
  fi

  # ----------------- #
  # ExI. Clean folder #
  # ----------------- #

  CleanFiles ${subf} ${KEEPF} "Final"

done

# Final format

if [ ${VERBOSE} -ge 1 ]; then
  printf "Finish:\n"
  date "+    date: %d/%m/%Y" 
  date "+    time: %H:%M:%S"
  printf "\n"
fi