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
    echo "Usage: ./$0 -i <INPUT_FOLDER> -s <SEARCH_TARGET> -t <THREADS>"
    printf "\n"
    echo "  -i  | --input_folder     Folder containing input genomes and reads. The software assumes that the folder contains sub-folders for each strain. For more details, execute <pipeline --folder_structure> (REQUIRED)."
    echo "  -s  | --search_target    Genomic region of interest in fasta format, e.g., 16S (REQUIRED)."
    echo "  -l  | --len_deviation    Total base-pairs for the haplotypes to deviate from the target length upstream and downstream (defaut: 100 bp)."
    echo "  -x  | --extension        Reference file extension (default: fasta)."
    echo "  -t  | --threads          Number of threads (default: 1)."
    echo "  -mn | --min_memory       Minimum memory required (default: 4GB)."
    echo "  -mx | --max_memory       Maximum memory required (default: 8GB)."
    echo "  -k  | --keep_files       Keep temporary files [0: none, 1: BAM's, or 2: All] (default: 0)."
    echo "  -f  | --force            Force re-computation of computed samples (default: False)."
    echo "  -h  | --help             Display help."
    echo "  -c  | --citation         Display citation."
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
          └── sub-folder_1
              ├── strain_1_R1.fastq.gz
              ├── strain_1_R2.fastq.gz
              └── strain_1.fasta
          └── sub-folder_2
              ├── strain_2_R1.fastq.gz
              ├── strain_2_R2.fastq.gz
              └── strain_2.fasta"
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
THREADS=1
EXTENSION='fasta'
MIN_MEM=4
MAX_MEM=8
BPDEV=100
FORCE=0
KEEPF=0

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
    -l | --len_deviation)
      shift
      BPDEV=$1
      shift
      ;;
    -x | --extension)
      shift
      EXTENSION=$1
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
    -f | --force)
      shift
      FORCE=1
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
  printf "\nExecution halted by user.\n"
  printf "\nExecution halted by user.\n" >> 01-Logs/log_${subf}.txt
  exit
}

### Checks

# Input folder
if [[ ! -d ${INPUT_FOLDER} ]]; then
  echo "ERROR: Folder ${INPUT_FOLDER} missing."
  exit
fi

# Fasta input
if [[ ! -e ${INPUT_FOLDER}/${subf}/${subf}.${EXTENSION} ]]; then
  echo "ERROR: File ${INPUT_FOLDER}/${subf}/${subf}.${EXTENSION} missing, please use the \"-e/--extension\" argument if the extension is not \"fasta\"."
  exit
fi

# Keep file
if [[ ! "$KEEPF" =~ ^[0-9]+$ || ${KEEPF} -gt 2 || ${KEEPF} -lt 0 ]]; then
  echo "ERROR: Keep temporary files should be a number between 0 and 2 [0: none, 1: BAM's, or 2: All]."
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
  output_dir=$5
  threads=$6
  min_mem=$7
  max_mem=$8

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

  # Index BAM
  # samtools index -@ ${threads} -b ${filtered_bam}

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
  plot-vcfstats -p ${plot_dir} -s ${comp_fn} # (ERROR HERE!)
}


### Execution

#-------------- #
# I. Haplotypes #
#-------------- #

for subf in $(ls ${INPUT_FOLDER}); do

  ## CHECK: Avoid re-computation
  if [[ -f 70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta && ${FORCE} == 0 ]]; then
    printf "\nWARNING: computation finished for ${subf}. To re-run include the -f/--force argument."
    continue
  elif [[ -f 70-Fingerprints/${subf}/${subf}_all_haplotypes.fasta && ${FORCE} == 1 ]]; then
    rm 10-Blast/${subf}.tsv
    # Remove results for sample
    for fld in 11-Sequences 20-Alignment 30-VariantCalling 40-Phasing 50-Haplotypes 60-Kallisto 70-Fingerprints; do
      rm -r ${fld}/${subf}
    done
  fi

  printf "\nSample: $subf\n"

  ## ------------------------------------
  ## Target recovery from contigs (BLAST)
  ## ------------------------------------

  printf "Performing: Mapping; "

  # Create folder
  mkdir -p 10-Blast 11-Sequences/${subf} 01-Logs
  # Blastn
  blastn -query ${INPUT_FOLDER}/${subf}/${subf}.${EXTENSION} -subject ${SEARCH_TARGET} -strand both -outfmt "6 std qseq" > 10-Blast/${subf}.tsv
  printf ">${subf}\n$(cat 10-Blast/${subf}.tsv | head -n 1 | cut -f 13 | sed 's/-//g')" > 11-Sequences/${subf}/${subf}.fasta

  # CHECK: Absent target match (Perhaps, remove small hits!)
  if [ $(cat 11-Sequences/${subf}/${subf}.fasta | wc -l) -lt 1 ]; then
    printf "\nWARNING: No target was found for ${subf}. Computation will be skipped.\n"
    continue
  fi

  ## --------------------------------------------------
  ## Recover pair-end reads for target (BWA & Samtools)
  ## --------------------------------------------------

  # Create folder
  mkdir -p 20-Alignment/${subf}

  printf "Alignment; " 

  # Alignment
  printf "### BWA mapping ###\n\n" > 01-Logs/log_${subf}.txt
  bwa-mem2 index 11-Sequences/${subf}/${subf}.fasta &>> 01-Logs/log_${subf}.txt
  bwa-mem2 mem 11-Sequences/${subf}/${subf}.fasta ${INPUT_FOLDER}/${subf}/${subf}_R1.fastq.gz ${INPUT_FOLDER}/${subf}/${subf}_R2.fastq.gz -t ${THREADS} 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.sam

  # Sam to BAM
  samtools view -b 20-Alignment/${subf}/${subf}.sam -@ ${THREADS} 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.bam
  # Sort BAM (Coordinate) for Variant Call
  samtools sort -o 20-Alignment/${subf}/${subf}.sort.bam -O bam 20-Alignment/${subf}/${subf}.bam -@ ${THREADS} 2>> 01-Logs/log_${subf}.txt
  # Obtain BAM of mapped reads (properly pair)
  samtools view -b -q 30 -f 0x2 20-Alignment/${subf}/${subf}.bam 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.mapped.bam

  printf "Reads recovery; "

  # Obtain Fastq's
  printf "\n\n### Reads recovery ###\n\n" >> 01-Logs/log_${subf}.txt
  samtools collate 20-Alignment/${subf}/${subf}.mapped.bam 20-Alignment/${subf}/${subf}.collate 2>> 01-Logs/log_${subf}.txt
  samtools fastq -1 20-Alignment/${subf}/${subf}_R1.fastq -2 20-Alignment/${subf}/${subf}_R2.fastq -s 20-Alignment/${subf}/${subf}_leftover.fastq 20-Alignment/${subf}/${subf}.collate.bam 2>> 01-Logs/log_${subf}.txt
  gzip -f 20-Alignment/${subf}/${subf}_*.fastq

  ## --------------------------------------------------
  ## Rebuild target from reads (SPAdes, BWA & Samtools)
  ## --------------------------------------------------

  # CHECK: Absent R1/R2 for de novo assembly
  if [[ ! -f 20-Alignment/${subf}/${subf}_R1.fastq.gz || ! -f 20-Alignment/${subf}/${subf}_R2.fastq.gz ]]; then
    printf "\nWARNING: No target reads were recovered for ${subf}. Computation will be skipped.\n"
    continue
  fi

  printf "Contig re-build; "

  # SPAdes (Unambigous nucleotides assembly)
  printf "\n\n### Contig re-build ###\n\n" >> 01-Logs/log_${subf}.txt
	spades.py -1 20-Alignment/${subf}/${subf}_R1.fastq.gz -2 20-Alignment/${subf}/${subf}_R2.fastq.gz -t ${THREADS} -o 20-Alignment/${subf}/spades -m ${MAX_MEM} &>> 01-Logs/log_${subf}.txt
  cp 20-Alignment/${subf}/spades/contigs.fasta 20-Alignment/${subf}/${subf}.fasta

  # Align
  printf "\n\n### Contig re-alignment (BWA) ###\n\n" >> 01-Logs/log_${subf}.txt
  bwa-mem2 index 20-Alignment/${subf}/${subf}.fasta &>> 01-Logs/log_${subf}.txt
  bwa-mem2 mem 20-Alignment/${subf}/${subf}.fasta 20-Alignment/${subf}/${subf}_R1.fastq.gz 20-Alignment/${subf}/${subf}_R2.fastq.gz -t ${THREADS} 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.rebuild.sam
	
  samtools view -b 20-Alignment/${subf}/${subf}.rebuild.sam -@ ${THREADS} 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.rebuild.bam
  samtools sort -o 20-Alignment/${subf}/${subf}.rebuild.sort.bam -O bam 20-Alignment/${subf}/${subf}.rebuild.bam -@ ${THREADS} 2>> 01-Logs/log_${subf}.txt
  samtools view -b -q 30 -f 0x2 20-Alignment/${subf}/${subf}.rebuild.bam 2>> 01-Logs/log_${subf}.txt > 20-Alignment/${subf}/${subf}.rebuild.mapped.bam 

  ## -------------------------------------------
  ## Variant Calling & Phasing (GATK & BCFTools)
  ## -------------------------------------------

  # CHECK: Absent rebuilt BAM
  if [[ ! -f 20-Alignment/${subf}/${subf}.rebuild.sort.bam || ! -f 20-Alignment/${subf}/${subf}.fasta ]]; then
    printf "\nWARNING: No target reads were recovered for ${subf}. Computation will be skipped.\n"
    continue
  fi

  # Create folders
  mkdir -p 30-VariantCalling/${subf}
  mkdir -p 30-VariantCalling/${subf}/mapped_filtered 30-VariantCalling/${subf}/genotyped 30-VariantCalling/${subf}/reference 30-VariantCalling/${subf}/variants

  printf "Variant calling; "

  # Variant call
  printf "\n\n### Variant calling ###\n\n" >> 01-Logs/log_${subf}.txt
  variantCalling 20-Alignment/${subf}/${subf}.rebuild.sort.bam 20-Alignment/${subf}/${subf}.fasta 20-Alignment/${subf}/${subf}.dict 30-VariantCalling/${subf}/mapped_filtered/${subf}.filtered.bam 30-VariantCalling/${subf} ${THREADS} ${MIN_MEM} ${MAX_MEM} &>> 01-Logs/log_${subf}.txt

  # CHECK: Absent variant file
  if [ ! -f 30-VariantCalling/${subf}/variants/${subf}.vcf.gz ]; then
    printf "\nWARNING: VCF file missing for ${subf}. Computation will be skipped.\n"
    continue
  fi

  # CHECK: No variants recovered
  if [ $(zcat 30-VariantCalling/${subf}/variants/${subf}.vcf.gz | grep -v "#" | wc -l) != 0 ]; then
 
    # Create folder
    mkdir -p 40-Phasing/${subf}

    printf "Phasing; "

    # Phasing
    printf "\n\n### Phasing ###\n\n" >> 01-Logs/log_${subf}.txt
    bash 00-Scripts/genome_phase_MOD.sh -s ${subf} -t ${THREADS} -r 20-Alignment/${subf}/${subf}.fasta -v 30-VariantCalling/${subf}/variants/${subf}.vcf.gz -i 20-Alignment -o 40-Phasing/${subf} &>> 01-Logs/log_${subf}.txt

    # CHECK: Absent haplotypes
    if [[ ! -f 40-Phasing/${subf}/${subf}_assembly_h1.fasta || ! -f 40-Phasing/${subf}/${subf}_assembly_h2.fasta ]]; then
      printf "\nWARNING: No haplotypes were recovered for ${subf}. Computation will be skipped.\n"
      continue
    fi

    # INCLUDE HERE HAPLOTYPE LENGTH CLEAN!!! (For small matches!)

    # Concatenate haplotypes and rename headers
    mkdir -p 50-Haplotypes/${subf}
    hnum=$(cat 40-Phasing/${subf}/${subf}_assembly_h*.fasta | sed 's/^> />/g' | grep "^>" | wc -l)

    for n in $(seq ${hnum})
    do
      if [ ${n} -eq "1" ]; then
        cat 40-Phasing/${subf}/${subf}_assembly_h*.fasta | sed 's/_[0-9]_length_.*//g' | sed -z "s|NODE|seq_h${n}|${n}" > tmp
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
      printf "\nWARNING: No haplotypes were recovered for ${subf}. Computation will be skipped.\n"
      continue
    fi

    # Kallisto (Only applied to strains with more that one haplotype)
    if [ $(grep "^>" 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta | wc -l) -gt "1" ]
    then
      # Create folder
      mkdir -p 60-Kallisto

      printf "Abundance ratio; "

      # Kallisto
      printf "\n\n### Kallisto ###\n\n" >> 01-Logs/log_${subf}.txt
      kallisto index -i 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta &>> 01-Logs/log_${subf}.txt
      kallisto quant -i 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx -o 60-Kallisto/${subf} ${INPUT_FOLDER}/${subf}/${subf}_R1.fastq.gz ${INPUT_FOLDER}/${subf}/${subf}_R2.fastq.gz -t ${THREADS} &>> 01-Logs/log_${subf}.txt

      # Filter haplotypes
      cp 60-Kallisto/${subf}/abundance.tsv 60-Kallisto/${subf}/abundance.orig.tsv
      Rscript 00-Scripts/filterHaplotypes.R -i 60-Kallisto/${subf}/abundance.tsv &>> 01-Logs/log_${subf}.txt
    fi
  fi

  # Merge Kallisto output
  # cat 60-Kallisto/*/abundance.tsv > 60-Kallisto/kallisto_output.tsv

  # --------------- #
  # II. Copy Number #
  # --------------- #

  printf "Copy number; "

  # Log
  printf "\n\n### Target Copy Number ###\n\n" >> 01-Logs/log_${subf}.txt

  # Header
  printf "Strain\tGenome_length\tGenome_nbases\tTarget_length\tTarget_nbases\tCopy_number\n" > 60-Kallisto/${subf}/copy_number.tsv

  # Variables
  # Get length of assembly
  lgen=$(cat ${INPUT_FOLDER}/${subf}/${subf}.fasta | grep -v "^>" | tr -d "\n" | wc -c)
  # Get number of bases in assembly reads
  braw=$(zcat ${INPUT_FOLDER}/${subf}/${subf}_R[12].fastq.gz | paste - - - - | cut -f 2 | tr -d "\n" | wc -c)

  if [ $(grep "^>" 20-Alignment/${subf}/${subf}.fasta | wc -l) -gt 1 ]
  then
    # Get length target (16S) - Select sequences ± BPDEV from target
    tl=$(grep -v "^>" ${SEARCH_TARGET} | wc -c)
    max_tl=$((tl+${BPDEV}))
    min_tl=$((tl-${BPDEV}))
    # Filter sequences within length range (array)
    ids=($(grep "^>" 20-Alignment/${subf}/${subf}.fasta | awk -F "_" -v min_tl="$min_tl" '($4 > min_tl)' | awk -F "_" -v max_tl="$max_tl" '($4 < max_tl)' | sed 's/>//g'))

    # Get length of longest recovered target
    l16S=$(for i in ${ids[*]}; do echo $i | cut -d "_" -f 4; done | awk -v max=0 'NR == FNR {if($1 > max) {max = $1}} END {print max}')
    # Number of bases in selected reads
    b16S=$(samtools view 20-Alignment/${subf}/${subf}.rebuild.sort.bam | awk -v var="${ids[*]}" 'BEGIN{split(var, arr); for (i in arr) names[arr[i]]} $3 in names' | cut -f 10 | tr -d "\n" | wc -c)

  else
    # Get length of unique recovered target
    l16S=$(cat 20-Alignment/${subf}/${subf}.fasta | grep -v "^>" | tr -d "\n" | wc -c)
    # Get number of bases of unique recovered target
    b16S=$(zcat 20-Alignment/${subf}/${subf}_R[12].fastq.gz | paste - - - - | cut -f 2 | tr -d "\n" | wc -c)
  fi

  # Compute ratio (round <1 to 1)
  cnum=$(echo "(${b16S}/${l16S}) / (${braw}/${lgen})" | bc -l)
  if (( $(echo "${cnum} < 0.5" | bc -l) )); then
    cnum="1"
  elif (( $(echo "${cnum} > 25.5" | bc -l) )); then
    cnum="25"
  fi

  # File
  printf "${subf}\t${lgen}\t${braw}\t${l16S}\t${b16S}\t${cnum}\n" >> 60-Kallisto/${subf}/copy_number.tsv

  # ---------------- #
  # III. Integration #
  # ---------------- #

  # CHECK: Absent kallisto output
  if [ ! -f 60-Kallisto/${subf}/abundance.tsv ]; then
    printf "\nWARNING: Missing kallisto output for ${subf}. Computation will be skipped.\n"
    continue
  fi

  printf "Integration; "
  # Create folder
  printf "\n\n### Integration ###\n\n" >> 01-Logs/log_${subf}.txt

  # Integrate step I and II
  Rscript 00-Scripts/Integration.R -r 60-Kallisto/${subf}/abundance.tsv -c 60-Kallisto/${subf}/copy_number.tsv -i ${subf} &>> 01-Logs/log_${subf}.txt

  # ---------------- #
  # IV. Fingerprints #
  # ---------------- #

  # CHECK: Absent integration output
  if [ ! -f 60-Kallisto/${subf}/integration.tsv ]; then
    printf "\nWARNING: Missing integration output for ${subf}. Computation will be skipped.\n"
    continue
  fi

  printf "Fingerprint\n"
  printf "\n\n### Fingerprint ###\n\n" >> 01-Logs/log_${subf}.txt

  # Separate fasta into multiples fasta
  seqkit split -i 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta -O 70-Fingerprints/${subf} &>> 01-Logs/log_${subf}.txt

  # Rename haplotypes files
  for i in $(ls -d 70-Fingerprints/${subf}/*); do mv ${i} $(echo ${i} | sed "s/clean_${subf}_haplotypes.part_//g"); done &>> 01-Logs/log_${subf}.txt

  # Keep haplotypes filtered in integration
  keep=($(tail -n +2 60-Kallisto/${subf}/integration.tsv | cut -f 1))
  for i in $(ls 70-Fingerprints/${subf} | sed 's/.fasta//g'); do
    if [[ ! "${keep[*]}" =~ "${i}" ]]
    then
      rm -f 70-Fingerprints/${subf}/${i}.fasta
    fi
  done

  # Concatenate haplotypes
  cat 60-Kallisto/${subf}/integration.tsv | awk '{print $1 "/" $(NF)}' | tail -n +2 > 60-Kallisto/${subf}/tmp.tsv
  for h in $(cat 60-Kallisto/${subf}/tmp.tsv); do
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

  # ----------------- #
  # ExI. Clean folder #
  # ----------------- #

  # Clean folder
  if [ ${KEEPF} == 0 ]; then
    rm 20-Alignment/${subf}/${subf}.bam 20-Alignment/${subf}/${subf}.*.bam
    rm 30-VariantCalling/${subf}/mapped_filtered/*.ba*
    rm -r 40-Phasing/${subf}/mapped
  fi
    if [ ${KEEPF} -lt 2 ]; then
    # 11-Sequences
    rm 11-Sequences/${subf}/${subf}.fasta.*
    # 20-Alignment
    rm -rf 20-Alignment/${subf}/${subf}.sam 20-Alignment/${subf}/${subf}.dict 20-Alignment/${subf}/${subf}.fasta.* 20-Alignment/${subf}/spades 20-Alignment/${subf}/${subf}.rebuild.sam
    # 30-VariantCalling
    rm -rf 30-VariantCalling/${subf}/genotyped 30-VariantCalling/${subf}/mapped_filtered 30-VariantCalling/${subf}/reference  30-VariantCalling/${subf}/variants/tmp_* 30-VariantCalling/${subf}/variants/db_workspace 30-VariantCalling/${subf}/variants/${subf}
    # 40-Phasing
    rm -rf 40-Phasing/${subf}/${subf}_reads.txt 
    # 50-Haplotypes
    rm 50-Haplotypes/${subf}/clean_${subf}_haplotypes.fasta.idx
    # 60-Kallisto
    rm -rf 60-Kallisto/${subf}/abundance.h5 60-Kallisto/${subf}/run_info.json 60-Kallisto/${subf}/tmp.tsv
    # 70-Fingerprints
    rm 70-Fingerprints/${subf}/${subf}_all_haplotypes.tmp.fasta
  fi

done

# Final format
printf "\n"