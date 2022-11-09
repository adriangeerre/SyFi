#!/bin/bash
#Author: Petros Skiadas (p.skiadas@uu.nl)
#Created on: October 19 2021
#Modified by: Adrián Gómez Repollés (adrian.gomez@mbg.au.dk)
#Modified on: June 16 2022
#Version: 1.2.0

#parse input options
while getopts "hs:t:r:v:i:o:" OPTION; do
    case $OPTION in
        s)
            sample=${OPTARG};;
        t)
            threads=${OPTARG};;
        r)
            ref_fasta=${OPTARG};;
        v)
            vcf_file=${OPTARG};;
        i)
            illumina_dir=${OPTARG};;
        o)
            output_dir=${OPTARG};;
        h)
            echo -e "${usage}"
            exit 0;;
    esac
done

#check input
if [ -z ${sample} ];then
    echo "You must specify a sample name (-s)">&2
    inputincomplete=true
fi

if [ -z ${threads} ]; then
    echo "You must specify a number of threads (-t)">&2
    inputincomplete=true
fi

if [ -z ${ref_fasta} ]; then
    echo "You must specify the path for the reference fasta file (-r)">&2
    inputincomplete=true
fi

if [ -z ${vcf_file} ]; then
    echo "You must specify the path for the vcf file (-v)">&2
    inputincomplete=true
fi

if [ -z ${illumina_dir} ]; then
    echo "You must specify the folder for the illumina reads per sample (-i)">&2
    inputincomplete=true
fi

if [ -z ${output_dir} ]; then
    echo "You must specify the path for the output directory (-o)">&2
    inputincomplete=true
fi

if [ ${inputincomplete} ]; then
    echo ""
    echo "ERROR: input is incomplete"
    exit 1
fi

#bam
illumina_bam="30-VariantCalling/${sample}/mapped_filtered/${sample}.filtered.bam"
if [ ! -e ${illumina_bam} ]; then
    echo "File ${illumina_bam} not found. Please, provide the file."
    exit 1
fi

#output directories
if [ ! -d ${output_dir} ]; then
    mkdir -p ${output_dir}
fi
#mapping
output_mapped=${output_dir}/mapped/
if [ ! -d ${output_mapped} ]; then
    mkdir -p ${output_mapped}
fi

#expected output from mapping
output_illumina=${output_dir}/mapped/${sample}_illumina.bam
#expected output files from whatshap
output_vcf=${output_dir}/${sample}_phased.vcf
output_gtf=${output_dir}/${sample}_phased.gtf
output_list=${output_dir}/${sample}_reads.txt

#mapping for short reads
bwa_mapping () {
    # Variables
    fasta=${1}
    threads=${2}
    output_illumina=${3}

    # Paired-end reads
    R1="${illumina_dir}/${sample}/${sample}*_R1.fastq.gz"
    R2="${illumina_dir}/${sample}/${sample}*_R2.fastq.gz"
    
    # Alignment
    # bwa-mem2 index ${fasta} # Already done in the variant calling
    # samtools faidx ${fasta} # Already done in the variant calling
    bwa-mem2 mem -M -t ${threads} ${fasta} ${R1} ${R2} | samtools sort -@ ${threads} -m 1G - -o ${output_illumina}
    samtools index -b ${output_illumina}
}


printf "\n--- Mapping ---\n"


##map to reference if no bam files
if ! { [ -f ${illumina_bam} ] && [ -f ${output_illumina} ]; }; then
    bwa_mapping ${ref_fasta} ${threads} ${output_illumina}
fi

printf "\n--- Phasing ---\n"

##phase variants and output vcf file
if [ ! -f ${output_vcf} ]; then
    whatshap phase -o ${output_vcf} --reference ${ref_fasta} ${vcf_file} ${illumina_bam} --ignore-read-groups --sample ${sample} --indels --output-read-list ${output_list} 
fi

printf "\n--- VCF to GTF ---\n"

##visualise
#convert vcf to gtf
if [ ! -f ${output_gtf} ]; then
    whatshap stats --gtf ${output_gtf} ${output_vcf}
fi

printf "\n--- Compress and Index ---\n"


#zip and index vcf
if [ ! -f ${output_vcf}.gz ]; then
    bgzip -c ${output_vcf} > ${output_vcf}.gz
    tabix -p vcf ${output_vcf}.gz
fi

echo "\n--- VCF to Fasta ---\n"

##phased referece
bcftools consensus -H 1 -f ${ref_fasta} ${output_vcf}.gz > ${output_dir}/${sample}_assembly_h1.fasta
bcftools consensus -H 2 -f ${ref_fasta} ${output_vcf}.gz > ${output_dir}/${sample}_assembly_h2.fasta
