#!/bin/bash
#Author: Petros Skiadas (p.skiadas@uu.nl)
#Created on: October 19 2021

#tools
# source /opt/conda/anaconda3/etc/profile.d/conda.sh

#parse input options
while getopts "hs:t:r:v:o:" OPTION; do
    case $OPTION in
        s)
            sample=${OPTARG};;
        t)
            threads=${OPTARG};;
        r)
            ref_fasta=${OPTARG};;
        v)
            vcf_file=${OPTARG};;
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

if [ -z ${output_dir} ]; then
    echo "You must specify the path for the output directory (-o)">&2
    inputincomplete=true
fi

if [ ${inputincomplete} ]; then
    echo "input is incomplete"
    exit 1
fi

#########
# ERASE #
#########
SAMTOOLS="/home/au614901/Software/miniconda3/envs/NGSTools/bin/samtools"
BCFTOOLS="/home/au614901/Software/miniconda3/envs/BCFTools/bin/bcftools"

##input reads or bam files
#reads
illumina_dir="00-Data"
#bam
illumina_bam="30-VariantCalling/${sample}/mapped_filtered/${sample}.filtered.bam"

#output directories
if [ ! -d ${output_dir} ]; then
    mkdir ${output_dir}
fi
#mapping
output_mapped=${output_dir}/mapped/
if [ ! -d ${output_mapped} ]; then
    mkdir ${output_mapped}
fi

#expected output from mapping
output_illumina=${output_dir}/mapped/${sample}_illumina.bam
#expected output files from whatshap
output_vcf=${output_dir}/${sample}_phased.vcf
output_gtf=${output_dir}/${sample}_phased.gtf
output_list=${output_dir}/${sample}_reads.txt

#mapping for short reads
bwa_mapping () {
    
    R1="${illumina_dir}/${sample}/${sample}*_R1.fastq.gz"
    R2="${illumina_dir}/${sample}/${sample}*_R2.fastq.gz"
    
    bwa-mem2 index ${1}
    ${SAMTOOLS} faidx ${1}
    bwa-mem2 mem -M -t ${threads} ${1} ${R1} ${R2} | ${SAMTOOLS} sort -@ 8 -m 1G - -o ${output_illumina}
    ${SAMTOOLS} index -b ${output_illumina}
}


##map to reference if no bam files
if ! { [ -f ${illumina_bam} ] && [ -f ${output_illumina} ]; }; then
    bwa_mapping ${ref_fasta}
fi

#activate conda enviroment with whatshap installation
# conda activate base

##phase variants and output vcf file
if [ ! -f ${output_vcf}* ]; then
    whatshap phase -o ${output_vcf} --reference ${ref_fasta} ${vcf_file} ${illumina_bam} \
                   --ignore-read-groups --sample ${sample} --indels --output-read-list ${output_list}
fi

##visualise
#convert vcf to gtf
if [ ! -f ${output_gtf} ]; then
    whatshap stats --gtf ${output_gtf} ${output_vcf}
fi

#zip and index vcf
if [ ! -f ${output_vcf}.gz ]; then
    bgzip -c ${output_vcf} > ${output_vcf}.gz
    tabix -p vcf ${output_vcf}.gz
fi

# conda deactivate

##phased referece
${BCFTOOLS} consensus -H 1 -f ${ref_fasta} ${output_vcf}.gz > ${output_dir}/${sample}_assembly_h1.fasta
${BCFTOOLS} consensus -H 2 -f ${ref_fasta} ${output_vcf}.gz > ${output_dir}/${sample}_assembly_h2.fasta
