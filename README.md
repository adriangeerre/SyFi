# Pipeline

## Dependencies

The pipeline depends on:

- Blast
- Bwa-Mem2
- Samtools
- GATK
- BCFtools
- SeqKit
- Kallisto
- R

## Build conda environment

```
conda create -n GijsPipeline
conda activate GijsPipeline
conda install -c bioconda blast samtools bwa-mem2 spades bcftools seqkit kallisto
conda install -c conda-forge r-base
```


conda install -c bioconda blast #=2.12.0
conda install -c bioconda samtools=1.13 # Check VariantCalling for correction
conda install -c bioconda bwa-mem2=2.21
conda install -c bioconda spades #=3.15.4
conda install -c bioconda bcftools=1.15.1
conda install -c conda-forge gsl
conda install -c bioconda seqkit #=2.2.0
conda install -c bioconda kallisto #=0.44.0
conda install -c conda-forge r-base