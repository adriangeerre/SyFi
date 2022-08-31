# Pipeline

This pipeline uses Illumina reads, contigs and a sequence target (e.g., 16S) to obtain the target haplotypes abundances ratio.

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

## Installation

### Build conda environment

Conda:

The conda environment is supplemented in the repository. You can create the environment using "conda env create -f GijsPipeline.yml". Otherwise, you can try creating your own environment using running the following code:

```
conda create -n GijsPipeline
conda activate GijsPipeline
conda install -c bioconda blast samtools bwa-mem2 spades bcftools seqkit kallisto
conda install -c conda-forge r-base
```

### Executable software

For the installation of GATK, we downloaded the pre-compile software from their Github site. In the following code we use the $SOFTWARE_FOLDER_PATH variable to define a potential software folder. Please, modify the code with your own folder path.

[GATK](https://github.com/broadinstitute/gatk/releases):
```
cd $SOFTWARE_FOLDER_PATH
wget https://github.com/broadinstitute/gatk/archive/refs/tags/4.2.6.1.tar.gz
tar -xvzf 4.2.6.1.tar.gz
echo 'export PATH="$SOFTWARE_FOLDER_PATH/gatk-4.2.6.1:$PATH"' >> $HOME/.bashrc
```
