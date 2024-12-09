# SyFi

## Pipeline

SyFi is divided into three sequential modules:

Main: This pipeline uses Illumina reads, contigs and a sequence target (e.g., 16S) to obtain the target haplotypes abundances ratio.
Amplicon: This pipeline retrieves amplicon fingerprints from the gene fingerprints from the first module using in silico primers.
Quant: This pipeline takes the results from the first module and quantify the fingerprint abundance from amplicon sequencing data

## Dependencies

The pipeline depends on:

- [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [Bwa-Mem2](https://github.com/bwa-mem2/bwa-mem2)
- [Samtools](http://www.htslib.org/)
- [GATK](https://github.com/broadinstitute/gatk)
- [Picard](https://github.com/broadinstitute/picard)
- [BCFtools](https://samtools.github.io/bcftools/)
- [Whatshap](https://whatshap.readthedocs.io/en/latest/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [Seqtk](https://github.com/lh3/seqtk)
- [Kallisto](https://pachterlab.github.io/kallisto/about)
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [R](https://www.r-project.org/)
- [Salmon](https://combine-lab.github.io/salmon/)
- [Qiime2](https://qiime2.org/)

## Installation

### Build conda environment

__Conda:__

The conda environment is supplemented in the repository. You can create the environment using `mamba env create -f SyFi.yml`. Otherwise, you can try creating your own environment using running the following code:

```
conda env create -n SyFi --file https://data.qiime2.org/distro/core/qiime2-2022.11-py38-linux-conda.yml

conda activate SyFi

mamba install -c conda-forge bioconda::salmon bioconda::spades=3.15.5 bioconda::whatshap=1.7 bioconda::bcftools=1.16 bioconda::bedtools=2.30.0 bioconda::bwa-mem2=2.2.1 bioconda::kallisto=0.48.0 bioconda::picard=2.27.5 bioconda::seqkit=2.3.1 bioconda::seqtk=1.3 bioconda::gatk

```
SyFi's second module (SyFi amplicon) makes use of Qiime2 to *in silico* extract the amplicon sequences from the SyFi-generated haplotypes to build an amplicon fingerprint. For that reason Qiime2 is first installed in the SyFi environment before installing all other remaining packages. 

### Executable software

For the installation of GATK, we downloaded the pre-compile software from their Github site. In the following code we use the *{SOFTWARE_FOLDER_PATH}* variable to define a potential software folder. Please, modify the code with your own folder path.

__JAVA:__

```
cd {SOFTWARE_FOLDER_PATH}
wget "https://download.oracle.com/java/19/latest/jdk-19_linux-x64_bin.tar.gz"
tar -xvzf jdk-19_linux-x64_bin.tar.gz
echo 'export PATH="{SOFTWARE_FOLDER_PATH}/jdk-19.0.1/bin:$PATH"' >> $HOME/.bashrc
```

More information on Java (jdk) [here](https://www.oracle.com/java/technologies/jdk-script-friendly-urls/).

### Download

Download the latest package release. Please, modify the code with your own folder path.

```
cd {SOFTWARE_FOLDER_PATH}
wget https://github.com/adriangeerre/SyFi/releases/download/beta/SyFi.<latest>.zip
unzip SyFi.<latest>.zip
echo 'export PATH="{SOFTWARE_FOLDER_PATH}/SyFi_{version}/:$PATH"' >> $HOME/.bashrc
```

### Usage

```
Usage: ./SyFi.sh <MODULE>

Sequential modules:
  main:      perform fingerprint identification from microbiome data.
  amplicon:  retrieve amplicon fingerprints from gene fingerprints using in silico primers.
  quant:     pseudoalign amplicon sequecing data to firgerprints.

Other:
  help:      display this help message.
  citation:  display citation.
  structure: display folder structure for execution.

```

### Input

SyFi assumes that the genomes and reads are organized in sub-folders inside of the input folder (-i | --input-folder). Each sub-folder should contain the genome (.fasta) and the reads (.fastq.gz). 

```
For example:

  input_folder/
            └── strain_1
                ├── strain_1_R1.fastq.gz
                ├── strain_1_R2.fastq.gz
                └── strain_1.fasta
            └── strain_2
                ├── strain_2_R1.fastq.gz
                ├── strain_2_R2.fastq.gz
                └── strain_2.fasta
            └── ...

```

SyFi loops through the samples of the folder and runs the steps in sequential order. It will run each sample one time and categorize it in *Success*, *Skipped* or *Failed*. Once, it runs over all samples, the option "-f | --force" must be used to re-run the sample through SyFi steps.

### Output

The default (minimum; k=0) output of SyFi consist of:

- 10-Blast/*{strain}*.tsv
- 11-Sequences/*{strain}*/*{strain}*.fasta
- 20-Alignment/*{strain}*/*{strain}*.fasta
- 20-Alignment/*{strain}*/*{strain}*.fastq.gz
- 30-VariantCalling/*{strain}*/variants/*{strain}*.vcf.gz
- 40-Phasing/*{strain}*/*{strain}*_assembly_h *{strain}*.fasta
- 40-Phasing/*{strain}*/*{strain}*_phased.vcf.gz
- 50-haplotypes/*{strain}*/clean\_*{strain}*\_haplotypes.fasta
- 60-Integration/*{strain}*/abundance.tsv
- 60-Integration/*{strain}*/copy_number.tsv
- 60-Integration/*{strain}*/integration.tsv
- 70-Fingerprints/*{strain}*/*{strain}*_all_haplotypes.fasta
- 70-Fingerprints/*{strain}*/seq_h *{number}*.fasta
- 80-Pseudoalignment/*{read_sample}*/quant.sf
- 90-Output/copy_number.tsv
- 90-Output/raw_output_table.txt
- 90-Output/norm_output_table.txt
