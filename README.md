# SyFi

## Pipeline

This pipeline uses Illumina reads, contigs and a sequence target (e.g., 16S) to obtain the target haplotypes abundances ratio.

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

Where *{ANACONDA_PATH}* is the path of your **own anaconda/miniconda** installation.

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
wget https://github.com/adriangeerre/SyFi/releases/download/beta/SyFi_beta.zip
unzip SyFi_beta.zip
echo 'export PATH="{SOFTWARE_FOLDER_PATH}/SyFi_{version}/:$PATH"' >> $HOME/.bashrc
```

### Usage

```
SyFi.sh -i <INPUT_FOLDER> -s <SEARCH_TARGET> -t <THREADS>

REQUIRED:
# Input
  -i | --input-folder     Folder containing input genomes and reads. The software assumes that the folder contains sub-folders for each strain. For more details, execute <pipeline --folder-structure>.
  -s | --search-target    Genomic region of interest in fasta format, e.g., 16S.

OPTIONAL:
# Haplotype deviation:
  -l | --len-deviation    Total base-pairs for the haplotypes to deviate from the target length upstream and downstream (defaut: 100 bp).
  -c | --cutoff            Maximum ratio deviation between haplotypes per sample. This parameter defined how much can an haplotype deviate from the minimum haplotype ratio (default: 25).

# Input extension:
  --fasta-extension        Reference file extension (default: fasta).
  --fastq-extension        Illumina reads file extension (default: fastq.gz).

# Computation:
  -t | --threads          Number of threads (default: 1).
  -m | --memory           Memory in GBs (default: 8GB).

# Output options:
  -k | --keep-files       Keep temporary files [0: Minimum, 1: BAM's, or 2: All] (default: 0).
  -v | --verbose          Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2).
  -f | --force            Force re-computation of computed samples [0: None, 1: All, 2: Skipped, or 3: Failed] (default: 0).

# Display:
  -h | --help             Display help.
  --citation               Display citation.
  --folder-structure       Display required folder structure.
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
