# SyFi

## Pipeline

This pipeline uses Illumina reads, contigs and a sequence target (e.g., 16S) to obtain the target haplotypes abundances ratio.

## Dependencies

The pipeline depends on:

- [Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [Bwa-Mem2](https://github.com/bwa-mem2/bwa-mem2)
- [Samtools](http://www.htslib.org/)
- [GATK](https://github.com/broadinstitute/gatk)
- [BCFtools](https://samtools.github.io/bcftools/)
- [Whatshap](https://whatshap.readthedocs.io/en/latest/)
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [Seqtk](https://github.com/lh3/seqtk)
- [Kallisto](https://pachterlab.github.io/kallisto/about)
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [R](https://www.r-project.org/)

## Installation

### Build conda environment

__Conda:__

The conda environment is supplemented in the repository. You can create the environment using `conda env create -f SyFi.yml`. Otherwise, you can try creating your own environment using running the following code:

```
conda create -n SyFi -y
conda activate SyFi
conda install -c bioconda blast bwa-mem2 spades bcftools seqkit kallisto whatshap=0.17 python=3.6.13 tabix bedtools seqtk
conda install -c conda-forge r-base r-optparse

# Correct Samtools libcrypto.so error:
cd $ANACONDA_PATH/envs/SyFi/lib
ln -s libcrypto.so.1.1 libcrypto.so.1.0.0
```

Where $ANACONDA_PATH is the path of your own conda/miniconda installation.

### Executable software

For the installation of GATK, we downloaded the pre-compile software from their Github site. In the following code we use the $SOFTWARE_FOLDER_PATH variable to define a potential software folder. Please, modify the code with your own folder path.

__JAVA:__

```
cd $SOFTWARE_FOLDER_PATH
wget "https://download.oracle.com/java/19/latest/jdk-19_linux-x64_bin.tar.gz"
tar -xvzf jdk-19_linux-x64_bin.tar.gz
echo 'export PATH="$SOFTWARE_FOLDER_PATH/jdk-19.0.1/bin:$PATH"' >> $HOME/.bashrc
```

More information on Java (jdk) [here](https://www.oracle.com/java/technologies/jdk-script-friendly-urls/).


__GATK:__
```
cd $SOFTWARE_FOLDER_PATH
wget https://github.com/broadinstitute/gatk/archive/refs/tags/4.2.6.1.tar.gz
tar -xvzf 4.2.6.1.tar.gz
echo 'export PATH="$SOFTWARE_FOLDER_PATH/gatk-4.2.6.1:$PATH"' >> $HOME/.bashrc
```

For more information on how to install GATK refeer to [GATK-Build](https://github.com/broadinstitute/gatk#building).

### Download

Download the latest package release. Please, modify the code with your own folder path.

```
cd $SOFTWARE_FOLDER_PATH
wget ...
tar -xvzf latest_release.tar.gz
echo 'export PATH="$SOFTWARE_FOLDER_PATH/SyFi_<version>/bin:$PATH"' >> $HOME/.bashrc
```

### Usage

```
./SyFi.sh -i <INPUT_FOLDER> -s <SEARCH_TARGET> -t <THREADS>

REQUIRED:
# Input
  -i  | --input_folder     Folder containing input genomes and reads. The software assumes that the folder contains sub-folders for each strain. For more details, execute <pipeline --folder_structure>.
  -s  | --search_target    Genomic region of interest in fasta format, e.g., 16S.

OPTIONAL:
# Haplotype deviation:
  -l  | --len_deviation    Total base-pairs for the haplotypes to deviate from the target length upstream and downstream (defaut: 100 bp).

# Input extension:
  --fasta-extension        Reference file extension (default: fasta).
  --fastq-extension        Illumina reads file extension (default: fastq.gz).

# Computation:
  -t  | --threads          Number of threads (default: 1).
  -mn | --min_memory       Minimum memory required in GB (default: 4GB).
  -mx | --max_memory       Maximum memory required in GB (default: 8GB).

# Output options:
  -k  | --keep_files       Keep temporary files [0: None, 1: BAM's, or 2: All] (default: 0).
  -v  | --verbose          Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2).
  -f  | --force            Force re-computation of computed samples [0: None, 1: All, 2: Skipped, or 3: Failed] (default: 0).

# Display:
  -h  | --help             Display help.
  -c  | --citation         Display citation.
  --folder_structure       Display required folder structure.
```

### Tips and tricks

**Important:** The target haplotypes recovered from the illumina reads might differ from the target/s found directly in the reference genome/MAGs, for example, by using Blast. This is because the genome/MAG might have mask the haplotype in a consensus sequence. Thus, from the illumina reads one might recover information from multiple populations. In other words, the consensus sequence "ACGTACGT" might be coming from a) "ACGTACGT" and b) "ACCTACGT" reads in the population. Given that, the "G/C" variant is masked when obtaining the consesus genome/MAG (in this case we have a G), the direct count of 16S haplotypes (in this case 1; a) from the reference could be different that the population 16S haplotypes (in this case 2; a and b).

**Important**: Low quality genomes could be skipped because of: missing or incomplete targets.
