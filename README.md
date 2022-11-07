# Syfi

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
- [Kallisto](https://pachterlab.github.io/kallisto/about)
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [R](https://www.r-project.org/)

## Installation

### Build conda environment

__Conda:__

The conda environment is supplemented in the repository. You can create the environment using `conda env create -f GijsPipeline.yml`. Otherwise, you can try creating your own environment using running the following code:

```
conda create -n SyFi -y
conda activate SyFi
conda install -c bioconda blast bwa-mem2 spades bcftools seqkit kallisto whatshap=0.17 python=3.6.13 tabix
conda install -c conda-forge r-base r-optparse

# Correct Samtools libcrypto.so error: (Whatshap still incompatible!)
# cd anaconda/envs/GijsPipeline/lib
# ln -s libcrypto.so.1.1 libcrypto.so.1.0.0
```

### Executable software

For the installation of GATK, we downloaded the pre-compile software from their Github site. In the following code we use the $SOFTWARE_FOLDER_PATH variable to define a potential software folder. Please, modify the code with your own folder path.

__GATK:__
```
cd $SOFTWARE_FOLDER_PATH
wget https://github.com/broadinstitute/gatk/archive/refs/tags/4.2.6.1.tar.gz
tar -xvzf 4.2.6.1.tar.gz
echo 'export PATH="$SOFTWARE_FOLDER_PATH/gatk-4.2.6.1:$PATH"' >> $HOME/.bashrc
```

__Samtools:__

Unfortunately, the softwares Samtools and Whatshap can not be installed in the same conda environment. After testing long, we found major incompatibilities between the packages among the different available version. Thus, we decided to install samtools and htslib manually because whatshap only provide the conda installation.

The compilation was tested in an Ubuntu 20.04 machine. You might require sudo access to run some steps of the installation. The htslib include also tabix and bgzip, you can decide to include them or not in the conda environment.

```
# Download source
cd $SOFTWARE_FOLDER_PATH
wget https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2

# Decompress
bzip2 samtools-1.16.tar.bz2
tar -xf samtools-1.16.tar
bzip2 htslib-1.16.tar.bz2
tar -xf htslib-1.16.tar

# Required libraries
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install libncurses5-dev libncursesw5-dev zlib1g-dev libbz2-dev liblzma-dev

# Build source (extracted from http://www.htslib.org/download/)
cd samtools-1.16
./configure --prefix=$SOFTWARE_FOLDER_PATH/Samtools-1.16
make
sudo make install

cd htslib-1.16
./configure --prefix=$SOFTWARE_FOLDER_PATH/HtsLib-1.16
make
sudo make install

# Once install correct Owner:Group of the folders

# Add to path
echo 'export PATH="$SOFTWARE_FOLDER_PATH/Samtools-1.16/bin:$PATH"' >> $HOME/.bashrc
echo 'export PATH="$SOFTWARE_FOLDER_PATH/HtsLib-1.16/bin:$PATH"' >> $HOME/.bashrc

# Test instalation in a new terminal
samtools view
```

### Download

Download the latest package release. Please, modify the code with your own folder path.

```
cd $SOFTWARE_FOLDER_PATH
wget ...
tar -xvzf latest_release.tar.gz
```

### Usage

```
./pipeline.sh -i <INPUT_FOLDER> -s <SEARCH_TARGET> -t <THREADS>
```

### Tips and tricks

**Important:** The target haplotypes recovered from the illumina reads might differ from the target/s found directly in the reference genome/MAGs, for example by using Blast. This is because the reference/MAG might have mask the haplotype in a consensus sequence. Thus, from the illumina reads one might recover information from multiple populations.

The option "-f" should always placed last!
Low quality genomes could be skipped because of: missing or incomplete targets.
