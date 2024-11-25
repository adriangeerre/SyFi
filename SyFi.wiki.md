# SyFi

Markdown rules:

Url: [Markdown rules!](https://www.markdownguide.org/basic-syntax)

**Bold**
*Italics*
> Quote!

1. Order List
2. Order List
    1. SubList

- UnOrder
- UnOrder
    - UnOrder

`Line of code!`

    block
    of
    code

___

## Table of content:

- [Introduction](#introduction)
- [Usage](#usage)
- [Potential situation](#potential-situations)

### Introduction

...

### Usage

`./SyFi.sh -i <INPUT_FOLDER> -s <SEARCH_TARGET> -t <THREADS>`

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

### Length deviation

Explain here how it works and why it is implemented ...

### Ratio deviation (cutoff)

...

### Potential situations

When running SyFi you might encounter different situations for your strains:

1. Success: The workflow was completed and that the fingerprint was recovered.
2. Skipped: The workflow was stopped but it might finish if the running options are modified. This happens when the recovered target is too small in comparison with the provided target (argument `-l/--len-deviation`)
3. Failed: The workflow was stopped by a major issue. This happens when there is no target to be recovered and thus, the workflow can not continue.

The situations are recovered in the file *progress.txt* which it is used to define the strains to run through the workflow. If the strain has a label attached to it, then the argument `-f/--force` must be used to re-run the workflow. If there is no label, the workflow will be run over the strain. This implementation avoids re-computation of results and provides continuity.


### Issues found (for us):

**Integration 1**

I has an issue were the integration failed because the resulting integration file was empty. This happened when I have multiple haplotypes. The issue was a bad while loop in *Integration.R*.

After that, I found out that the ratio could be low (<0.5) causing values to be odd. Therefore, when we encounter a ratio that it is below 1 we automatically upgrade it to 1 because if we are able to compute both the abundance and the copy number is because, at least, there is 1 target copy. However, to inform of this modification, we added the column *adjusted_values* (Yes|No) to determine if the ratio was below 1 (when Yes).

**Integration 2**

When multiple haplotypes are present, the script integration.R works in the following way:
- Using the abundance (kallisto) and the copy number (syfi), we compute the proportion per haplotype by dividing the copy number by the haplotype divisible (or the rounded abundance ratio).
- If the proportion is below 0.5 we remove the haplotype with the smallest abundance ratio and recompute the proportion until it is equal or above 0.5 or there is only one haplotype left.
- In the case the remaining haplotype is below 0.5, we round the proportion to 1 and we define the column "adjusted values" as "Yes". In the column the value "No" means that the values are the original results from the calculations and they have not been adjusted by us.
- If the proportion is above 0.5, we keep the haplotypes and no further processing is performed.

**Modification of progress.txt**

For the moment, the file *progress.txt* is not hidden. That means that the user could modify the file, affecting the execution of the software. This is a double-side sword because users can take advantage of this file to force re-computation of samples but at the same time, if they are not good enough, there modification can cause issues. For example, if you have run a sample through SyFi and later on you remove the corresponding line in *progress.txt*, the sample will be re-computed but the previous files will not be erased. This might cause, issues like files with multiple lines and intermediate errors that will end in Skipped or Failed outcomes.

**SPAdes does not produce a contigs fasta**

I have seen in P1_A8 and P2_G4 that the software crashes (seqtk) because SPAdes does not produce a fasta file. This seems to be caused by a low number of recovered reads (288 and 44, repectively). I have add a check and a label in the progress file ("SPAdes did not assembly any target sequence").

### SyFi steps:

1. Mapping:
    1. Run blastn using the reference fasta file as the query and the target as the subject
    2. Extract the largest hit (or one out of many) from the blast output
2. Alignment:
    1. Index target fasta file (Step 1.2)
    2. Map input paired-end reads against index
    3. Subset high quality and proper paired reads of the alignment output
3. Read recovery:
    1. Extract paired-end reads from subset (Step 2.1)
4. Contig re-build:
    1. _De novo_ assembly of extracted paired-end reads
    2. Minimum size select the _de novo_ contigs
    3. Run Blast against the target to determine the sides (left and right) flanking regions
    4. Define size of blast hit without flanking regions
    5. Select the largest target recovered
5. Contig re-alignment:
    1. Align reads against largest recovered target
6. Variant Calling:
    1. Remove alignment duplicates (Step 5.1)
    2. Haplotype caller
    3. Join genotyping
    4. Define if sample has 1 or more haplotypes (2 ways: unique or multiple)
7. Phasing
    1. Perform the phasing using the reference fasta, variants file and alignment file
    2. Call the consensus haplotypes of allele 1 and 2 individually
8. Abundance ratio:
    1. Index the clean haplotypes (Step 7.X)
    2. Quantify the abundance of each haplotypes using the input peared-end reads
    3. Filter the abundance using the cutoff
9. Target Copy Number:
    1. Define length of assembly
    2. Define number of bases in reads
    3. Define length of largest recovered target
    4. Define number of bases in reads within the recovered target
    5. Compute target's copy number
10. Integration:
    1. Merge haplotype abundance and target copy number
    2. Compute proportion and haplotype divisible (Modes: unique or multiple)
    3. Recursively, remove the haplotype with the lowest ratio if proportion is below 0.5 and if there is more than one haplotype left
    4. Adjust values when there is one haplotype and proportion is below 0.5
11. Fingerprint:
    1. Concatenate all haplotypes with their defined haplotype ratio using 10xN as delimiter
12. Clean Files
13. Create Summary




### Troubleshooting

**Important:** If when running SyFi all the strains shows "WARNING: No target reads were recovered for *{strain}*. Computation will be skipped.", samtools might not be working properly. Please, check if `samtools --version` returns the proper output. If not, check that you run the libcrypto correction mentioned above. Otherwise, try obtaining a working samtools software in your system.

**Important:** If when running SyFi all the strains shows "WARNING: VCF file missing for *{strain}*. Computation will be skipped.", GATK might not be working properly. Please, check if `gatk --version` returns the proper output. If not, install gatk in your system or download the pre-compile version from this repository.

**Important:** The target haplotypes recovered from the illumina reads might differ from the target/s found directly in the reference genome/MAGs, for example, by using Blast. This is because the genome/MAG might have mask the haplotype in a consensus sequence. Thus, from the illumina reads one might recover information from multiple populations. In other words, the consensus sequence "ACGTACGT" might be coming from a) "ACGTACGT" and b) "ACCTACGT" reads in the population. Given that, the "G/C" variant is masked when obtaining the consesus genome/MAG (in this case we have a G), the direct count of 16S haplotypes (in this case 1; a) from the reference could be different that the population 16S haplotypes (in this case 2; a and b).

**Important**: Low quality genomes could be skipped because of: missing or incomplete targets.

**Important**: The file "progress.txt" contains the information for the software to track the computation, the removal of the file will provoke the re-run of all samples even if they were finished (success).
