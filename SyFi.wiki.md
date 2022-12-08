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

**Integration**

I has an issue were the integration failed because the resulting integration file was empty. This happened when I have multiple haplotypes. The issue was a bad while loop in *Integration.R*.

After that, I found out that the ratio could be low (<0.5) causing values to be odd. Therefore, when we encounter a ratio that it is below 1 we automatically upgrade it to 1 because if we are able to compute both the abundance and the copy number is because, at least, there is 1 target copy. However, to inform of this modification, we added the column *adjusted_values* (Yes|No) to determine if the ratio was below 1 (when Yes).

**Modification of progress.txt**

For the moment, the file *progress.txt* is not hidden. That means that the user could modify the file, affecting the execution of the software. This is a double-side sword because users can take advantage of this file to force re-computation of samples but at the same time, if they are not good enough, there modification can cause issues. For example, if you have run a sample through SyFi and later on you remove the corresponding line in *progress.txt*, the sample will be re-computed but the previous files will not be erased. This might cause, issues like files with multiple lines and intermediate errors that will end in Skipped or Failed outcomes.