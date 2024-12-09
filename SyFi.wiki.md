# SyFi

## Table of content:

- [Introduction](#introduction)
- [SyFi steps](#syfi-steps)
- [Usage](#usage)
- [Potential situations](#potential-situations)
- [Explanation of summary and progress text files](#explanation-of-summary-and-progress-text-files)
- [SyFi runtime](#syfi-runtime)
- [SyFi validation dataset](#syfi-validation-dataset)
- [SyFi future implementations](#syfi-future-implementations)
- [Troubleshooting](#troubleshooting)

### Introduction

**SyFi (SynCom Fingerprinting)** is a bioinformatics workflow designed to enhance the identification and quantification of SynCom isolates in complex datasets. To understand root microbiome assembly, dynamics, and functioning, the plant microbiome field has witnessed the emergence of increasingly complex SynCom datasets that mimick the natural root microbiome complexity. These increasingly complex SynComs often consist of microbial isolates with highly similar marker genes, making it difficult for existing bioinformatics methods to distinguish between them accurately. SyFi overcomes this limitation by creating a "fingerprint" of the sequenced marker gene for each SynCom isolate (SyFi main) and using the generated fingerprints as a reference for pseudoalignment of SynCom data reads (SyFi quant), enabling better resolution of isolates with identical or highly similar marker genes.

*SyFi main* 

The first module of SyFi (SyFi main) requires at least three files to generate fingerprints of a SynCom isolate: 1. the genome of the SynCom isolate, 2. the corresponding genomic reads that were used to assemble the genome, and 3. a target fasta sequence of the marker gene (e.g. the *16S rRNA* gene). SyFi aligns the target to the bacterial genome and extracts the sequence with the highest sequence identity score (BLAST v2.13.0). The genomic reads are mapped to the sequence using BWA with default parameters (v2.2.1), and filtered using Samtools (v1.16.1). The target sequence is then reassembled using SPAdes with default parameters (v3.15.5) to ascertain non-ambiguous bases in the marker sequence. The SPAdes-generated target sequence is subjected to length thresholds and trimmed accordingly (using Samtools v1.16.1). 

The cleaned marker sequence is processed through a variant calling pathway (Picard algorithm (v2.27.5) in GATK software (v3.8)), which finds SNPs, insertions and deletions in the marker sequence. When variants are found in the sequence, Whatshap (v1.7) then investigates the co-occurrence of these variants in the marker sequence, which will lead to multiple haplotypes (Mode 1). When multiple haplotypes are find, either by GATK-Whatshap (Mode 1) or already directly by SPAdes (Mode 2), Kallisto (v0.48.0) is employed to pseudoalign the target reads to both haplotypes and estimate the abundance of each haplotype in the genome. Finally, the copy number of each haplotype is calculated by comparing the marker sequence coverage to the total genomic coverage. The haplotype sequences are then concatenated according to their copy number to generate the fingerprint. 

SyFi fingerprints can best be generated on an entire gene and not a fragment (or amplicon) of a gene. Using *in silico* primers, and the second module 'SyFi amplicon', we allow the user to extract the amplicon fingerprint from the original fingerprint generated with SyFi main. 

![SyFi_Fig1](https://github.com/user-attachments/assets/505caf2c-0c23-41b4-a34e-022ff6b39952)

*SyFi quant*

The second module of SyFi (SyFi quant) requires the metagenomic reads from the SynCom dataset and the generated marker fingerprints. SyFi first concatenates all fingerprints into one large fasta file. Using the pseudoalignment tool Salmon (v1.4.0), the metagenomic reads are pseudoaligned to the fingerprints to quantify the SynCom isolates. Finally, the marker copy numbers that were calculated in the first module are used to normalize the Salmon counts. 

The pseudoalignment tool Salmon is used at default settings with the exception of the '--minScoreFraction', which is set at 0.95. This high value allows Salmon to pseudoalign the reads more accurately to highly similar sequences. This is evident when pseudoaligning reads to the fingerprints of two isolates in which both harbor five marker copies, from which one strain has one copy with biological variations that makes it different from the other four copies. Even though the strains are highly similar (sharing four identical marker sequences), the biological variations in the last marker sequence allow distinction between the isolates within a complex SynCom dataset.

![SyFi_Fig2-1](https://github.com/user-attachments/assets/5ff800c5-f5b1-46b6-a0d8-1095870e65e4)

### SyFi steps

*SyFi main*

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
12. Clean files
13. Create Summary

*SyFi amplicon*

1. Import all haplotypes a qiime2 object (qza)
2. Extract amplicon region using *in silico* primers
3. Export extracted amplicon haplotypes into fasta format
4. Fingerprint:
   1. Concatenate all amplicon haplotypes with their defined haplotype ratio (from SyFi main) using 10xN as delimiter
5. Clean files

*SyFi quant*

1. Concatenation
   1. Concatenate fingerprints in the SyFi main output folder (70-Fingerprints or 71-Amplicon) into one fasta file
2. Salmon pseudoalignment
   1. Building Salmon index on fingerprints
   2. Pseudoaligning sample reads to fingerprint index
   3. Processing Salmon output into one microbiome table
3. Copy number normalization
   1. Extracting copy numbers from summary.tsv (SyFi main)
   2. Divide isolate counts by copy number
4. Clean files


### Usage

*SyFi main*

`./SyFi.sh main -i <INPUT_FOLDER> -s <SEARCH_TARGET> -t <THREADS>`

    REQUIRED:
    # Input
    -i  | --input_folder     Folder containing input genomes and reads. The software assumes that the folder contains sub-folders for each strain. For more details, execute <pipeline --folder_structure>.
    -s  | --search_target    Genomic region of interest in fasta format, e.g., 16S.

    OPTIONAL:
    # Haplotype deviation:
    -l  | --len_deviation    Total base-pairs for the haplotypes to deviate from the target length upstream and downstream (default: 100 bp).
    -c  | --cutoff    Maximum ratio deviation between haplotypes per sample. This parameter defined how much can an haplotype deviate from the minimum haplotype ratio (default: 25).
    
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

*SyFi amplicon*

`./SyFi.sh amplicon -i <INPUT_FOLDER> -f <FORWARD_PRIMER> -r <REVERSE_PRIMER> -n <MINIMUM_LENGTH> -x <MAXIMUM_LENGTH>`

    REQUIRED:
    # Input
    -i  | --input_folder     The 70-Fingerprint folder containing SyFi output executed on full length target.
    -f  | --forward_primer    Forward primer for extracting amplicon of interest.
    -r  | --reverse_primer    Reverse primer for extracting amplicon of interest.
    -n  | --minimum_length    Minimum length of the amplicon of interest.
    -x  | --maximum_length    Maximum length of the amplicon of interest.

    # Output options
    -v  | --verbose    Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2).
    
    # Display:
    -h  | --help    Display help.
    --citation    Display citation.
    


*SyFi quant*

`./SyFi.sh quant -i <READ_FOLDER> -s <FINGERPRINT_FOLDER> -t <THREADS>`
    
    REQUIRED:
    # Input
    -i  | --read_folder    Folder containing the sample reads. The software assumes that the folder contains sub-folders for each sample containing the fastq.gz files of either single end reads or paired end read.
    -f  | --fingerprint_folder    SyFi main output folder that contains the fingerprints or SyFi amplicon output folder that contains the amplicon fingerprints. If the user wants to exclude bacterial members, a personalized folder can be created with the SynCom members of interest instead.

    OPTIONAL:
    # Read type
    -r  | --read_type    Paired or single end reads [paired or single] (default: paired).
  
    # Pseudoalignment percentage identity score:
    -m  | --minscorefraction    Percentage identity score that Salmon uses for pseudoaligning metagenomic reads to the SyFi-generated fingerprints (default: 0.95).
 
    # Output options:
    -k  | --keep_files    Keep temporary files [0: Minimum, 1: Salmon output, or 2: All] (default: 0).
    -v  | --verbose    Verbose mode [0: Quiet 1: Samples, or 2: All] (default: 2).
 
    # Display:
    -h  | --help                             Display help.
    --citation                               Display citation.
    --folder_structure                       Display required folder structure.


### Length deviation

Before variants are called in the SPAdes-assembled target sequence of the SynCom isolate, the sequence is subjected to a length deviation threshold. This parameter can be tweaked, so the user can customize this according to the length variability that the marker sequence may have across the bacterial kingdom. 

Reasons to change the length deviation parameter:

-Decreasing the length deviation ascertains more clean and comparable marker gene fingerprints. This way we only include complete target genes and target genes with unclear boundaries. 
-Increasing the length deviation is necessary when the length of the marker sequence can vary substantially between bacterial taxa. Another reason to increase the length deviation is when working with a big proportion of incomplete genome assemblies. SyFi might not be able to find many clean and complete target sequences from incomplete genome assemblies when the length deviation is too strict. A more lenient length deviation may overcome this, though may also increase the risk of incorporating marker sequences deriving from contaminating contigs. 

### Ratio deviation (cutoff)

When calculating the proportion of marker sequence haplotypes in the bacterial genome, a threshold is implemented to avoid unrealistic marker sequence copy numbers that may derive from biological contamination or technical errors. E.g. contaminating reads or PCR errors may lead to the construction of a haplotype that isn't not part of the bacterial genome; a biological or technical artifact. Due to the rarity of these errors among the marker sequence reads, the haplotype proportion may indicate how the actual marker sequence haplotype is a 100-fold more present than the haplotype deriving from these artifacts. A ratio deviation cutoff provides a treshold (e.g. the default is 25) to ensure SyFi doesn't provide a copy number of 100, but instead removes the rare (technical error-derived) haplotype and recalculates the proportion to a biological more meaningful amount.

### --minScoreFraction (Salmon pseudoaligment sequence identity threshold)

Salmon implements a minimum sequence identity threshold of 0.65 as default for pseudoalignment of reads to the index. Since the SyFi-generated fingerprints may vary on only a few nucleotides difference, we increased the default to 0.95. This increase led to greater accuracy in identifying and quantifying SynCom isolates (see manuscript Figure S5). Increasing this parameter even further in SyFi may lead to an enhanced distinction of SynCom isolates that have near-identical fingerprints, though the number of pseudoaligned reads may decrease (see manuscript Figure S6) and therewith the number of identified SynCom isolates. Decreasing the parameter may lead to more pseudoaligned reads, though the accuracy in identifying SynCom isolates may decrease (see manuscript Figure S5). 

### Potential situations

When running SyFi you might encounter different situations for your strains:

1. Success: The workflow was completed and that the fingerprint was recovered.
2. Skipped: The workflow was stopped but it might finish if the running options are modified. This happens when the recovered target is too small in comparison with the provided target (argument `-l/--len-deviation`)
3. Failed: The workflow was stopped by a major issue. This happens when there is no target to be recovered and thus, the workflow can not continue.

The situations are recovered in the file *progress.txt* which it is used to define the strains to run through the workflow. If the strain has a label attached to it, then the argument `-f/--force` must be used to re-run the workflow. If there is no label, the workflow will be run over the strain. This implementation avoids re-computation of results and provides continuity.


### Explanation of summary and progress text files

One of the steps in SyFi main is the copy number and haplotype number calculation, which are employed to build the fingerprint. The Rscript Integration.R is composed of a number of if statements that will lead to the most accurate fingerprint. Many intermediate steps in the SyFi main pipeline, including the intermediate steps taken in that Rscript, are saved in Summary.tsv. Here below is described how to interpret each column in the Summary.tsv file.

*Summary.tsv*

Isolate - Name of the SynCom isolate.

Input_folder - Name of the folder with all SynCom isolates' genomes and genomic read files.

Target_file - name of the target gene sequence file used as input for SyFi (the query).

Target_length - length of the target sequence in Target_file.

Recovered_target_length - length of the target sequence that was retrieved from the SynCom isolate's genome and rebuild with SPAdes. 

Length_deviation - the length deviation parameter in basepairs used in SyFi main (default 100). If the recovered_target_length is larger than 'Target_length + Length_deviation', or smaller than 'Target_length - Length_deviation', SyFi main will be stopped for that SynCom isolate since it was unable to find a proper target.

Recovered_reads - number of reads that map the recovered target (two numbers since there are paired reads).

Number_SNPs - number of variations found in the recovered target using the GATK software.

Cutoff - the maximum amount of target sequence copies that is allowed in the bacterial genome (default 25). If a SynCom isolate with only one haplotype exceeds this number, it is reset to 25. If a SynCom isolate with more than one haplotype exceeds this number, the least frequent occuring haplotype sequence is removed and the copy number is recalculated.

Number_haplotypes - number of unique target sequences found in the SynCom isolate's genome (output of WhatsHap).

Copy_number - the target sequence copy number, which is calculated by dividing the target gene coverage by the genome coverage. 

Haplotype_ratio - the calculated occurrences of each haplotype in the genome. Kallisto is implemented to map the target gene reads to the haplotypes. The number of pseudoaligned reads to each haplotype is divided by the amount of pseudoaligned reads to the haplotype with the fewest amounts of reads. These values are rounded to provide the proportions to which the haplotypes occur in the genome. If the proportion is below 0.5, the haplotype is removed. Similarly, when the total count of all haplotypes (the rounded values) exceed the cutoff value, the haplotype with the fewest amount of reads is removed, and the haplotype proportions are recalculated. This process continues until the sum of the proportions do not exceed the cut-off value or only one haplotype is left.

Modified_output - in the event that only one haplotype remains (after calculating the haplotype ratio) that has a copy number < 0.5, we set this copy number to 1. In that event 'Modified_output' will display 'Yes' instead of 'No', since the output was modified.

*progress.txt*

This text file contains the information for each SynCom isolate whether building a fingerprint was succesful or not (only for SyFi main). In the event it failed, the reason is provided why. 

Be aware that when you later add new SynCom isolates to the input folder for SyFi main to process, SyFi main will use progress.txt to skip SynCom isolates that have already been processed. Removing the lines with specific SynCom isolates will lead to SyFi main recomputing the fingerprints again. 

### SyFi runtime

SyFi's runtime heavily depends on the size of the genomic read files provided and the threads utilized. A SynCom that is composed of 100 isolates will take ~1-2 hours to run when the genomic read files are of the size between 50-100 MB and 10 threads are used. Larger SynComs and larger genomic read files will cause SyFi's runtime to increase, though this might increase the runtime from a couple of hours to half a day or a whole day, and does not lead to months of increased computation.

### SyFi validation dataset

To provide evidence for SyFi's improved accuracy in disentangling complex SynCom datasets, we ran SyFi's Module 1 on a collection of 447 Arabidopsis-derived bacterial genomes (NCBI Project numbers PRJNA1138681, PRJNA1139421 (Genomes), and PRJNA1131834 (Genomic reads)). Subsequently, we used SyFi's Module 2 on a complex SynCom dataset that inoculated these 447 bacterial isolates on Arabidopsis, Barley, and Lotus roots by pseudoalignment of *16S rRNA* V3-V4 and V5-V7 amplicon reads to the SyFi-generated fingerprints (PRJNA1191388) and compared the SynCom isolate counts to a shotgun metagenomics-sequenced dataset of the same samples (PRJNA1131994).

The output of these results can be found in the 'SyFi validation' folder on this GitHub. This folder contains the generated datasets that were compared and the Rscripts that were used to generate the figures used in the SyFi manuscript.

### SyFi future implementations

SyFi is more accurate than current benchmark methodologies to identify and quantify SynCom isolates, however, we are considering the following implementations to enhance SyFi's performance even further:

- Whatshap has been developed to extract both allele sequences of a gene in diploid organisms. Ergo, the phasing tool in SyFi is restricted to find only one other haplotype per sequence. When SPAdes assembles only one  sequence, whatshap can maximally find two haplotypes, while it can find maximally four haplotypes when SPAdes is able to already build two haplotypes, etc. Future implementations of SyFi may include Whatshap polyphase, a tool that can extract multiple haplotypes of a gene as it has been developed for polyploid organisms.

- As indicated in Appendix 3, SyFi may overestimate the marker gene copy number due to a higher GC content in the marker gene sequence as compared to the rest of the genome. Future implementations may include a copy number normalization based on this marker sequence vs genome GC content. 

- It can still occur that SynCom isolates have identical fingerprints. When generating the Salmon index in SyFi module 2 (SyFi quant), a nonredundant index is created. That means that SyFi may indicate the abundances of one strain, though they might also belong to the strain with an identical fingerprint. We are considering adding another plug-in that indicates which isolates have identical fingerprints and are lost in the second module of SyFi. If the user would like this information in the current version of SyFi, they may implement tools like VSEARCH to investigate which fingerprints are identical.

- Kallisto and Salmon are similar pseudoalignment tools with one crucial difference for SyFi. Salmon allows the changing of the minimum sequence identity score (--minScoreFraction), which we have found to be crucial for SyFi's accuracy when quantifying the SynCom isolates (Figure S5 in the manuscript). In the first implementation of SyFi, we still use kallisto to investigate haplotype abundances in the bacterial genomes. We aim to replace kallisto in this first module by Salmon to maintain consistency in used tools.

- Using different target sequences as input for SyFi module 1 (SyFi main) may lead to different outputs. We are considering including a marker sequence database for the most commonly used marker sequences (e.g. *16S rRNA*, *rpoB*, *recA*, and *gyrB*) spanning the entire bacterial kingdom. This database will be used to taxonomically annotate the bacterial genome and use that information to use the marker sequence of the  phylogenetically closest relative as input for SyFi main. 

### Troubleshooting

**Important:** If when running SyFi all the strains shows "WARNING: No target reads were recovered for *{strain}*. Computation will be skipped.", samtools might not be working properly. Please, check if `samtools --version` returns the proper output. If not, check that you run the libcrypto correction mentioned above. Otherwise, try obtaining a working samtools software in your system.

**Important:** If when running SyFi all the strains shows "WARNING: VCF file missing for *{strain}*. Computation will be skipped.", GATK might not be working properly. Please, check if `gatk --version` returns the proper output. If not, install gatk in your system or download the pre-compile version from this repository.

**Important:** The target haplotypes recovered from the illumina reads might differ from the target/s found directly in the reference genome/MAGs, for example, by using Blast. This is because the genome/MAG might have mask the haplotype in a consensus sequence. Thus, from the illumina reads one might recover information from multiple populations. In other words, the consensus sequence "ACGTACGT" might be coming from a) "ACGTACGT" and b) "ACCTACGT" reads in the population. Given that, the "G/C" variant is masked when obtaining the consesus genome/MAG (in this case we have a G), the direct count of 16S haplotypes (in this case 1; a) from the reference could be different that the population 16S haplotypes (in this case 2; a and b).

**Important**: Low quality genomes could be skipped because of: missing or incomplete targets.

**Important**: The file "progress.txt" contains the information for the software to track the computation, the removal of the file will provoke the re-run of all samples even if they were finished (success).
