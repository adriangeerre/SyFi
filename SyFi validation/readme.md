# SyFi validation

The SyFi validation folder contains all the files to recreate the figures in the manuscript. The folder 'Rscripts' contains the R code to generate and save the main figures 2, 3 and 4, the supplemental figures S1-S7, supplemental table S5, and the figure in Appendix 3. SyFi_full.R contains all scripts together, while the other scripts are Rscripts of the individual figures or tables.

The other files and folders are called into the scripts to generate these. Downloading the entire 'SyFi validation' folder and setting the location of this folder into the working_directory line in the Rscript will allow the user to easily generate all figures. 

What the individual files indicate:

## Main folder

16S_table.txt - An overview of the number of 16S copies and unique haplotypes found in the 447 bacterial isolates found by SyFi. 

AtSC_taxonomy_GTDB.tsv - GTDB (version v2.3.0) database-inferred taxonomic annotation of the 447 bacterial isolates 

No_clusters.txt - The number of distinguishable SynCom isolates from the 447 isolates on the V3-V4 and V5-V7 level when clustering the amplicon sequences on 95, 97, 98, 99, and 100% sequence identity (VSEARCH version 1.2.13) or clustering the SyFi-generated V3-V4 or V5-V7 fingerprints.

metadata.txt - Metadata of the 44 samples used for the SyFi validation. Columns display the host-species on which the SynCom was inoculated (Plant), the inoculum used, the plant compartment that was harvested, the time point at which the sample was harvested, the nutrient condition under which the host was grown, and the biological replicate number. The metadata corresponds to the V3-V4, V5-V7, and the shotgun metagenomics-sequenced dataset as they are all of the same 44 sample.

reads_percentage.txt - the proportion of reads that pseudoaligned/matched to/with the SynCom isolates for each dataset. Datasets include both V3-V4 and V5-V7 datasets in which the SynCom's V3-V4 or V5-V7 sequences are matched with the dataset's ASVs on 95, 97, 98, 99, and 100% sequence identity as well as the SyFi Module 2 data in which reads were pseudoaligned against the V3-V4 and V5-V7 datasets.

strains_with_fingerprint.txt - A list of the isolates for which SyFi was able to fully construct a 16S rRNA fingerprint

## Subfolders

16S_BLAST - Folder with the number of 16S rRNA copies (16S_blast) and 16S rRNA haplotypes (16S_haplotypes.txt) when directly aligning the target 16S rRNA sequence against the bacterial genomes of the 447 SynCom isolates. Which showcases how SyFi is able to find more haplotypes than a direct alignment can.

16S_and_fingerprint_fastas - Folder with the fasta sequences of the studied datasets in the manuscript. These include the directly blasted 16S sequences as well as the SyFi-generated 16S fingerprints of the human gut and plant root (AtSC) bacterial isolates. For the plant root dataset there are also the direct sequences and the fingerprints of the V3-V4 and V5-V7 amplicon regions.

comparison_datasets - Folder with all the V3-V4 and V5-V7 microbiome tables when pseudoaligning the SynCom data reads to the V3-V4 and V5-V7-SyFi-generated fingerprints, as well as when matching the SynCom isolates' V3-V4 and V5-V7 sequences with the SynCom datasets V3-V4 and V5-V7 ASVs with 95, 97, 98, 99, and 100% sequence identity. The datasets that contain 'Salmon' in the name are a comparative dataset in which the SynCom dataset's reads are pseudoaligned to the SynCom isolates' V3-V4 and V5-V7 sequences instead of the SyFi-generated fingerprints. The file AtSC_meta_norm.tsv is the shotgun-metagenomics dataset. How this was constructed can be found here: https://www.biorxiv.org/content/10.1101/2024.08.22.609090v1 .

comparison_datasets_presence_absence - Folder with the same datasets of 'comparison_datasets'. In these datasets, however, each dataset is subsetted to only include the 382 isolates for which SyFi was able to generate a fingerprint. These datasets were used to investigate the presence/absence accuracy of SyFi also shown in Table S5.

complete_genomes - Folder with files that show information regarding a set of 3621 closed genomes downloaded from NCBI. Information is shown regarding the GC content of the 16S rRNA gene and the full genome (16S_GC_genomes_full.txt), the 16S copy number (16S_copy_number.tsv), the number of 16S haplotypes (16S_haplotypes.tsv), and the taxonomy (taxonomy.tsv). 

genome_assembly - Folder with files on the genome assembly quality of the 447 SynCom isolates. The filenames indicate what information is contained, e.g. the GC content in the SyFi-generated 16S rRNA fingerprints (GC_16S_fin.txt), in the genomes (GC_content.txt), the genome size (genome_size.txt), the amount of basepairs in the genomic read files (Genome_total_base_pairs.txt), the L50 (length of the contig that, together with all larger contigs) covers 50% of the genome) (L50.txt), the checkm genome quality report (containing information like completeness, heterogeneity, and contamination), the length of the SyFi-generated 16S rRNA fingerprints (length_16S_fin.txt), and a list of the SynCom isolates (list.txt)

human_set - folder with files that show information on the 16S rRNA copy number and haplotype number of the human gut microbiome bacterial isolate set. Also included is the taxonomy of these strains and a list of dereplicated genomes (genomes more dissimilar > 99.99%) 

minscorefraction - Folder with the V3-V4 and V5-V7 datasets in which the pseudoalignment of reads to the SyFi-generated fingerprints was done at 65, 75, 85, and 95% minscorefraction when using Salmon in SyFi module 2. The AtSC_meta_norm.tsv file is the shotgun-metagenomics dataset of the same samples for comparison and evaluating SyFi's accuracy.

vsearch_clusters - Folder with files that indicates which isolates are clustered together based on their V3-V4 or V5-V7 amplicon sequence at 95, 97, 98, 99, and 100% sequence identity or when clustering at 100% sequence identity on the SyFi-generated fingerprints.

vsearch_clusters_presence_absence - Folder with the same files in the 'vsearch_clusters' folder. In these folders, however, only the amplicon sequences are taken from the 382 SynCom isolates for which SyFi was able to generate a complete 16S rRNA fingerprint. These files were used to investigate the presence/absence accuracy of SyFi also shown in Table S5.




