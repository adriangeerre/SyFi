# Running SyFi without genomic reads

In the event genomic read files are not present to run SyFi main, the possibility exists to still run SyFi amplicon and SyFi quant. Instead of calling variants to find marker sequence variants, they can be directly extracted from the genomic sequence. These directly-extracted reference sequences, however, may be less accurate as compared to SyFi-generated fingerprints that require the genomic read files.

To skip SyFi main, the following steps can be taken to build references for SyFi quant:

## Preparation 

The following folders and files need to be created by the user:
- A text file that contains a list of the bacterial genome names
- A folder with the bacterial genomes with *.fasta as suffix. These do not require the embedded structure as is shown in SyFi main
- A file with the marker sequence that will be used from the direct alignment (e.g. the 16S rRNA sequence)

```
#Example commands for generation

#List of genomes
nano genome_list.txt

#folder with bacterial genomes
mkdir genomes
while read line ; do mv "$line".fasta genomes/. ; done<genome_list.txt

#sequence for alignment
nano 16S_sequence.fasta
```

## Building the reference without SyFi main

With the SyFi environment activated, we require the 60-Integration folder, the 70-Fingerprints folder, and the Summary.tsv file that SyFi main generates

```
#Make a 70-Fingerprints folder for putting in the references and the SyFi main Integration folder for extracting amplicons
mkdir 70-Fingerprints
mkdir 60-Integration

#Generate header of Summary.tsv file
printf "#Isolate\tInput_folder\tTarget_file\tTarget_length\tRecovered_mean_target_length\tLength_deviation\tRecovered_reads\tNumber_SNPs\tCutoff\tNumber_haplotypes\tCopy_number\tHaplotype_ratio\tModified_output\n" >> Summary.tsv
```

Finally the following code loops through the genomes, extracts all hits from the alignment, and builds a fingerprint from these. Additionally, for each bacterial genome an integration.tsv file (necessary for SyFi amplicon) is build that contains the occurrence of each marker sequence hit (in this case they are all 1). Summary.tsv will contain the 16S copy number by counting all alignment hits (neccesary for SyFi quant).

```
#blast the marker sequence (16S rRNA) to each bacterial genome, extract all hits as fasta's

while read -r line; do
    # Make one folder per genome
    mkdir -p 70-Fingerprints/${line}
    mkdir -p 60-Integration/${line}

    # Blast
    blastn -query genomes/${line}.fasta -subject 16S_query.fasta -strand both -outfmt "6 std qseq" > 70-Fingerprints/${line}/${line}.tsv

    # Prepare integration.tsv file per isolate
    printf "target_id\tlength\ttarget_length\teff_length\test_counts\ttpm\tratio\tratio_round\tcopy_number\thaplotype_divisible\tproportion\tproportion_round\tfinal_output\tper_haplotype\tadjusted_values\n" > 60-Integration/${line}/integration.tsv

    # Concatenate all blast hits together into one reference
    printf ">$line\n" > 70-Fingerprints/${line}/${line}_all_haplotypes.fasta

    count=$(cut -f1 70-Fingerprints/${line}/${line}.tsv | wc -l)

    for size in $(seq 1 ${count}); do
        # Extract the nth line of BLAST result and append to fasta
        head -n ${size} 70-Fingerprints/${line}/${line}.tsv | tail -n1 | cut -f13 | tr -d ' \n' | sed 's/-//g' >> 70-Fingerprints/${line}/${line}_all_haplotypes.fasta
        printf "NNNNNNNNNN" >> 70-Fingerprints/${line}/${line}_all_haplotypes.fasta

        # Save each hit to its own FASTA file
        printf ">SEQ_H$size\n" > 70-Fingerprints/${line}/seq_h${size}.fasta
        head -n ${size} 70-Fingerprints/${line}/${line}.tsv | tail -n1 | cut -f13 | tr -d ' \n' | sed 's/-//g' >> 70-Fingerprints/${line}/seq_h${size}.fasta

	# Add to integration.tsv file for amplicon extraction
	printf "seq_h$size\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t1\tNA\n" >> 60-Integration/${line}/integration.tsv
	
    done

    # Calculate copy number and append to Summary.tsv
    number=$(cut -f1 70-Fingerprints/${line}/${line}.tsv | wc -l)
    printf "$line\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t$number\n" >> Summary.tsv

    #remove intermediate blast
    rm 70-Fingerprints/${line}/${line}.tsv

done < genome_list.txt
```

Using these commands, the folders and files are generated for SyFi amplicon and SyFi quant without using genomic reads

```

#Run SyFi amplicon if necessary, e.g. with the V3-V4 primer sequences
bash SyFi.sh amplicon -i 70-Fingerprints/ -fp CCTACGGGNGGCWGCAG -rp GACTACHVGGGTATCTAATCC -n 300 -x 500
   
#Run SyFi quant with the microbiome reads
bash SyFi.sh quant -i Microbiome_read_folder -f 71-Amplicon/
```
