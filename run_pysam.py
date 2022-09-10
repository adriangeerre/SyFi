#!/usr/bin/env python

import argparse, os

# Parse arguments
parser = argparse.ArgumentParser(description="Run Samtools in python.")

parser.add_argument('-s', '--sam', help="SAM file (input).", metavar="STR", type=str, required=True)
parser.add_argument('-t', '--threads', help="Threads.", metavar="STR", type=str, required=True)

parser.parse_args()


def after_alingment(sam, name, threads):
    # Sam to BAM
    pysam.view("-b", "-@", str(args.threads), "-o", name + ".bam", sam)

    # Sort BAM (Coordinate) for Variant Call
    pysam.sort("-o", name + ".sort.bam", "-O", "bam", name + ".bam", "-@" , threads)

    # Obtain BAM of mapped reads (properly pair)
    pysam.view("-b", "-q", "30", "-f", "0x2", "-@", threads, "-o", name + ".mapped.bam", name + ".sort.bam", catch_stdout=False)

    print("Reads recovery; ")

    # Obtain Fastq's
    print("\n\n### Reads recovery ###\n\n")
    pysam.collate(name + ".sort.bam", name + ".collate")
    pysam.fastq("-1", name + "_R1.fastq", "-2", name + "_R2.fastq", "-s", name + "_leftover.fastq", name + ".collate.bam")


# Execution
sam = str(args.sam)
name = sam.replace(".sam", "")
after_alingment(sam = sam, name = name, threads = threads)

# -------------------------------------- #

# PySam

# Properly pair reads
#pysam.view("-b", "-q", "30", "-f", "0x2", "-@", threads, "-o", "strain.bam", "strain.sam", catch_stdout=False)

# SamToBam and Sort
#pysam.sort("-o", "strain.sort.bam", "strain.bam")

# Collate
#pysam.collate("strain.sort.bam", "strain.collate")

# Extract reads
#pysam.fastq("-1", "strain_R1.fastq", "-2", "strain_R2.fastq", "-s", "strain_leftover.fastq", "strain.collate.bam") # leftover are empty