#!/usr/bin/env python
#----------------------------------------------------------------------------
# Created By    : Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 13/06/2022
# version       : '1.0'
# ---------------------------------------------------------------------------
# """ Script to run samtools in the fingerprint pipeline. """ 
# ---------------------------------------------------------------------------

# Import
import pysam
import argparse
from distutils import util

# Parse arguments
parser = argparse.ArgumentParser(description="Run Samtools in python.")

parser.add_argument('-s', '--sam', help="SAM file (input).", type=str, required=True)
parser.add_argument('-e', '--extract_reads', help="Extract target reads from BAM file (values = true/yes/1 - false/no/0).", type=util.strtobool, required=True)
parser.add_argument('-t', '--threads', help="Threads.", type=str, required=True)

args = parser.parse_args()

# Function
def after_alingment(sam, name, threads, extract):
    # Sam to BAM
    pysam.view("-b", "-@", threads, "-o", name + ".bam", sam,  catch_stdout=False)

    # Sort BAM (Coordinate) for Variant Call
    pysam.sort("-o", name + ".sort.bam", "-O", "bam", name + ".bam", "-@" , threads)

    # Obtain BAM of mapped reads (properly pair)
    pysam.view("-b", "-q", "30", "-f", "0x2", "-@", threads, "-o", name + ".mapped.bam", name + ".sort.bam", catch_stdout=False)

    # Obtain Fastq's
    if extract:
        print("\n\n### Reads recovery ###\n\n")
        pysam.collate(name + ".sort.bam", name + ".collate")
        pysam.fastq("-1", name + "_R1.fastq", "-2", name + "_R2.fastq", "-s", name + "_leftover.fastq", name + ".collate.bam")


# Execution
sam = str(args.sam)
name = sam.replace(".sam", "")
threads = str(args.threads)
extract = args.extract_reads
after_alingment(sam = sam, name = name, threads = threads, extract = extract)

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