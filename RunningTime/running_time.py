# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 14/03/2023
# version       : '1.0'
# ---------------------------------------------------------------------------
""" Python script to compute the running time of SyFi. This workflow comes
after running Bowtie2 and obtaining the data stats.""" 
# ---------------------------------------------------------------------------
# Imports 
# ---------------------------------------------------------------------------

import os
import sys
import time
import datetime
import argparse
import subprocess

def params():
    # Parser
    parser = argparse.ArgumentParser(prog='SyFi Running time', description='Python script to measure SyFi\'s running time.')
    
    # Arguments Indiv
    parser.add_argument('-i', '--isolate', dest='iso', action='store', help='Name of bacterial isolate', required=True)
    parser.add_argument('-g', '--genomelen', dest='glen', action='store', help='Genome length', required=True)
    parser.add_argument('-a', '--non_target_reads', dest='nt', action='store', help='Number of non-target reads', required=True)
    parser.add_argument('-b', '--target_reads', dest='t', action='store', help='Number of target reads', required=True)
    parser.add_argument('-l', '--readslen', dest='rlen', action='store', help='Average read length', required=True)
    parser.add_argument('-c', '--coverage', dest='cov', action='store', help='Coverage', required=True)
    parser.add_argument('-r', '--reps', dest='rep', action='store', help='Number of data points to compute (replicates)', required=True)
    parser.add_argument('-t', '--threads', dest='threads', action='store', help='Threads (default: %(default)s)', default=1, type=int)

    # Args
    args = parser.parse_args()

    return args

# Funtions
#---------

def define_proportions(data: dict, bact: str, cov: int):
	# Number of reads required
	nr = round((data[bact][0] * cov) / data[bact][3])

	# Total number of reads
	total = data[bact][1] + data[bact][2]

	# Ratio
	ratio =  nr / total

	# Required vs Original
	if ratio == 1:

		return None # Run SyFi directly

	elif ratio < 1:
		# Subset 
		ntr = round(data[bact][1] * ratio) # All required reads minus 16S reads
		tr = data[bact][2] # All 16S reads

		return [ntr,tr, 'subset']

	elif ratio > 1:
		# Multiply
		ntr = round(data[bact][1] * ratio)
		tr = round(data[bact][2] * ratio)

		# MODIFY THIS TO ACCOUNT FOR times dataset + subset!

		return [ntr,tr, 'multiply']

def exec(data, cov, rep, threads):
	# Variables
	bact = str(list(data.keys())[0])
	name=f"{bact}_{cov}_{rep}_{threads}"

	# Check if already run

	if os.path.exists(f"running_time/{name}.tsv"): sys.exit("Computation was done previously.")

	# Progress
	print(f"Running: {bact} {cov} {rep} {threads}")

	# Create {bact}_{cov}
	if not os.path.isdir(f"temp_files/{name}/{name}"): os.makedirs(f"temp_files/{name}/{name}")

	# Define proportions
	ps = define_proportions(data, bact, cov)

	# Check mode
	if ps[2] == "subset":
		# Define reads
		num_ntr = ps[0] // 2

		# Extract non-target reads (R1 & R2)
		cmd=f"seqtk sample -s100 03-DivideReads/{bact}/{bact}.non_target.1.fastq {num_ntr} > temp_files/{name}/{name}/{bact}.subset.non_target.1.fastq; \
		seqtk sample -s100 03-DivideReads/{bact}/{bact}.non_target.2.fastq {num_ntr} > temp_files/{name}/{name}/{bact}.subset.non_target.2.fastq"
		subprocess.call(cmd, shell=True)

		# Merge reads
		cmd=f"cat 03-DivideReads/{bact}/{bact}_R1.target.fastq temp_files/{name}/{name}/{bact}.subset.non_target.1.fastq > temp_files/{name}/{name}/{name}_R1.fastq; rm temp_files/{name}/{name}/{bact}.subset.non_target.1.fastq"
		subprocess.call(cmd, shell=True)
		cmd=f"cat 03-DivideReads/{bact}/{bact}_R2.target.fastq temp_files/{name}/{name}/{bact}.subset.non_target.2.fastq > temp_files/{name}/{name}/{name}_R2.fastq; rm temp_files/{name}/{name}/{bact}.subset.non_target.2.fastq"
		subprocess.call(cmd, shell=True)

		# Copy fasta to {bact}_{cov}
		cmd=f"cp 00-RTData/{bact}/{bact}.fasta temp_files/{name}/{name}/{name}.fasta"
		subprocess.call(cmd, shell=True)

		# Run SyFi
		today = datetime.date.today().strftime("%d-%m-%Y")
		start = time.time()
		cmd=f"mkdir -p .logs; bash SyFi.sh -i temp_files/{name} -s target.fna -t {threads} --fastq-extension fastq > .logs/{name}_{today}.log"
		subprocess.call(cmd, shell=True)
		end = time.time()
		elapsed = end - start

		# Clean SyFi files
		if os.path.exists("progress.txt"): os.remove("progress.txt")

		# Create output
		if not os.path.isdir("running_time"): os.mkdir("running_time")

		# Save for table
		if not os.path.exists(f"running_time/{name}.tsv"):
			f = open(f"running_time/{name}.tsv","w")
			f.write(f"{bact}\t{cov}\t{rep}\t{start}\t{end}\t{elapsed}\t{threads}\n")
			f.close()

# Execution
#----------

# Parameters
args = params()

# Variables
data = {}
data[args.iso] = [int(args.glen), int(args.nt), int(args.t), int(args.rlen)]

# Call execution
exec(data = data, cov = int(args.cov), rep = args.rep, threads = args.threads)
