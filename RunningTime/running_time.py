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
import time
import shutil
import subprocess

# Funtions
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

# Loop
def exec(data, covs, reps, threads):
	# Table
	t = []

	# Loop isolates
	for bact in data.keys():
		for cov in covs:
			for rep in range(1,reps+1):	
				# Create folder
				if not os.path.isdir("tmp"): os.makedirs("tmp")

				# Define proportions
				ps = define_proportions(data, bact, cov)

				# Check mode
				if ps[2] == "subset":
					# Define reads
					num_ntr = ps[0] // 2

					# Extract non-target reads (R1 & R2)
					cmd=f"mkdir -p tmp/{bact}_{cov}; \
					seqtk sample -s100 03-DivideReads/{bact}/{bact}.non_target.1.fastq {num_ntr} > tmp/{bact}_{cov}/{bact}.subset.non_target.1.fastq; \
					seqtk sample -s100 03-DivideReads/{bact}/{bact}.non_target.2.fastq {num_ntr} > tmp/{bact}_{cov}/{bact}.subset.target.2.fastq"
					subprocess.call(cmd, shell=True)

					# Merge reads
					cmd=f"cat 03-DivideReads/{bact}/{bact}.target.1.fastq tmp/{bact}_{cov}/{bact}.subset.target.1.fastq > tmp/{bact}_{cov}/{bact}_{cov}_R1.fastq"
					subprocess.call(cmd, shell=True)
					cmd=f"cat 03-DivideReads/{bact}/{bact}.target.2.fastq tmp/{bact}_{cov}/{bact}.subset.target.2.fastq > tmp/{bact}_{cov}/{bact}_{cov}_R2.fastq"
					subprocess.call(cmd, shell=True)

					# Copy fasta to folder
					cmd=f"cp 00-RTData/{bact}/{bact}.fasta tmp/{bact}_{cov}.fasta"
					subprocess.call(cmd, shell=True)

					# Run SyFi
					start = time.time()
					cmd=f"bash SyFi.sh -i tmp/ -s target.fna -t {threads}"
					subprocess.call(cmd, shell=True)
					end = time.time()
					
					# Save for table
					t.append([bact, cov, rep, start, end])

				# Remove tmp folder
				if os.path.isdir("tmp"): shutil.rmtree("tmp")

	return t


# Variables
data = {'Marinithermus': [2269167,25943791,583,76], 'Marinomonas': [3899940,29132357,238,76], 'Olsenella': [2051896,19937742,1776,76], 'Rhizobium': [4854518,36782323,501,76]}
covs = [30,50,80,100,150,200]
reps=10
threads=4
