#!/bin/bash

# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Gijs Selten, Florian Lamouche and Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 06/04/2022 - 07/11/2022
# version       : '2.0'
# ---------------------------------------------------------------------------
# Pipeline (bash) to perform fingerprint identification from microbiome data.
# ---------------------------------------------------------------------------

function logo() {
	# Colors
	bold=$(tput bold)
	normal=$(tput sgr0)

	# Logo
	echo ""
	echo "${bold}/¯¯¯¯/|¯¯¯| |¯¯¯\ |¯¯¯¯||¯¯¯¯¯¯¯¯¯||¯¯¯¯¯¯¯¯|${normal}"
	echo "${bold}\____\|___| '\_  \|   _||    |\___||_/|  |\_|${normal}"
	echo "${bold}|¯¯¯¯|\  ¯¯\   |     |  |    ¯¯|   |¯\|  |/¯|${normal}"
	echo "${bold}|____|/___/    /____/   |___|¯¯    |________|${normal}"
	echo ""
}

function usage() {
	# Colors
	bold=$(tput bold)
	normal=$(tput sgr0)

	# Logo
	logo help

	# Usage
	echo "Usage: ./$0 <MODULE>"

	# Modules
	printf "\n"
	echo "${bold}Sequential modules:${normal}"
	echo "  main:      perform fingerprint identification from microbiome data."
	echo "  amplicon:  retrieve amplicon fingerprints from gene fingerprints using in silico primers."
	echo "  quant:     pseudoalign amplicon firgerprints to sequencing data."
	printf "\n"
	echo "${bold}Other:${normal}"
	echo "  help:      display this help message."
	echo "  citation:  display citation."
	echo "  structure: display folder structure for execution."
	printf "\n"
}

function folder_structure() {
	# Colors
	bold=$(tput bold)
	normal=$(tput sgr0)

	# Main
	printf "\n"
	echo "${bold}Main module${normal}"
	echo "The pipeline assumes that the genomes and reads are organized in sub-folders inside of the input folder. Each sub-folder should contain the genome (.fasta) and the reads (.fastq.gz)."
	echo "For example:"
	echo ""
	echo "input_folder/
  └── strain_1
    ├── strain_1_R1.fastq.gz
    ├── strain_1_R2.fastq.gz
    └── strain_1.fasta
  └── strain_2
    ├── strain_2_R1.fastq.gz
    ├── strain_2_R2.fastq.gz
    └── strain_2.fasta"
	printf "\n"

	# Quant
	printf "\n"
	echo "${bold}Quant module${normal}"
	echo "The pipeline assumes that the sample reads are organized in sub-folders inside of the input folder (-i | --read_folder). Each sub-folder should contain the reads either in single end or paired end format (.fastq.gz)."
	echo "For example:"
	echo ""
	echo "read_folder/
  └── sample_1
    ├── sample_1_R1.fastq.gz
    └── sample_1_R2.fastq.gz
  └── samplen_2
    ├── sample_2_R1.fastq.gz
    └── sample_2_R2.fastq.gz"
	printf "\n"
}

function citation()
{
	logo help
	printf "\n"
	echo "Thank you for using SyFi. Please, cite:"
	echo ""
	echo "<Include citation>"
	printf "\n"
}

# Module execution functions
function main_module() {
	SyFi_main.sh "$@"
}

function amplicon_module() {
	SyFi_amp.sh "$@"
}

function quant_module() {
	SyFi_quant.sh "$@"
}

# display whether the obligatory files are called
if [[ $# -lt 1  && $1 != "-h" && $1 != "--help" && $1 != "help" && $1 != "citation" && $1 != "structure" ]]; then
	usage
	exit 2
fi

#Get parameters
while [[ "$1" > 0 ]]; do
	case $1 in
		main)
			shift
			main_module "$@"
			exit
			;;
		amplicon)
			shift
			amplicon_module "$@"
			exit
			;;
		quant)
			shift
			quant_module "$@"
			exit
			;;
		structure)
			folder_structure
			exit
			;;
		citation)
			citation
			exit
			;;
		help)
			usage
			exit
			;;
		*)
			echo $1
			echo "ERROR: Invalid argument."
			exit
			;;
	esac
done

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT

function ctrl_c() {
	# Colors
	red=$(tput setaf 1)
	normal=$(tput sgr0)
	printf "\n${red}Execution halted by user.${normal}\n"
	exit
}