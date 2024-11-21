#!/bin/bash

# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By    : Gijs Selten, Florian Lamouche and Adrián Gómez Repollés
# Email         : adrian.gomez@mbg.au.dk
# Created Date  : 06/04/2022 - 07/11/2022
# version       : '1.2'
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
if [[ $# -lt 1  && $1 != "-h" && $1 != "--help" && $1 != "help" && $1 != "--citation" ]]; then
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