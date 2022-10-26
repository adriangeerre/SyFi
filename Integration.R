
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			    : adrian.gomez@mbg.au.dk
# Created Date	: 13/06/2022
# version 		  : '1.0'
# ---------------------------------------------------------------------------
# """ Script to integrate step I (haplotype ratio) and II (copy number). """ 
# ---------------------------------------------------------------------------

# Libraries
# ---------
library(optparse)

# Parameters
# ----------
option_list = list(
  make_option(c("-r", "--haplo_ratio"), action="store", default=NA, type='character', help="Input path for haplotype ratio file (abundance.tsv)."),
  make_option(c("-c", "--copy_number"), action="store", default=NA, type='character', help="Input path for copy number file (copy_number.tsv)."),
  make_option(c("-i", "--isolate"), action="store", default=NA, type='character', help="Isolate name."))
opt = parse_args(OptionParser(option_list=option_list))

# Data
# ----
abun <- read.delim(opt$haplo_ratio, header=T)
cnum <- read.delim(opt$copy_number, header=T)

# Integration
# -----------
# First computation
abun$ratio_round <- round(abun$ratio)
abun$copy_number <- cnum$Copy_number
abun$haplotype_divisible <- sum(abun$ratio_round)
abun$proportion <- abun$copy_number / abun$haplotype_divisible
abun$proportion_round <- round(abun$proportion)
abun$final_output <- unique(abun$proportion_round) * abun$haplotype_divisible
abun$per_haplotype <- unique(abun$proportion_round) * abun$ratio_round

# Check and recompute, if necessary
while (unique(abun$proportion < 0.5) | nrow(abun) == 0) {
  # Remove haplotype with lowest ratio
  min_ratio <- min(abun$ratio)
  abun <- abun[!abun$ratio == min_ratio,]

  # Re-compute
  abun$haplotype_divisible <- sum(abun$ratio_round)
  abun$proportion <- abun$copy_number / abun$haplotype_divisible
  abun$proportion_round <- round(abun$proportion)
  abun$final_output <- unique(abun$proportion_round) * abun$haplotype_divisible
  abun$per_haplotype <- unique(abun$proportion_round) * abun$ratio_round
}

# Save
# ----
write.table(abun, file=paste("60-Kallisto/", opt$i, "/integration.tsv", sep=""), quote=F, col.names=T, row.names=F, sep="\t")
