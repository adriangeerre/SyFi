# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			    : adrian.gomez@mbg.au.dk
# Created Date	: 13/06/2022
# version 		  : '1.0'
# ---------------------------------------------------------------------------
# """ Script to integrate step I (haplotype ratio) and II (copy number) if 
#  there are multiple haplotypes. If not, it will produce the integration 
#  with the copy number only.""" 
# ---------------------------------------------------------------------------

# Libraries
# ---------
library(optparse)

# Parameters
# ----------
option_list = list(
  make_option(c("-r", "--haplo_ratio"), action="store", default=NA, type='character', help="Input path for haplotype ratio file (abundance.tsv)."),
  make_option(c("-c", "--copy_number"), action="store", default=NA, type='character', help="Input path for copy number file (copy_number.tsv)."),
  make_option(c("-i", "--isolate"), action="store", default=NA, type='character', help="Isolate name."),
  make_option(c("-m", "--mode"), action="store", default=NA, type='character', help="Mode."))
opt = parse_args(OptionParser(option_list=option_list))

if (opt$mode == "unique") {
  # Data
  # ----
  cnum <- read.delim(opt$copy_number, header=T)
  
  # Integration
  # -----------
  inte <- data.frame(target_id = "seq_h1", total_length=cnum$Target_length, target_length=cnum$Target_length, eff_length="-", est_counts="-", tpm="-", ratio=1, ratio_round=1, copy_number=cnum$Copy_number, haplotype_divisible=1, proportion=cnum$Copy_number, proportion_round=round(as.numeric(cnum$Copy_number)), final_output=round(as.numeric(cnum$Copy_number)), per_haplotype=round(as.numeric(cnum$Copy_number)), adjusted_values = "No")
  
  # Save
  # ----
  write.table(inte, file=paste("60-Integration/", opt$i, "/integration.tsv", sep=""), quote=F, col.names=T, row.names=F, sep="\t")
  
  
} else if (opt$mode == "multiple") {
  # Data
  # ----
  abun <- read.delim(opt$haplo_ratio, header=T)
  cnum <- read.delim(opt$copy_number, header=T)
  
  # Integration
  # -----------
  # First computation
  abun$target_length <- cnum$Target_length
  abun <- abun[, c("target_id", "length", "target_length", "eff_length", "est_counts", "tpm", "ratio")]
  abun$ratio_round <- round(abun$ratio)
  abun$copy_number <- cnum$Copy_number
  abun$haplotype_divisible <- sum(abun$ratio_round)
  abun$proportion <- abun$copy_number / abun$haplotype_divisible
  abun$proportion_round <- round(abun$proportion)
  abun$final_output <- unique(abun$proportion_round) * abun$haplotype_divisible
  abun$per_haplotype <- unique(abun$proportion_round) * abun$ratio_round
  abun$adjusted_values <- "No"
  
  # Check and recompute, if necessary
  enter_while = TRUE
  while (unique(abun$proportion) < 0.5) {
    enter_while = FALSE
    # Remove haplotype with lowest ratio
    min_ratio <- min(abun$ratio)
    if (nrow(abun) <= 1) {
      break
    } else {
      abun <- abun[!abun$ratio == min_ratio,]
    }
    
    # Re-compute
    abun$haplotype_divisible <- sum(abun$ratio_round)
    abun$proportion <- abun$copy_number / abun$haplotype_divisible
    abun$proportion_round <- round(abun$proportion)
    abun$final_output <- unique(abun$proportion_round) * abun$haplotype_divisible
    abun$per_haplotype <- unique(abun$proportion_round) * abun$ratio_round
    abun$adjusted_values <- "No"
  }
  
  if (enter_while == FALSE & unique(abun$proportion) < 0.5) {
    # Select maximum ratio haplotype
    max_ratio <- max(abun$ratio)
    abun <- abun[abun$ratio == max_ratio,]
    
    # Transform ratio to 1
    abun$ratio <- 1
    abun$ratio_round <- 1
    
    # Transform copy number to 1 if below, otherwise leave it
    if (abun$copy_number < 1) {
      abun$copy_number <- 1
    }
    
    # Recompute values
    abun$haplotype_divisible <- sum(abun$ratio_round)
    abun$proportion <- abun$copy_number / abun$haplotype_divisible
    abun$proportion_round <- round(abun$proportion)
    abun$final_output <- unique(abun$proportion_round) * abun$haplotype_divisible
    abun$per_haplotype <- unique(abun$proportion_round) * abun$ratio_round
    abun$adjusted_values <- "Yes"
  }
  
  # Check and save
  if (unique(abun$proportion) >= 0.5) {
    # Save
    #-----
    write.table(abun, file=paste("60-Integration/", opt$i, "/integration.tsv", sep=""), quote=F, col.names=T, row.names=F, sep="\t")
  }
}
