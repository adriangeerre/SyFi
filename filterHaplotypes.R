# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  	: Adrián Gómez Repollés
# Email			    : adrian.gomez@mbg.au.dk
# Created Date	: 13/06/2022
# version 		  : '1.0'
# ---------------------------------------------------------------------------
# """ Script to compute the haplotype ratio from the Kallisto output. The ra-
# tio will be used to filter out unwanted haplotypes.""" 
# ---------------------------------------------------------------------------

# Libraries
# ---------
library(optparse)

# Parameters
# ----------
option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="Input path with kallisto output file."),
  )
opt = parse_args(OptionParser(option_list=option_list))

# Data
# ----
df <- read.delim(opt$i, header=T)

# Transform
# --------
cutoff = 25

while(TRUE) {
  minimum <- min(df$est_counts)
  if (minimum == 0) {
    df <- df[df$est_counts != minimum,]
  }
  v <- df$est_counts / minimum
  if (max(v) > cutoff) {
    df <- df[df$est_counts != minimum,]
  } else {
    # Ratio
    df$ratio <- df$est_counts / minimum
    break
  }
}

# Save
# ----
write.table(df, file=opt$i, quote=F, col.names=T, row.names=F, sep="\t")
