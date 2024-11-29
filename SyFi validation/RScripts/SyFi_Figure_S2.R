#Required R packages
library("ggplot2") #Version 3.4.2
library("reshape2") #Version 1.4.4
library("ggpubr") #Version 0.6.0
library("ggrepel") #Version 0.9.3
library("tidyheatmaps") #Version 0.1.0
library("grid") #Version 4.4.2
library("patchwork") #Version 1.1.2

#Set working_directory and results directory
working_directory <- ""
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory, "results/", sep="")

###Figure S2 - 16S copy number and haplotypes using direct BLAST =====

#These scripts generate the 16S copy number and haplotype distribution of the 447 AtSC-isolates when blasting directly from the genome assembly instead of using SyFi module 1.
table_16S <- read.table(paste(working_directory, "16S_BLAST/16S_blast.txt", sep = ""), header =F, sep = "\t")
table_16S_hap <- read.table(paste(working_directory, "16S_BLAST/16S_haplotypes.txt", sep =""), header =F,  sep = "\t")
colnames(table_16S) <- c("isolate", "X16S_copy_number")
colnames(table_16S_hap) <- c("isolate", "X16S_haplotypes")

table_16S$Haplotype <- table_16S_hap$X16S_haplotypes[match(table_16S$isolate, table_16S_hap$isolate)]

taxonomy <- read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header =T,  sep = "\t")

table_16S_2 <- table_16S[table_16S$X16S_copy_number != 0,]

plot_16S_CN <- ggplot(table_16S_2, aes(x=X16S_copy_number)) +
  geom_bar() +
  ggtitle("16S copy number - blast") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S copy number") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_CN

plot_16S_hap <- ggplot(table_16S_2, aes(x=Haplotype)) +
  geom_bar() +
  ggtitle("16S haplotypes - blast") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S haplotypes") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_hap

figure <- ggarrange(plot_16S_CN, plot_16S_hap, labels = c("A", "B"),font.label = list(size = 18), ncol =2)

pdf(paste(results.dir,"Figure_S2.pdf", sep=""), width=16, height=6)
print(figure)
dev.off()


