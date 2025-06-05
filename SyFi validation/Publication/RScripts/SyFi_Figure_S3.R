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

###Figure S3 - 16S haplotypes using direct BLAST =====

#Human
table_16S_human <- read.table(paste(working_directory, "16S_BLAST/16S_human_blast.txt", sep = ""), header =F, sep = "\t")
colnames(table_16S_human) <- c("isolate","Haplotype")

table_16S_human_2 <- table_16S_human[table_16S_human$Haplotype != 0,]

plot_16S_hap_human <- ggplot(table_16S_human_2, aes(x=Haplotype)) +
  geom_bar() +
  ggtitle("16S rRNA haplotypes - Human gut bacteria - BLAST") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S rRNA haplotypes") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_hap_human

#These scripts generate the 16S copy number and haplotype distribution of the 447 AtSC-isolates when blasting directly from the genome assembly instead of using SyFi module 1.
table_16S_plant <- read.table(paste(working_directory, "16S_BLAST/16S_plant_blast.txt", sep = ""), header =F, sep = "\t")
colnames(table_16S_plant) <- c("isolate", "Haplotype")

table_16S_plant_2 <- table_16S_plant[table_16S_plant$Haplotype != 0,]

plot_16S_hap_plant <- ggplot(table_16S_plant_2, aes(x=Haplotype)) +
  geom_bar() +
  ggtitle("16S rRNA haplotypes - Plant root - BLAST") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S rRNA haplotypes") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_hap_plant

figure <- ggarrange(plot_16S_hap_human,plot_16S_hap_plant, labels = c("A", "B"),font.label = list(size = 18), ncol =2, nrow =1)

pdf(paste(results.dir,"Figure_S3.pdf", sep=""), width=18, height=6)
print(figure)
dev.off()


