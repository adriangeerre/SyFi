#Required R packages
library("ggplot2") #Version 3.4.2
library("reshape2") #Version 1.4.4
library("ggpubr") #Version 0.6.0
library("ggrepel") #Version 0.9.3
library("tidyheatmaps") #Version 0.1.0
library("grid") #Version 4.4.2
library("patchwork") #Version 1.1.2

#Set working_directory and results directory
working_directory <- "/home/gijs_selten/SyFi_scripts_2024/"
dir.create(paste(working_directory, "results", sep = ""))
results.dir <- paste(working_directory, "results/", sep="")

###Figure S4 - 16S copy number and haplotypes in full dataset of complete genomes =====

#This script generates the figures showing the 16S copy number and haplotype distribution from closed bacterial genomes
table_16S <- read.table(paste(working_directory, "complete_genomes/16S_copy_number.tsv", sep = ""), header =T, sep = "\t")
table_16S_hap <- read.table(paste(working_directory,"complete_genomes/16S_haplotypes.tsv", sep =""), header =T,  sep = "\t")
table_16S$Haplotype <- table_16S_hap$X16S_haplotypes[match(table_16S$isolate, table_16S_hap$isolate)]

taxonomy <- read.table(paste(working_directory, "complete_genomes/taxonomy.tsv", sep =""), header =T,  sep = "\t")

table_16S_2 <- table_16S[table_16S$X16S_copy_number != 0,]

#Family plot
table_16S_2$Family <- taxonomy$Family[match(table_16S_2$isolate, taxonomy$Isolate)]
table_16S_2 <- na.omit(table_16S_2)

table_per_family <- as.data.frame(unique(table_16S_2$Family))
table_per_family <- table_per_family[!table_per_family$`unique(table_16S_2$Family)` %in% c("Bacteroidales", "Thermodesulfobacteriota", "SKNC01", "GCA-2696645", "UBA5962", "DSM-45221", "Burkholderiales"),]
table_per_family <- data.frame(table_per_family)
colnames(table_per_family) <- "Family"
table_per_family$average_CN <- "NA"
family_vector <- as.vector(table_per_family$Family)
table_per_family$size_family <- "NA"


for (family in family_vector) {
  table_sub <- table_16S_2[table_16S_2$Family == paste(family),]
  average <- unique(ave(table_sub$X16S_copy_number))
  table_per_family$average_CN[table_per_family$Family == paste(family)] <- average
  length <- length(table_sub$Family)
  table_per_family$size_family[table_per_family$Family == paste(family)] <- length
}

table_per_family$length_family <- paste(table_per_family$Family," (n = ", table_per_family$size_family,")",sep="")

table_per_family_2 <- table_per_family[order(as.numeric(table_per_family$average_CN),decreasing=T),]
ordered_vector <- as.vector(table_per_family_2$length_family)

table_16S_2$length_family <- table_per_family$length_family[match(table_16S_2$Family,table_per_family$Family)]
table_16S_2$length_family <- factor(table_16S_2$length_family, levels = ordered_vector)
table_4 <- table_16S_2[c("isolate", "X16S_copy_number", "Haplotype", "Family", "length_family")]

table_5 <- melt(table_4)
table_5 <- na.omit(table_5)
table_6 <- table_5[table_5$variable == "X16S_copy_number",]

plot_family <- ggplot(table_6, aes(x=length_family, y=value, color=variable)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90,size = 6, hjust = 1)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  ylab("16S copy number") + 
  xlab("Bacterial family") +
  ggtitle("16S copy number from complete NCBI genomes")+
  scale_color_discrete(name = "16S") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
plot_family

pdf(paste(results.dir,"Figure_S4.pdf", sep=""), width=30, height=8)
print(plot_family)
dev.off()


