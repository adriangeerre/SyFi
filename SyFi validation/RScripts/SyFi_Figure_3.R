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

###Figure 3 - 16S distribution - SyFi output =====
table <- read.table(paste(working_directory,"16S_table.txt", sep = ""), header=T, sep="\t")
taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep = ""), header=T,sep="\t",quote="\"", fill = FALSE)

table$Family <- taxonomy$family[match(table$isolate, taxonomy$isolate)]
table_2 <- table[table$Family != "",]

#Figure 3A - 16S copy number distribution of the 447 bacterial isolates
plot_16S_CN <- ggplot(table_2, aes(x=Copy_number)) +
  geom_bar() +
  ggtitle("16S copy number") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S copy number") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_CN

#Figure 3B - 16S haplotype number distribution of the 447 bacterial isolates
plot_16S_hap <- ggplot(table_2, aes(x=Number_haplotypes)) +
  geom_bar() +
  ggtitle("16S haplotypes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S haplotypes") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_hap

#Plot A and B together
figure <- ggarrange(plot_16S_CN, plot_16S_hap,
                    labels = c("A", "B" ),font.label = list(size = 24),ncol = 2, nrow = 1, widths = 15, heights = 8)
figure

#Figure 3C - 16S copy number and haplotype distribution per bacterial family according to GTDB taxonomy

#Excluding isolates without a family taxonomic assignment
table_3 <- table_2[table_2$Family != "Unassigned",]

table_per_family <- as.data.frame(unique(table_3$Family))
table_per_family$average_CN <- "NA"
family_vector <- as.vector(table_per_family$`unique(table_3$Family)`)
table_per_family$size_family
table_per_family$size_family <- "NA"

#Calculating the number of isolates per family and the average 16S copy number
for (family in family_vector) {
  table_sub <- table_3[table_3$Family == paste(family),]
  average <- unique(ave(table_sub$Copy_number))
  table_per_family$average_CN[table_per_family$`unique(table_3$Family)` == paste(family)] <- average
  length <- length(table_sub$Family)
  table_per_family$size_family[table_per_family$`unique(table_3$Family)` == paste(family)] <- length
}

table_per_family$length_family <- paste(table_per_family$`unique(table_3$Family)`," (n = ", table_per_family$size_family,")",sep="")

#Order families from highest 16S copy number to lowest
table_per_family_2 <- table_per_family[order(as.numeric(table_per_family$average_CN),decreasing=T),]
ordered_vector <- as.vector(table_per_family_2$length_family)

table_3$length_family <- table_per_family$length_family[match(table_3$Family,table_per_family$`unique(table_3$Family)`)]
table_3$length_family <- factor(table_3$length_family, levels = ordered_vector)
table_4 <- table_3[c("isolate", "Copy_number", "Number_haplotypes", "Family", "length_family")]
colnames(table_4)[colnames(table_4) == "Copy_number"] <- "Copy number"
colnames(table_4)[colnames(table_4) == "Number_haplotypes"] <- "Haplotypes"

table_5 <- melt(table_4)

#Figure 3C - 16S copy number and haplotyp distribution across bacterial families
plot_family <- ggplot(table_5, aes(x=length_family, y=value, color=variable)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, size = 14, hjust = 1)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  ylab("16S copy number") + 
  xlab("Bacterial family") +
  ggtitle("16S copy number per family")+
  scale_color_discrete(name = "16S") +
  theme(plot.title = element_text(hjust = 0.5))
plot_family

#Fusing plot 3A, 3B, and 3C together
figure_family <- ggarrange(plot_family,
                           labels = c("C" ),font.label = list(size = 24),ncol = 1, nrow = 1, widths = 15, heights = 8)
figure_combined <- ggarrange(figure,figure_family, ncol =1, nrow =2, widths = 15, heights = c(1,2))
figure_combined

pdf(paste(results.dir,"Figure_3.pdf", sep=""), width=13, height=13)
print(figure_combined)
dev.off()
