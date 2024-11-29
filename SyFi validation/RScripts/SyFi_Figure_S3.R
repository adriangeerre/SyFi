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

###Figure S3 - 16S copy number and haplotypes comparison - SyFi vs complete genomes =====

#This script generates the figures showing the 16S copy number and haplotype distribution from closed bacterial genomes
table_16S <- read.table(paste(working_directory, "complete_genomes/16S_copy_number.tsv", sep = ""), header =T, sep = "\t")
table_16S_hap <- read.table(paste(working_directory,"complete_genomes/16S_haplotypes.tsv", sep =""), header =T,  sep = "\t")
table_16S$Haplotype <- table_16S_hap$X16S_haplotypes[match(table_16S$isolate, table_16S_hap$isolate)]

taxonomy <- read.table(paste(working_directory, "complete_genomes/taxonomy.tsv", sep =""), header =T,  sep = "\t")

table_16S_2 <- table_16S[table_16S$X16S_copy_number != 0,]

#Figure S3A
plot_16S_CN <- ggplot(table_16S_2, aes(x=X16S_copy_number)) +
  geom_bar() +
  ggtitle("16S copy number - complete genomes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S copy number") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_CN

#Figure S3B
plot_16S_hap <- ggplot(table_16S_2, aes(x=Haplotype)) +
  geom_bar() +
  ggtitle("16S haplotypes - complete genomes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S haplotypes") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_hap

#Family plot - Figure S3C
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

#Overlay with our own SynCom and SyFi results for comparison
taxonomy_AtSC = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

table_7 <- table_6[table_6$Family %in% taxonomy_AtSC$family,]

plot_family_2 <- ggplot(table_7, aes(x=length_family, y=value, color=variable)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90,size = 14, hjust = 1)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  ylab("16S copy number") + 
  xlab("Bacterial family") +
  ggtitle("16S copy number from complete NCBI genomes")+
  scale_color_discrete(name = "16S") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")
plot_family_2

figure <- ggarrange(plot_16S_CN, plot_16S_hap,
                    labels = c("A", "B" ),font.label = list(size = 24),ncol = 2, nrow = 1, widths = 15, heights = 8)

figure_family <- ggarrange(plot_family_2,
                           labels = c("C" ),font.label = list(size = 24),ncol = 1, nrow = 1, widths = 15, heights = 8)
figure_combined <- ggarrange(figure,figure_family, ncol =1, nrow =2, widths = 15, heights = c(1,2))
figure_combined

#This produces figures S3A, B, and C. 
pdf(paste(results.dir,"Figures_S3ABC.pdf", sep=""), width=13, height=13)
print(figure_combined)
dev.off()

#Figures 3D and E are generated using tidyheatmaps with the commands below
table_7_2 <- table_5[table_5$Family %in% taxonomy_AtSC$family,]

table_8 <- table_7_2[c("isolate", "Family", "variable", "value")]

new_data <- data.frame()

for (family in unique(table_8$Family)){
  table_8_sub <- table_8[table_8$Family == paste(family),]
  table_8_sub_av_16S <- sum(table_8_sub$value[table_8_sub$variable == "X16S_copy_number"])/length(table_8_sub$value[table_8_sub$variable == "X16S_copy_number"])
  table_8_sub_av_hap <- sum(table_8_sub$value[table_8_sub$variable != "X16S_copy_number"])/length(table_8_sub$value[table_8_sub$variable != "X16S_copy_number"])
  
  new_data <- rbind(new_data, data.frame(t(data.frame(c(paste(family), table_8_sub_av_16S, table_8_sub_av_hap)))))
  
}
row.names(new_data) <- NULL
colnames(new_data) <- c("Family", "16S_copy_number", "haplotypes")

#Adding the 16S copy number and haplotype distribution of the AtSC (447 isolates)
table_9 <- read.table(paste(working_directory, "16S_table.txt", sep =""), header=T, sep="\t")
taxonomy = read.table(paste(working_directory, "AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

table_9$Family <- taxonomy$family[match(table_9$isolate, taxonomy$isolate)]
table_10 <- table_9[table_9$Family != "",]
table_11 <- melt(table_10)

new_data_2 <- data.frame()

for (family in unique(table_11$Family)){
  table_8_sub <- table_11[table_11$Family == paste(family),]
  
  table_8_sub_av_16S <- sum(table_8_sub$value[table_8_sub$variable == "Copy_number"])/length(table_8_sub$value[table_8_sub$variable == "Copy_number"])
  table_8_sub_av_hap <- sum(table_8_sub$value[table_8_sub$variable != "Copy_number"])/length(table_8_sub$value[table_8_sub$variable != "Copy_number"])
  new_data_2 <- rbind(new_data_2, data.frame(t(data.frame(c(paste(family), table_8_sub_av_16S, table_8_sub_av_hap)))))
  
}
row.names(new_data_2) <- NULL
colnames(new_data_2) <- c("Family", "16S_copy_number", "haplotypes")

new_data$Database <- "NCBI"
new_data_2$Database <- "SyFi"  

table_12 <- rbind(new_data, new_data_2)
table_12_2 <- table_12[c(1,2,4)]
table_12_3 <- table_12[c(1,3,4)]

colnames(table_12_2) <- c("Family", "value", "Database")
colnames(table_12_3) <- c("Family", "value", "Database")

table_12_2$variable <- "16S_copy_number"
table_12_3$variable <- "haplotypes"
table_13 <- rbind(table_12_2,table_12_3)

ordered_vector_2 <- table_per_family_2$Family
ordered_vector_3 <- ordered_vector_2[ordered_vector_2 %in% table_13$Family]

table_13$Family <- factor(table_13$Family, levels = ordered_vector_3)
table_14 <- na.omit(table_13)
#Group samples based on the category (in this case Compartment and Plant_species)
table_15 <- dplyr::arrange(table_14,Database, .by_group = TRUE)
table_16 <- dplyr::arrange(table_15, variable, .by_group = TRUE)
table_16$variable_2 <- paste(table_16$variable, table_16$Database, sep = "_")
table_16$value <- as.numeric(table_16$value)
#Sphingomonadaceae is missing in the AtSC -447 isolates so it is removed for the heatmap (no comparison)
table_16 <- table_16[table_16$Family != "Sphingomonadaceae",]

table_16_1 <- table_16[table_16$variable == "16S_copy_number",]
table_16_2 <- table_16[table_16$variable != "16S_copy_number",]

table_16_1$variable <- gsub("16S_copy_number", "16S copy number", table_16_1$variable)
table_16_2$variable <- gsub("haplotypes", "16S haplotypes", table_16_2$variable)

#Figure S3D
Heatmap <- tidy_heatmap(table_16_1,rows = "Family", columns = "Database", values = "value", annotation_col = c(Database, variable),scale = "none", gaps_col = variable, fontsize_row = 6, main = "SyFi vs NCBI")

#Figure S3E
Heatmap_2 <- tidy_heatmap(table_16_2,rows = "Family", columns = "Database", values = "value", annotation_col = c(Database, variable),scale = "none", gaps_col = variable, fontsize_row = 6, main = "SyFi vs NCBI")

ggsave("Figure_S3D.png", plot = Heatmap, width = 12, height = 10, units = "cm", device = "png", dpi = 300, path = results.dir)
ggsave("Figure_S3E.png", plot = Heatmap_2, width = 12, height = 10, units = "cm", device = "png", dpi = 300, path = results.dir)



