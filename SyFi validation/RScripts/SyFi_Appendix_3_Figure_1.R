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

###Appendix 3 - Figure 1 - 16S copy number supplemental figure =====
#comparing contamination levels with 16S occurences
table <- read.table(paste(working_directory,"16S_table.txt", sep =""), header=T, sep="\t")
taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

table$Family <- taxonomy$family[match(table$isolate, taxonomy$isolate)]
table_2 <- table[table$Family != "",]

completeness <- read.table(paste(working_directory,"genome_assembly/checkm.tsv",sep =""), header = T, sep = "\t")

table_2$Contamination <- completeness$Contamination[match(table_2$isolate, completeness$Isolate)]
table_2$Heterogeneity <- completeness$Strain_heterogeneity[match(table_2$isolate, completeness$Isolate)]

table_3 <- table_2[, c(1,3,5)]
table_4 <- table_2[, c(1,3,6)]

colnames(table_3) <- c("Isolate", "Copy_number", "value")
colnames(table_4) <- c("Isolate", "Copy_number", "value")

table_3$Metric <- "Contamination"
table_4$Metric <- "Heterogeneity"

table_5 <- rbind(table_3, table_4)

plot_occ <- ggscatter(table_5, x = "Copy_number", y = "value", color = "Metric") + 
  ggtitle("Contamination and heterogeneity") + 
  stat_cor(aes(color = Metric)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("(%)") + 
  xlab("16S copy number") +
  theme(legend.position = "right") +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18))
plot_occ

#GC content
table <- read.table(paste(working_directory,"16S_table.txt", sep = ""), header=T, sep="\t")
taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""),header=T,sep="\t",quote="\"", fill = FALSE)

table$Family <- taxonomy$family[match(table$isolate, taxonomy$isolate)]
table_2 <- table[table$Family != "",]

#This was the 16S fingerprints
fin_length <- read.table(paste(working_directory,"genome_assembly/length_16S_fin.txt", sep =""), header = F, sep = "\t")
GC_length <- read.table(paste(working_directory,"genome_assembly/GC_16S_fin.txt", sep =""), header = F, sep = "\t")

table_2$length <- fin_length$V2[match(table_2$isolate, fin_length$V1)]
table_2$GC_length <- GC_length$V2[match(table_2$isolate, GC_length$V1)]
table_2$GC_content <- table_2$GC_length/table_2$length * 100

#In comparison with Genomic GC content
GC_content <- read.table(paste(working_directory,"genome_assembly/GC_content.txt", sep = ""),header =F, sep = "\t")
Genome_content <- read.table(paste(working_directory,"genome_assembly/Genome_size.txt", sep =""), header =F, sep = "\t")

table_2$GC_content_genome <- GC_content$V2[match(table_2$isolate, GC_content$V1)]
table_2$Genome_content <- Genome_content$V2[match(table_2$isolate, Genome_content$V1)]
table_2$GC_content_perc <- table_2$GC_content_genome/table_2$Genome_content *100

table_2$Diff <- table_2$GC_content/table_2$GC_content_perc
#Getting coefficient
model <- lm(Copy_number ~ Diff, data = table_2)

plot_GC <- ggscatter(table_2, x = "Copy_number", y = "Diff") + 
  ggtitle("GC content - SyFi ") + 
  stat_cor(label.x = 1.2) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("GC content 16S/Genome") + 
  xlab("16S copy number") +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18))
plot_GC


###Adding complete genomes
genomes_GC <- read.table(paste(working_directory,"complete_genomes/16S_GC_genomes_full.txt", sep =""), header = T, sep = "\t")
genomes_GC$Diff <- (genomes_GC$X16S_GC/genomes_GC$X16S)/(genomes_GC$genome_GC/genomes_GC$genome)

table_16S <- read.table(paste(working_directory,"complete_genomes/16S_copy_number.tsv", sep =""),header =T, sep = "\t")

genomes_GC$Copy_number <- table_16S$X16S_copy_number[match(genomes_GC$Isolate, table_16S$isolate)]
genomes_GC_2 <- genomes_GC[genomes_GC$Copy_number != 0,]

model_2 <- lm(Copy_number ~ Diff, data = genomes_GC_2)

plot_GC_comp <- ggscatter(genomes_GC_2, x = "Copy_number", y = "Diff") + 
  ggtitle("GC content - NCBI") + 
  stat_cor(label.x = 1.2) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("GC content 16S/Genome") + 
  xlab("16S copy number") +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18))
plot_GC_comp

#Overestimation
over <- data.frame(seq(0.5,1.5,0.1))
colnames(over) <- "GC_content"
over$CN_genome <- over$GC_content*1.5 + 2.9
over$CN_own <- over$GC_content*15.1 -8.7
over$difference <- over$CN_own/over$CN_genome

model_2 <- lm(difference ~ GC_content, data = over)
text <- c("16S_CN SyFi = 15.1* diff GC - 8.7","16S_CN NCBI = 1.5* diff GC + 2.9", "Diff SyFi-NCBI = 16S_CN SyFi/16S_CN NCBI")
tl <- textGrob(paste(strwrap(text, 40), collapse="\n"), hjust=0, x=0,gp = gpar(col = "black", fontsize = 10))

over_2 <- over[,c(1,2,4)]
over_3 <- over[,c(1,3,4)]
over_2$Group <- "NCBI"
over_3$Group <- "SyFi"

colnames(over_2) <- c("GC_content", "CN", "difference", "Group")
colnames(over_3) <- c("GC_content", "CN", "difference", "Group")

over_4 <- rbind(over_2, over_3)

plot_line <- ggplot() + 
  geom_line(aes(x = over$GC_content, y = over$difference),color = "black",size=1.2) +
  geom_line(aes(x = over_4$GC_content, y = over_4$CN, color = over_4$Group),size=1.2) +
  ggtitle("Overestimation of 16S copy number") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("overestimation (difference slope SyFi-NCBI)") + 
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") + 
  geom_vline(xintercept = 0.85, color = "grey40", linetype = "dashed") + 
  xlab("diff GC content 16S vs genome") +
  labs(color="Genomes") +
  scale_y_continuous(name = "Overestimation",sec.axis = dup_axis(name="16S copy number")) +
  theme(axis.text.x = element_text(size = 12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=18)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_line

plot_line_2 <- plot_line + inset_element(tl, left = 0.05, bottom = 0.78, right =0.75, top = 1) 

figure <- ggarrange(plot_occ, plot_GC,plot_GC_comp, plot_line_2, labels = c("A", "B", "C", "D"),font.label = list(size = 18), ncol =2, nrow =2, widths = c(6,5))

pdf(paste(results.dir,"Appendix_Figure_1.pdf", sep=""), width=14, height=14)
print(figure)
dev.off()


