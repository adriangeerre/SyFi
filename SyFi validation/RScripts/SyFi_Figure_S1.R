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

###Figure S1 - Genome assembly effect on SyFi output =====

#Figure S1A - Genome completeness plot
#Loading the files
completeness <- read.table(paste(working_directory,"genome_assembly/checkm.tsv", sep =""), header = T, sep = "\t")
strains_with_fingerprint <- as.vector(read.table(paste(working_directory, "strains_with_fingerprint.txt", sep =""), header =F, sep = "\t"))

tax_df = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)
rownames(tax_df) <- tax_df$isolate
tax_df_2 <- tax_df %>% dplyr::select (-isolate)
colnames(tax_df_2)=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "SynCom")

tax_df_AtSC <- tax_df_2[tax_df_2$SynCom == "AtSC",]

tax_df_AtSC$Completeness <- completeness$Completeness[match(row.names(tax_df_AtSC), completeness$Isolate)]

tax_df_AtSC$Fingerprint <- "Absent"
tax_df_AtSC$Fingerprint[row.names(tax_df_AtSC) %in% strains_with_fingerprint$V1] <- "Present"

plot_comp <- ggplot(tax_df_AtSC, aes(x=Fingerprint , y=log10(Completeness), color = Fingerprint)) + 
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Genome completeness") + 
  geom_jitter(aes(size = 3), cex =3)+
  theme(legend.position="top")+ 
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size =12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18)) +  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("Genome completeness (%)") + 
  xlab("Fingerprint")
plot_comp_2 <- plot_comp + stat_compare_means(comparisons = list(c("Present", "Absent")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
plot_comp_2

#Figure S1B and C - Genome coverage and GC content
coverage <- read.table(paste(working_directory,"genome_assembly/Genome_total_base_pairs.txt",sep =""), header =F, sep = "\t")
tax_df_AtSC$coverage <- coverage$V2[match(row.names(tax_df_AtSC), coverage$V1)]
GC_content <- read.table(paste(working_directory,"genome_assembly/GC_content.txt", sep =""), header =F, sep = "\t")
Genome_content <- read.table(paste(working_directory,"genome_assembly/Genome_size.txt", sep=""), header =F, sep = "\t")

tax_df_AtSC$GC_content <- GC_content$V2[match(row.names(tax_df_AtSC), GC_content$V1)]
tax_df_AtSC$Genome_size <- Genome_content$V2[match(row.names(tax_df_AtSC), Genome_content$V1)]
tax_df_AtSC$GC_perc <- tax_df_AtSC$GC_content/tax_df_AtSC$Genome_size * 100

tax_df_AtSC$coverage_2 <- log10(tax_df_AtSC$coverage/tax_df_AtSC$Genome_size)

tax_df_AtSC_2 <- na.omit(tax_df_AtSC)

plot_comp_3 <- ggplot(tax_df_AtSC_2, aes(x=Fingerprint , y=coverage_2, color = Fingerprint)) + 
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Genome coverage") + 
  geom_jitter(aes(size = 3), cex =3)+
  theme(legend.position="top")+ 
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size =12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18)) +  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("log10(Coverage)") + 
  xlab("Fingerprint")
plot_comp_4 <- plot_comp_3 + stat_compare_means(comparisons = list(c("Present", "Absent")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
plot_comp_4


plot_GC <- ggplot(tax_df_AtSC, aes(x=Fingerprint , y=GC_perc, color = Fingerprint)) + 
  geom_boxplot(outlier.shape = NA)+
  ggtitle("Genomic GC content") + 
  geom_jitter(aes(size = 3), cex =3)+
  theme(legend.position="top")+ 
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size =12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18)) +  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("GC content (%)") + 
  xlab("Fingerprint")
plot_GC_2 <- plot_GC + stat_compare_means(comparisons = list(c("Present", "Absent")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
plot_GC_2

#Figure S1D - L50
L50 <- read.table(paste(working_directory,"genome_assembly/L50.txt", sep =""), header =F, sep = "\t")

tax_df_AtSC$L50 <- log10(L50$V2[match(row.names(tax_df_AtSC), L50$V1)])
tax_df_AtSC$L50[tax_df_AtSC$L50 == 0] <- NA
tax_df_AtSC$L50[tax_df_AtSC$L50 == "-Inf"] <- NA

tax_df_AtSC_3 <- na.omit(tax_df_AtSC)

plot_L50 <- ggplot(tax_df_AtSC_3, aes(x=Fingerprint , y=L50, color = Fingerprint)) + 
  geom_boxplot(outlier.shape = NA)+
  ggtitle("L50") + 
  geom_jitter(aes(size = 3), cex =3)+
  theme(legend.position="top")+ 
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size =12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18)) +  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("log10(L50)") + 
  xlab("Fingerprint")
plot_L50_2 <- plot_L50 + stat_compare_means(comparisons = list(c("Present", "Absent")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
plot_L50_2

#Put all four panels together and print the output
all <- ggarrange(plot_comp_2, plot_comp_4, plot_GC_2, plot_L50_2,labels = c("A", "B", "C", "D"),font.label = list(size = 18), ncol =2, nrow =2)

pdf(paste(results.dir,"Figure_S1.pdf", sep=""), width=10, height=8)
print(all)
dev.off()


