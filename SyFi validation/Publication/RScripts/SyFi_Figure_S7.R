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

###Figure S7 - Percentage matched reads =====
reads <- read.table(paste(working_directory,"reads_percentage.txt", sep =""), header = T, sep ="\t")
meta <- read.table(paste(working_directory,"metadata.txt", sep = ""), header = T, sep ="\t")

reads_2 <- melt(reads)
reads_2$Compartment <- meta$Compartment[match(reads_2$variable, meta$sample_id)]

reads_V3V4 <- reads_2[grep("V3V4", reads_2$Sample), ]
reads_V5V7 <- reads_2[grep("V5V7", reads_2$Sample), ]

reads_V3V4$Percentage <- reads_V3V4$value * 100
reads_V5V7$Percentage <- reads_V5V7$value * 100

#Order
reads_V3V4$Sample <- gsub("V3V4", "", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("V5V7", "", reads_V5V7$Sample)

reads_V3V4$Sample <- gsub("_", "%", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("_", "%", reads_V5V7$Sample)

reads_V3V4$Sample <- gsub("%%", ":", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("%%", ":", reads_V5V7$Sample)

reads_V3V4$Sample <- gsub("SyFi-", "SyFi", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("SyFi-", "SyFi", reads_V5V7$Sample)

reads_V3V4$Sample <- gsub("Mothur%", "Mothur", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("Mothur%", "Mothur", reads_V5V7$Sample)

reads_V3V4$Sample <- gsub("Salmon-", "Salmon", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("Salmon-", "Salmon", reads_V5V7$Sample)

reads_V3V4$Group <- "Qiime2-VSEARCH"
reads_V5V7$Group <- "Qiime2-VSEARCH"
reads_V3V4$Group[grep("Amplicon",reads_V3V4$Sample)] <- "Kraken2 - Amplicon"
reads_V5V7$Group[grep("Amplicon",reads_V5V7$Sample)] <- "Kraken2 - Amplicon"
reads_V3V4$Group[grep("Fingerprint",reads_V3V4$Sample)] <- "Kraken2 - Fingerprint"
reads_V5V7$Group[grep("Fingerprint",reads_V5V7$Sample)] <- "Kraken2 - Fingerprint"
reads_V3V4$Group[grep("SyFi",reads_V3V4$Sample)] <- "SyFi"
reads_V5V7$Group[grep("SyFi",reads_V5V7$Sample)] <- "SyFi"
reads_V3V4$Group[grep("Salmon",reads_V3V4$Sample)] <- "Salmon"
reads_V5V7$Group[grep("Salmon",reads_V5V7$Sample)] <- "Salmon"
reads_V3V4$Group[grep("Mothur",reads_V3V4$Sample)] <- "Mothur"
reads_V5V7$Group[grep("Mothur",reads_V5V7$Sample)] <- "Mothur"

reads_V3V4$Sample <- gsub("9", "VS-9", reads_V3V4$Sample)
reads_V3V4$Sample <- gsub("9VS-", "9", reads_V3V4$Sample)
reads_V3V4$Sample <- gsub("100", "VS-100", reads_V3V4$Sample)
reads_V3V4$Sample <- gsub(":Amplicon", "-A", reads_V3V4$Sample)
reads_V3V4$Sample <- gsub(":Fingerprint", "-F", reads_V3V4$Sample)
reads_V3V4$Sample <- gsub("Kraken", "KM-", reads_V3V4$Sample)

reads_V5V7$Sample <- gsub("9", "VS-9", reads_V5V7$Sample)
reads_V5V7$Sample <- gsub("9VS-", "9", reads_V5V7$Sample)
reads_V5V7$Sample <- gsub("100", "VS-100", reads_V5V7$Sample)
reads_V5V7$Sample <- gsub(":Amplicon", "-A", reads_V5V7$Sample)
reads_V5V7$Sample <- gsub(":Fingerprint", "-F", reads_V5V7$Sample)
reads_V5V7$Sample <- gsub("Kraken", "KM-", reads_V5V7$Sample)

reads_V3V4$Sample <- factor(reads_V3V4$Sample, levels = c("VS-95%", "VS-97%", "VS-98%","VS-99%", "VS-100%", "KM-2-A", "KM-3-A","KM-5-A","KM-2-F", "KM-3-F","KM-5-F","Mothur","Salmon","SyFi"))
reads_V5V7$Sample <- factor(reads_V5V7$Sample, levels = c("VS-95%", "VS-97%", "VS-98%","VS-99%", "VS-100%", "KM-2-A", "KM-3-A","KM-5-A","KM-2-F", "KM-3-F","KM-5-F","Mothur","Salmon","SyFi"))

reads_V3V4$Compartment <- gsub("RZ", "Rhizosphere", reads_V3V4$Compartment)
reads_V3V4$Compartment <- gsub("ES", "Endosphere", reads_V3V4$Compartment)
reads_V5V7$Compartment <- gsub("RZ", "Rhizosphere", reads_V5V7$Compartment)
reads_V5V7$Compartment <- gsub("ES", "Endosphere", reads_V5V7$Compartment)

reads_V3V4$Compartment <- factor(reads_V3V4$Compartment, levels = c("Input", "Endosphere", "Rhizosphere"))
reads_V5V7$Compartment <- factor(reads_V5V7$Compartment, levels = c("Input", "Endosphere", "Rhizosphere"))

reads_V3V4$Sample_2 <- gsub("-","_",reads_V3V4$Sample)
reads_V5V7$Sample_2 <- gsub("-","_",reads_V5V7$Sample)

#Statistics
#Stats
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(compartment in unique(reads_V3V4$Compartment)) {
  # Subset data for the current plant
  reads_V3V4_sub <- reads_V3V4[reads_V3V4$Compartment == paste(compartment), ]
  # Perform ANOVA
  fitAnova <- aov(Percentage ~ Sample_2, data=reads_V3V4_sub)

  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Sample_2$Letters
  
  # Store results
  anova_results[[compartment]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(compartment) {
  data.frame(Compartment = compartment, Sample_2 = names(anova_results[[compartment]]), Letters = anova_results[[compartment]])
}))
# Order factors based on your original setup
ltlbl_combined$Compartment <- factor(ltlbl_combined$Compartment, levels = c("Input", "Endosphere", "Rhizosphere"))
ltlbl_combined$Sample <- reads_V3V4$Sample[match(ltlbl_combined$Sample_2,reads_V3V4$Sample_2 )]
ltlbl_combined$Sample <- factor(ltlbl_combined$Sample, levels =  c("VS-95%", "VS-97%", "VS-98%","VS-99%", "VS-100%", "KM-2-A", "KM-3-A","KM-5-A","KM-2-F", "KM-3-F","KM-5-F","Mothur","Salmon","SyFi"))
ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Compartment, ltlbl_combined$Sample), ]

plot_V3V4 <- ggplot(reads_V3V4, aes(x=Sample, y=Percentage, colour =Group)) + 
  geom_violin(scale="width") +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(colour = Group),show.legend = T) +
  ylab("Mapped/pseudoaligned reads (%)") + 
  xlab("Metric") +
  ggtitle("Mapped/pseudoaligned reads to SynCom isolates") + 
  scale_color_manual(values =c("#d55e00","#cc79a7","#0072b2","#f0e442","#009e73","#8BB5E7") ) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_text(size = 14))+
  guides(colour=guide_legend(title="Method")) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle = 25, hjust =1), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(strip.text.x = element_text(size = 14)) +
  facet_wrap(~Compartment,scales = "free_y") +
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(Percentage)*1.05), show.legend=FALSE)
plot_V3V4

#Statistics V5V7

#Statistics
#Stats
anova_results_2 <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(compartment in unique(reads_V5V7$Compartment)) {
  # Subset data for the current plant
  reads_V3V4_sub <- reads_V3V4[reads_V3V4$Compartment == paste(compartment), ]
  # Perform ANOVA
  fitAnova <- aov(Percentage ~ Sample_2, data=reads_V3V4_sub)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Sample_2$Letters
  
  # Store results
  anova_results_2[[compartment]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined_2 <- do.call(rbind, lapply(names(anova_results_2), function(compartment) {
  data.frame(Compartment = compartment, Sample_2 = names(anova_results_2[[compartment]]), Letters = anova_results_2[[compartment]])
}))
# Order factors based on your original setup
ltlbl_combined_2$Compartment <- factor(ltlbl_combined_2$Compartment, levels = c("Input", "Endosphere", "Rhizosphere"))
ltlbl_combined_2$Sample <- reads_V3V4$Sample[match(ltlbl_combined_2$Sample_2,reads_V5V7$Sample_2 )]
ltlbl_combined_2$Sample <- factor(ltlbl_combined_2$Sample, levels =  c("VS-95%", "VS-97%", "VS-98%","VS-99%", "VS-100%", "KM-2-A", "KM-3-A","KM-5-A","KM-2-F", "KM-3-F","KM-5-F","Mothur","Salmon","SyFi"))
ltlbl_combined_2 <- ltlbl_combined_2[order(ltlbl_combined_2$Compartment, ltlbl_combined_2$Sample), ]

plot_V5V7 <- ggplot(reads_V5V7, aes(x=Sample, y=Percentage, colour =Group)) + 
  geom_violin(scale="width") +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(colour = Group),show.legend = T) +
  ylab("Mapped/pseudoaligned reads (%)") + 
  xlab("Metric") +
  #ggtitle("Mapped/pseudoaligned reads to SynCom isolates") + 
  scale_color_manual(values =c("#d55e00","#cc79a7","#0072b2","#f0e442","#009e73","#8BB5E7") ) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_text(size = 14))+
  guides(colour=guide_legend(title="Method")) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle = 25, hjust =1), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(strip.text.x = element_text(size = 14)) +
  facet_wrap(~Compartment,scales = "free_y") +
  stat_summary(geom = 'text', label = ltlbl_combined_2$Letters, fun.y = max, aes(y = max(Percentage)*1.05), show.legend=FALSE)
plot_V5V7

figure_1 <- ggarrange(plot_V3V4, plot_V5V7,
                      labels = c("A", "B" ),font.label = list(size = 28),ncol = 1, nrow = 2, common.legend = T)

pdf(paste(results.dir,"Figure_S7.pdf", sep=""), width=24, height=12)
print(figure_1)
dev.off()

