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

###Figure S6 - Percentage mapped reads=====
reads <- read.table(paste(working_directory,"reads_percentage.txt", sep =""), header = T, sep ="\t")
meta <- read.table(paste(working_directory,"metadata.txt", sep = ""), header = T, sep ="\t")

reads_2 <- melt(reads)
reads_2$Compartment <- meta$Compartment[match(reads_2$variable, meta$sample_id)]

reads_V3V4 <- reads_2[grep("V3V4", reads_2$Sample), ]
reads_V5V7 <- reads_2[grep("V5V7", reads_2$Sample), ]

reads_V3V4$Percentage <- reads_V3V4$value * 100
reads_V5V7$Percentage <- reads_V5V7$value * 100

#Order
reads_V3V4$Sample <- gsub("V3V4", "V3-V4", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("V5V7", "V5-V7", reads_V5V7$Sample)

reads_V3V4$Sample <- gsub("_", "% ", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("_", "% ", reads_V5V7$Sample)

reads_V3V4$Sample <- gsub("SyFi% ", "SyFi ", reads_V3V4$Sample)
reads_V5V7$Sample <- gsub("SyFi% ", "SyFi ", reads_V5V7$Sample)

reads_V3V4$Sample <- factor(reads_V3V4$Sample, levels = c("95% V3-V4", "97% V3-V4", "98% V3-V4","99% V3-V4", "100% V3-V4", "SyFi V3-V4"))
reads_V5V7$Sample <- factor(reads_V5V7$Sample, levels = c("95% V5-V7", "97% V5-V7", "98% V5-V7","99% V5-V7", "100% V5-V7", "SyFi V5-V7"))

reads_V3V4$Compartment <- gsub("RZ", "Rhizosphere", reads_V3V4$Compartment)
reads_V3V4$Compartment <- gsub("ES", "Endosphere", reads_V3V4$Compartment)
reads_V5V7$Compartment <- gsub("RZ", "Rhizosphere", reads_V5V7$Compartment)
reads_V5V7$Compartment <- gsub("ES", "Endosphere", reads_V5V7$Compartment)

reads_V3V4$Compartment <- factor(reads_V3V4$Compartment, levels = c("Input", "Endosphere", "Rhizosphere"))
reads_V5V7$Compartment <- factor(reads_V5V7$Compartment, levels = c("Input", "Endosphere", "Rhizosphere"))

plot_V3V4 <- ggplot(reads_V3V4, aes(x=Sample, y=Percentage)) + 
  geom_violin(scale="width") +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(colour = reads_V3V4$Compartment),show.legend = T) +
  ylab("Mapped/pseudoaligned reads (%)") + 
  xlab("Metric") +
  ggtitle("Mapped/pseudoaligned reads to SynCom isolates") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(strip.text.x = element_text(size = 14))+
  guides(colour=guide_legend(title="Sample Type")) + 
  facet_wrap(~Compartment,scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle = 25, hjust =1), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(strip.text.x = element_text(size = 14))
plot_V3V4

plot_V3V4_2 <- plot_V3V4 + stat_compare_means(comparisons = list(c("SyFi V3-V4", "100% V3-V4")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))

plot_V5V7 <- ggplot(reads_V5V7, aes(x=Sample, y=Percentage)) + 
  geom_violin(scale="width") +
  geom_jitter(shape=16, position=position_jitter(0.2), aes(colour = reads_V3V4$Compartment),show.legend = T) +
  ylab("Mapped/pseudoaligned reads (%)") + 
  xlab("Metric") +
  #ggtitle("Mapped reads to SynCom isolates") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  guides(colour=guide_legend(title="Sample Type")) +
  facet_wrap(~Compartment,scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle = 25, hjust =1), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(strip.text.x = element_text(size = 14))
plot_V5V7

plot_V5V7_2 <- plot_V5V7 + stat_compare_means(comparisons = list(c("SyFi V5-V7", "100% V5-V7")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))

figure_1 <- ggarrange(plot_V3V4_2, plot_V5V7_2,
                      labels = c("A", "B" ),font.label = list(size = 28),ncol = 1, nrow = 2, common.legend = T)

pdf(paste(results.dir,"Figure_S6.pdf", sep=""), width=18, height=12)
print(figure_1)
dev.off()

