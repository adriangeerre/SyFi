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

###Figure S1 - SyFi computing time =====
Time_table <- read.table(paste(working_directory,"Time.txt", sep = ""), header =T, sep = "\t")

Time_table$File.size..MB. <- as.numeric(Time_table$File.size..MB.)
Time_table$Duration <- as.numeric(Time_table$Duration)
Time_table$Threads <- factor(Time_table$Threads, levels =c("1","2","4","8","16","24"))

plot_occ <- ggscatter(Time_table, x = "Duration", y = "File.size..MB.", color = "Threads") + 
  ggtitle("SyFi main computing time") + 
  scale_color_manual(values =c("#d55e00","#cc79a7","#0072b2","#f0e442","#009e73","#8BB5E7") ) +
  ylab("Genomic read file size (Mbp)") + 
  xlab("Duration (seconds)") +
  labs(color = "Threads") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10), axis.title = element_text(size = 12), 
        axis.text.y = element_text(size=10), legend.title = element_text(size=12), 
        legend.text = element_text(size=10)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

plot_occ

pdf(paste(results.dir,"Figure_S1.pdf", sep=""), width=10, height=6)
print(plot_occ)
dev.off()
