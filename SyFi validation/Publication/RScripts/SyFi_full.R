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

###Figure 2 - 16S distribution - human gut dataset - SyFi output=====
table <- read.table(paste(working_directory,"human_set/SyFi_human_16S.txt",sep = ""), header=T, sep="\t")
taxonomy = read.table(paste(working_directory, "human_set/SyFi_human_taxonomy.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)
derep <- read.table(paste(working_directory, "human_set/derep_list_SyFi_human.txt", sep = ""))

table$Family <- taxonomy$Family[match(table$isolate, taxonomy$Identifier)]
table_2 <- table[table$Family != "",]
table_3 <- table_2[table_2$isolate %in% derep$V1,]

#Figure 3A
plot_16S_CN <- ggplot(table_3, aes(x=Copy_number)) +
  geom_bar() +
  ggtitle("16S copy number") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S copy number") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_CN

#Figure 2B
plot_16S_hap <- ggplot(table_3, aes(x=No_of_haplotypes)) +
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

#Figure 2C - 16S copy number and haplotype distribution per bacterial family according to GTDB taxonomy

#Excluding isolates without a family taxonomic assignment
table_4 <- table_3[table_3$Family != "Unassigned",]

table_per_family <- as.data.frame(unique(table_4$Family))
table_per_family$average_CN <- "NA"
family_vector <- as.vector(table_per_family$`unique(table_4$Family)`)
table_per_family$size_family
table_per_family$size_family <- "NA"

#Calculating the number of isolates per family and the average 16S copy number
for (family in family_vector) {
  table_sub <- table_4[table_4$Family == paste(family),]
  average <- unique(ave(table_sub$Copy_number))
  table_per_family$average_CN[table_per_family$`unique(table_4$Family)` == paste(family)] <- average
  length <- length(table_sub$Family)
  table_per_family$size_family[table_per_family$`unique(table_4$Family)` == paste(family)] <- length
}

table_per_family$length_family <- paste(table_per_family$`unique(table_4$Family)`," (n = ", table_per_family$size_family,")",sep="")

#Order families from highest 16S copy number to lowest
table_per_family_2 <- table_per_family[order(as.numeric(table_per_family$average_CN),decreasing=T),]
ordered_vector <- as.vector(table_per_family_2$length_family)

table_4$length_family <- table_per_family$length_family[match(table_4$Family,table_per_family$`unique(table_4$Family)`)]
table_4$length_family <- factor(table_4$length_family, levels = ordered_vector)
table_5 <- table_4[c("isolate", "Copy_number", "No_of_haplotypes", "Family", "length_family")]
colnames(table_5)[colnames(table_5) == "Copy_number"] <- "Copy number"
colnames(table_5)[colnames(table_5) == "Number_haplotypes"] <- "Haplotypes"

table_6 <- melt(table_5)

#Figure 2C - 16S copy number and haplotyp distribution across bacterial families
plot_family <- ggplot(table_6, aes(x=length_family, y=value, color=variable)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(face = "italic",angle = 90, size = 14, hjust = 1)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  ylab("16S copy number") + 
  xlab("Bacterial family") +
  ggtitle("16S copy number per family")+
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
  guides(colour = guide_legend(title = "16S")) +
  theme(plot.title = element_text(hjust = 0.5))
plot_family

#Fusing plot 3A, 3B, and 3C together
figure_family <- ggarrange(plot_family,
                           labels = c("C" ),font.label = list(size = 24),ncol = 1, nrow = 1, widths = 15, heights = 8)
figure_combined <- ggarrange(figure,figure_family, ncol =1, nrow =2, widths = 15, heights = c(1,2))
figure_combined

pdf(paste(results.dir,"Figure_2.pdf", sep=""), width=13, height=13)
print(figure_combined)
dev.off()

###Figure 3 - 16S distribution - plant root dataset - SyFi output =====
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
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
  guides(colour = guide_legend(title = "16S")) +
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


###Figure 4a - genus - NRMSE =====
#Function relative abundance (calculates relative abundance per sample from a long format table)
Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

#Taxonomy file
taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

#Vector of each dataset - V3V4 or V5V7 amplicons/fingerprints - Qiime2-VSEARCH clustering at 95%, 97%, 98%, 99%, 100% and SyFi (Salmon pseudoalignment of amplicon reads to SyFi-generated fingerprints)
vector <- c("95_V3V4", "97_V3V4", "98_V3V4", "99_V3V4", "100_V3V4", 
            "Salmon_V3V4", "SyFi_V3V4","Mothur_V3V4",
            "Kraken2:Amplicon_V3V4","Kraken3:Amplicon_V3V4","Kraken5:Amplicon_V3V4",
            "Kraken2:Fingerprint_V3V4","Kraken3:Fingerprint_V3V4","Kraken5:Fingerprint_V3V4",
            "95_V5V7", "97_V5V7", "98_V5V7", "99_V5V7", "100_V5V7", 
            "Salmon_V5V7","SyFi_V5V7", "Mothur_V5V7",
            "Kraken2:Amplicon_V5V7","Kraken3:Amplicon_V5V7","Kraken5:Amplicon_V5V7",
            "Kraken2:Fingerprint_V5V7","Kraken3:Fingerprint_V5V7","Kraken5:Fingerprint_V5V7")
#Create data frame from vector 
data_tog_3 <- data.frame()

#For loop that loads in each dataset of the vector
#It loads the same shotgun-metagenome-sequenced dataset every time
#And the clustering output that says which SynCom isolate amplicon sequences or fingerprints cluster together based on the implemented sequence identity clustering percentage
#For SyFi this percentage is set on 100%

#the loop will cluster the microbiome abundances of the dataset according to the clustering
#Similarly it will use the same clustering for the isolate counts in the shotgun metagenome-sequenced dataset
#Isolates with 0 counts in the loaded amplicon dataset or the shotgun metagenome-sequenced dataset are discarded

for (metric in vector) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory, "comparison_datasets/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert dataset for comparison (SyFi, or Qiime2-VSEARCH-generated table)
  otu_table = read.table(paste0(working_directory, "comparison_datasets/AtSC_", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep = ""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  taxonomy <- read.table(paste(working_directory, "AtSC_taxonomy_GTDB.tsv",sep= ""), sep = "\t", row.names =1, header =T)
  
  #Reformatting table of interest 
  no_of_genera <- length(na.omit(unique(taxonomy$genus)))
  
  #Create new table with the right dimensions
  new_table <- data.frame(matrix(NA, nrow = no_of_genera, ncol = length(row.names(samples_meta))+1))
  
  genera_all <- na.omit(unique(taxonomy$genus))
  
  for (number in 1:no_of_genera){
    number_2 <- genera_all[number]
    new_table[number,1] <- paste(number_2)
  }
  
  column_names <- c("Genus", paste(row.names(samples_meta),sep=""))
  colnames(new_table) <- column_names
  
  row.names(new_table) <- new_table$Genus
  new_table <- new_table %>% dplyr::select(-Genus)
  
  #Rename and reorder column names of otu_table
  colnames(otu_table) = gsub("_t3", "", colnames(otu_table))
  colnames(otu_table) = gsub("_t0", "", colnames(otu_table))
  colnames(otu_table) = gsub("_undil", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v5_v7", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v3_v4", "", colnames(otu_table))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_table_2 <- otu_table[column_names_order]
  otu_table <- otu_table_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:length(genera_all)){
    genus_pick <- genera_all[row]
    
    isolates <- na.omit(row.names(taxonomy)[taxonomy$genus == paste(genus_pick)])
    
    if (length(isolates) > 1){
      other_list <- row.names(otu_table) 
      the_isolate <- isolates[isolates %in% other_list]
    } else {
      the_isolate <- isolates
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_table[row.names(otu_table) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    new_table[row,] <- otu_table_3
  }
  
  #Create new table with the right dimensions for metagenome table
  meta_table <- data.frame(matrix(NA, nrow = no_of_genera, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_genera){
    number_2 <- genera_all[number]
    meta_table[number,1] <- paste(number_2)
  }
  
  column_names <- c("Genus", paste(row.names(samples_meta),sep=""))
  colnames(meta_table) <- column_names
  
  row.names(meta_table) <- meta_table$Genus
  meta_table <- meta_table %>% dplyr::select(-Genus)
  
  #Reshuffle columns of the metagenome data
  colnames(otu_meta) = gsub("_undil", "", colnames(otu_meta))
  colnames(otu_meta) = gsub("_HL_", "_HL_orig_", colnames(otu_meta))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_meta_2 <- otu_meta[colnames(otu_meta) %in% column_names_order]
  otu_meta <- otu_meta_2[column_names_order]
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:length(genera_all)){
    genus_pick <- genera_all[row]
    isolates <- na.omit(row.names(taxonomy)[taxonomy$genus == paste(genus_pick)])
    
    if (length(isolates) > 1){
      other_list <- row.names(otu_meta) 
      the_isolate <- isolates[isolates %in% other_list]
    } else {
      the_isolate <- isolates
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_meta[row.names(otu_meta) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    meta_table[row,] <- otu_table_3
  }
  
  #Melt tables
  meta_table$Genus <- row.names(meta_table)
  meta_table_2 <- melt(meta_table)
  new_table$Genus <- row.names(new_table)
  new_table_2 <- melt(new_table)
  
  #Make relative abundances
  meta_table_3 <- Hank_the_normalizer(meta_table_2,"variable","value")
  new_table_3 <- Hank_the_normalizer(new_table_2,"variable","value")
  
  meta_table_3$ID <- paste(meta_table_3$variable,meta_table_3$Genus,sep="_")
  new_table_3$ID <- paste(new_table_3$variable,new_table_3$Genus,sep="_")
  
  fused_table <- data.frame(matrix(NA, nrow = length(meta_table_3$ID), ncol = 0))
  fused_table$V1 <- meta_table_3$Rel
  row.names(fused_table) <- meta_table_3$ID
  fused_table$V2 <- new_table_3$Rel[match(row.names(fused_table),new_table_3$ID)]
  colnames(fused_table) <- c("metagenome_abundances", "metric_abundances")
  
  fused_table_2 <- na.omit(fused_table[fused_table$metric_abundances != 0,])
  fused_table_3 <- na.omit(fused_table_2[fused_table_2$metagenome_abundances != 0,])
  fused_table_3$Genus <- meta_table_3$Genus[match(row.names(fused_table_3), meta_table_3$ID)]
  
  data_tog_2 <- data.frame()
  
  for (genus in unique(fused_table_3$Genus)){
    fused_table_3_sub <- fused_table_3[fused_table_3$Genus == paste(genus),]
    fused_table_3_sub$NRMSE <- NA
    
    for (i in 1:length(fused_table_3_sub$metagenome_abundances)){
      sub_table <- fused_table_3_sub[i,]
      t <- sub_table$metagenome_abundances
      w <- sub_table$metric_abundances
      
      value <- (abs(w - t)^2)/((w+t)/2)
      fused_table_3_sub$NRMSE[i] <- value
    }
    
    #value_sub <- sqrt(mean(fused_table_3_sub$NRMSE))
    value_sub <- sqrt(1/length(fused_table_3_sub$NRMSE) * sum(fused_table_3_sub$NRMSE))
    
    value_sub_av <- sum(fused_table_3_sub$metagenome_abundances)/length(fused_table_3_sub$metagenome_abundances)
    
    data_tog <- data.frame(t(data.frame(c(paste(genus), value_sub,value_sub_av))))
    data_tog_2 <- rbind(data_tog_2, data_tog)
  }
  
  data_tog_2$metric <- paste(metric) 
  data_tog_3 <- rbind(data_tog_3, data_tog_2)
  
}

row.names(data_tog_3) <- NULL
colnames(data_tog_3) <- c("Genus", "NRMSE","Relative_abundance","Metric")

data_tog_3$Amplicon <- sapply(strsplit(data_tog_3$Metric, "_"),"[[",2)
data_tog_3$Method <- sapply(strsplit(data_tog_3$Metric, "_"),"[[",1)

data_tog_3$NRMSE <- as.numeric(data_tog_3$NRMSE)
data_tog_3$Relative_abundance <- as.numeric(data_tog_3$Relative_abundance)
data_tog_3$Group <- data_tog_3$Method
data_tog_3$Group[data_tog_3$Method %in% c("100","99","98","97","95")] <- "Qiime2-VSEARCH"
data_tog_3$Group[grep("Amplicon", data_tog_3$Method)] <- "Kraken2 - Amplicon"
data_tog_3$Group[grep("Fingerprint", data_tog_3$Method)] <- "Kraken2 - Fingerprint"

data_tog_3$Method_2 <- data_tog_3$Method
data_tog_3$Method_2 <- gsub("100","VS-100%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("99","VS-99%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("98","VS-98%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("97","VS-97%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("95","VS-95%",data_tog_3$Method_2)

data_tog_3$Group[grep("Fingerprint", data_tog_3$Method)] <- "Kraken2 - Fingerprint"
data_tog_3$Group[grep("Amplicon", data_tog_3$Method)] <- "Kraken2 - Amplicon"

data_tog_3$Method_2[grep("Kraken2:Fingerprint", data_tog_3$Method_2)] <- "KM-2-F"
data_tog_3$Method_2[grep("Kraken3:Fingerprint", data_tog_3$Method_2)] <- "KM-3-F"
data_tog_3$Method_2[grep("Kraken5:Fingerprint", data_tog_3$Method_2)] <- "KM-5-F"
data_tog_3$Method_2[grep("Kraken2:Amplicon", data_tog_3$Method_2)] <- "KM-2-A"
data_tog_3$Method_2[grep("Kraken3:Amplicon", data_tog_3$Method_2)] <- "KM-3-A"
data_tog_3$Method_2[grep("Kraken5:Amplicon", data_tog_3$Method_2)] <- "KM-5-A"

#Stats
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(metric in unique(data_tog_3$Amplicon)) {
  # Subset data for the current plant
  data_tog_3_sub <- data_tog_3[data_tog_3$Amplicon == paste(metric), ]
  
  # Perform ANOVA
  fitAnova <- aov(NRMSE ~ Method, data=data_tog_3_sub)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Method$Letters
  
  # Store results
  anova_results[[metric]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(metric) {
  data.frame(Metric = metric, Method = names(anova_results[[metric]]), Letters = anova_results[[metric]])
}))
# Order factors based on your original setup
anova_results$
  ltlbl_combined$Metric <- factor(ltlbl_combined$Metric, levels = c("V3V4","V5V7"))

ltlbl_combined$Method <- factor(ltlbl_combined$Method, levels = c("95","97","98","99","100","Kraken2:Amplicon","Kraken3:Amplicon",
                                                                  "Kraken5:Amplicon","Kraken2:Fingerprint","Kraken3:Fingerprint",
                                                                  "Kraken5:Fingerprint","Mothur","Salmon","SyFi"))
ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Metric, ltlbl_combined$Method), ]

data_tog_3$Method_2 <- factor(data_tog_3$Method_2, levels = c("VS-95%","VS-97%","VS-98%","VS-99%","VS-100%","KM-2-A","KM-3-A", "KM-5-A",
                                                              "KM-2-F","KM-3-F","KM-5-F","Mothur","Salmon","SyFi"))

data_tog_3$Amplicon <- gsub("V3V4","V3-V4", data_tog_3$Amplicon)
data_tog_3$Amplicon <- gsub("V5V7","V5-V7", data_tog_3$Amplicon)

plot_violin_genus <- ggplot(data_tog_3, aes(x = Method_2, y = NRMSE, colour = Group)) + 
  geom_violin(scale = "width") +
  geom_jitter(
    aes(size = Relative_abundance), 
    shape = 16, 
    position = position_jitter(0.2), 
    show.legend = TRUE
  ) +
  scale_color_manual(values =c("#d55e00","#cc79a7","#0072b2","#f0e442","#009e73","#8BB5E7") ) +
  ylab("NRMSE") + 
  xlab("Method") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, angle = 25, hjust = 1), 
    axis.title.y = element_text(size = 18), 
    axis.title.x = element_blank(), 
    axis.text.y = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 14), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 14)
  ) +
  ggtitle("Genus - Weighted accuracy") +
  guides(colour = guide_legend(title = "Method")) +
  guides(size = guide_legend(title = "Relative abundance")) +
  facet_wrap(~Amplicon, scales = "free_x") +
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(NRMSE)*1.05), show.legend=FALSE)

plot_violin_genus

###Figure 4b - species - NRMSE =====
#Create data frame from vector 
data_tog_3 <- data.frame()

#For loop that loads in each dataset of the vector
#It loads the same shotgun-metagenome-sequenced dataset every time
#And the clustering output that says which SynCom isolate amplicon sequences or fingerprints cluster together based on the implemented sequence identity clustering percentage
#For SyFi this percentage is set on 100%

#the loop will cluster the microbiome abundances of the dataset according to the clustering
#Similarly it will use the same clustering for the isolate counts in the shotgun metagenome-sequenced dataset
#Isolates with 0 counts in the loaded amplicon dataset or the shotgun metagenome-sequenced dataset are discarded

for (metric in vector) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory, "comparison_datasets/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert dataset for comparison (SyFi, or Qiime2-VSEARCH-generated table)
  otu_table = read.table(paste0(working_directory, "comparison_datasets/AtSC_", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep = ""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  taxonomy <- read.table(paste(working_directory, "species_clustering.txt",sep= ""), sep = "\t", row.names =1, header =T)
  
  #Reformatting table of interest 
  no_of_species <- length(na.omit(unique(taxonomy$primary_cluster)))
  
  #Create new table with the right dimensions
  new_table <- data.frame(matrix(NA, nrow = no_of_genera, ncol = length(row.names(samples_meta))+1))
  
  species_all <- na.omit(unique(taxonomy$primary_cluster))
  
  for (number in 1:no_of_species){
    number_2 <- species_all[number]
    new_table[number,1] <- paste(number_2)
  }
  
  column_names <- c("Species", paste(row.names(samples_meta),sep=""))
  colnames(new_table) <- column_names
  
  row.names(new_table) <- new_table$Species
  new_table <- new_table %>% dplyr::select(-Species)
  
  #Rename and reorder column names of otu_table
  colnames(otu_table) = gsub("_t3", "", colnames(otu_table))
  colnames(otu_table) = gsub("_t0", "", colnames(otu_table))
  colnames(otu_table) = gsub("_undil", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v5_v7", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v3_v4", "", colnames(otu_table))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_table_2 <- otu_table[column_names_order]
  otu_table <- otu_table_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:length(species_all)){
    species_pick <- species_all[row]
    
    isolates <- na.omit(row.names(taxonomy)[taxonomy$primary_cluster == paste(species_pick)])
    
    if (length(isolates) > 1){
      other_list <- row.names(otu_table) 
      the_isolate <- isolates[isolates %in% other_list]
    } else {
      the_isolate <- isolates
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_table[row.names(otu_table) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    new_table[row,] <- otu_table_3
  }
  
  #Create new table with the right dimensions for metagenome table
  meta_table <- data.frame(matrix(NA, nrow = no_of_species, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_species){
    number_2 <- species_all[number]
    meta_table[number,1] <- paste(number_2)
  }
  
  column_names <- c("Genus", paste(row.names(samples_meta),sep=""))
  colnames(meta_table) <- column_names
  
  row.names(meta_table) <- meta_table$Genus
  meta_table <- meta_table %>% dplyr::select(-Genus)
  
  #Reshuffle columns of the metagenome data
  colnames(otu_meta) = gsub("_undil", "", colnames(otu_meta))
  colnames(otu_meta) = gsub("_HL_", "_HL_orig_", colnames(otu_meta))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_meta_2 <- otu_meta[colnames(otu_meta) %in% column_names_order]
  otu_meta <- otu_meta_2[column_names_order]
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:length(species_all)){
    species_pick <- species_all[row]
    isolates <- na.omit(row.names(taxonomy)[taxonomy$primary_cluster == paste(species_pick)])
    
    if (length(isolates) > 1){
      other_list <- row.names(otu_meta) 
      the_isolate <- isolates[isolates %in% other_list]
    } else {
      the_isolate <- isolates
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_meta[row.names(otu_meta) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    meta_table[row,] <- otu_table_3
  }
  
  #Melt tables
  meta_table$Species <- row.names(meta_table)
  meta_table_2 <- melt(meta_table)
  new_table$Species <- row.names(new_table)
  new_table_2 <- melt(new_table)
  
  #Make relative abundances
  meta_table_3 <- Hank_the_normalizer(meta_table_2,"variable","value")
  new_table_3 <- Hank_the_normalizer(new_table_2,"variable","value")
  
  meta_table_3$ID <- paste(meta_table_3$variable,meta_table_3$Species,sep="_")
  new_table_3$ID <- paste(new_table_3$variable,new_table_3$Species,sep="_")
  
  fused_table <- data.frame(matrix(NA, nrow = length(meta_table_3$ID), ncol = 0))
  fused_table$V1 <- meta_table_3$Rel
  row.names(fused_table) <- meta_table_3$ID
  fused_table$V2 <- new_table_3$Rel[match(row.names(fused_table),new_table_3$ID)]
  colnames(fused_table) <- c("metagenome_abundances", "metric_abundances")
  
  fused_table_2 <- na.omit(fused_table[fused_table$metric_abundances != 0,])
  fused_table_3 <- na.omit(fused_table_2[fused_table_2$metagenome_abundances != 0,])
  fused_table_3$Species <- meta_table_3$Species[match(row.names(fused_table_3), meta_table_3$ID)]
  
  data_tog_2 <- data.frame()
  
  for (species in unique(fused_table_3$Species)){
    fused_table_3_sub <- fused_table_3[fused_table_3$Species == paste(species),]
    fused_table_3_sub$NRMSE <- NA
    
    for (i in 1:length(fused_table_3_sub$metagenome_abundances)){
      sub_table <- fused_table_3_sub[i,]
      t <- sub_table$metagenome_abundances
      w <- sub_table$metric_abundances
      
      value <- (abs(w - t)^2)/((w+t)/2)
      fused_table_3_sub$NRMSE[i] <- value
    }
    
    #value_sub <- sqrt(mean(fused_table_3_sub$NRMSE))
    value_sub <- sqrt(1/length(fused_table_3_sub$NRMSE) * sum(fused_table_3_sub$NRMSE))
    
    value_sub_av <- sum(fused_table_3_sub$metagenome_abundances)/length(fused_table_3_sub$metagenome_abundances)
    
    data_tog <- data.frame(t(data.frame(c(paste(genus), value_sub,value_sub_av))))
    data_tog_2 <- rbind(data_tog_2, data_tog)
  }
  
  data_tog_2$metric <- paste(metric) 
  data_tog_3 <- rbind(data_tog_3, data_tog_2)
  
}


row.names(data_tog_3) <- NULL
colnames(data_tog_3) <- c("Species", "NRMSE","Relative_abundance","Metric")

data_tog_3$Amplicon <- sapply(strsplit(data_tog_3$Metric, "_"),"[[",2)
data_tog_3$Method <- sapply(strsplit(data_tog_3$Metric, "_"),"[[",1)

data_tog_3$NRMSE <- as.numeric(data_tog_3$NRMSE)
data_tog_3$Relative_abundance <- as.numeric(data_tog_3$Relative_abundance)
data_tog_3$Group <- data_tog_3$Method
data_tog_3$Group[data_tog_3$Method %in% c("100","99","98","97","95")] <- "Qiime2-VSEARCH"
data_tog_3$Group[grep("Amplicon", data_tog_3$Method)] <- "Kraken2 - Amplicon"
data_tog_3$Group[grep("Fingerprint", data_tog_3$Method)] <- "Kraken2 - Fingerprint"

data_tog_3$Method_2 <- data_tog_3$Method
data_tog_3$Method_2 <- gsub("100","VS-100%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("99","VS-99%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("98","VS-98%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("97","VS-97%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("95","VS-95%",data_tog_3$Method_2)

data_tog_3$Group[grep("Fingerprint", data_tog_3$Method)] <- "Kraken2 - Fingerprint"
data_tog_3$Group[grep("Amplicon", data_tog_3$Method)] <- "Kraken2 - Amplicon"

data_tog_3$Method_2[grep("Kraken2:Fingerprint", data_tog_3$Method_2)] <- "KM-2-F"
data_tog_3$Method_2[grep("Kraken3:Fingerprint", data_tog_3$Method_2)] <- "KM-3-F"
data_tog_3$Method_2[grep("Kraken5:Fingerprint", data_tog_3$Method_2)] <- "KM-5-F"
data_tog_3$Method_2[grep("Kraken2:Amplicon", data_tog_3$Method_2)] <- "KM-2-A"
data_tog_3$Method_2[grep("Kraken3:Amplicon", data_tog_3$Method_2)] <- "KM-3-A"
data_tog_3$Method_2[grep("Kraken5:Amplicon", data_tog_3$Method_2)] <- "KM-5-A"

#Stats
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(metric in unique(data_tog_3$Amplicon)) {
  # Subset data for the current plant
  data_tog_3_sub <- data_tog_3[data_tog_3$Amplicon == paste(metric), ]
  
  # Perform ANOVA
  fitAnova <- aov(NRMSE ~ Method, data=data_tog_3_sub)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Method$Letters
  
  # Store results
  anova_results[[metric]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(metric) {
  data.frame(Metric = metric, Method = names(anova_results[[metric]]), Letters = anova_results[[metric]])
}))
# Order factors based on your original setup
ltlbl_combined$Metric <- factor(ltlbl_combined$Metric, levels = c("V3V4","V5V7"))

ltlbl_combined$Method <- factor(ltlbl_combined$Method, levels = c("95","97","98","99","100","Kraken2:Amplicon","Kraken3:Amplicon",
                                                                  "Kraken5:Amplicon","Kraken2:Fingerprint","Kraken3:Fingerprint",
                                                                  "Kraken5:Fingerprint","Mothur","Salmon","SyFi"))
ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Metric, ltlbl_combined$Method), ]

data_tog_3$Method_2 <- factor(data_tog_3$Method_2, levels = c("VS-95%","VS-97%","VS-98%","VS-99%","VS-100%","KM-2-A","KM-3-A", "KM-5-A",
                                                              "KM-2-F","KM-3-F","KM-5-F","Mothur","Salmon","SyFi"))

data_tog_3$Amplicon <- gsub("V3V4","V3-V4", data_tog_3$Amplicon)
data_tog_3$Amplicon <- gsub("V5V7","V5-V7", data_tog_3$Amplicon)

plot_violin_species <- ggplot(data_tog_3, aes(x = Method_2, y = NRMSE, colour = Group)) + 
  geom_violin(scale = "width") +
  geom_jitter(
    aes(size = Relative_abundance), 
    shape = 16, 
    position = position_jitter(0.2), 
    show.legend = TRUE
  ) +
  scale_color_manual(values =c("#d55e00","#cc79a7","#0072b2","#f0e442","#009e73","#8BB5E7") ) +
  ylab("NRMSE") + 
  xlab("Method") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, angle = 25, hjust = 1), 
    axis.title.y = element_text(size = 18), 
    axis.title.x = element_blank(), 
    axis.text.y = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 14), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 14)
  ) +
  ggtitle("Species - Weighted accuracy") +
  guides(colour = guide_legend(title = "Method")) +
  guides(size = guide_legend(title = "Relative abundance")) +
  facet_wrap(~Amplicon, scales = "free_x") +
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(NRMSE)*1.05), show.legend=FALSE)

plot_violin_species

###Figure 4c - strain - NRMSE =====
#Create data frame from vector 
data_tog_3 <- data.frame()

#For loop that loads in each dataset of the vector
#It loads the same shotgun-metagenome-sequenced dataset every time
#And the clustering output that says which SynCom isolate amplicon sequences or fingerprints cluster together based on the implemented sequence identity clustering percentage
#For SyFi this percentage is set on 100%

#the loop will cluster the microbiome abundances of the dataset according to the clustering
#Similarly it will use the same clustering for the isolate counts in the shotgun metagenome-sequenced dataset
#Isolates with 0 counts in the loaded amplicon dataset or the shotgun metagenome-sequenced dataset are discarded

for (metric in vector) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory, "comparison_datasets/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert dataset for comparison (SyFi, or Qiime2-VSEARCH-generated table)
  otu_table = read.table(paste0(working_directory, "comparison_datasets/AtSC_", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep = ""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  #Insert clusters file (which isolate is in which cluster based on amplicon sequence or fingerprint and at different sequence identities)
  if (metric %in% c("Salmon_V3V4", "Kraken2:Amplicon_V3V4","Kraken3:Amplicon_V3V4","Kraken5:Amplicon_V3V4")){
    metric_2 <- "100_V3V4"
  } else if (metric %in% c("Salmon_V5V7", "Kraken2:Amplicon_V5V7","Kraken3:Amplicon_V5V7","Kraken5:Amplicon_V5V7")){
    metric_2 <- "100_V5V7"
  } else if (metric %in% c("Kraken2:Fingerprint_V3V4","Kraken3:Fingerprint_V3V4","Kraken5:Fingerprint_V3V4")) {
    metric_2 <- "SyFi_V3V4"
  } else if (metric %in% c("Kraken2:Fingerprint_V5V7","Kraken3:Fingerprint_V5V7","Kraken5:Fingerprint_V5V7")) {
    metric_2 <- "SyFi_V5V7"
  } else {
    metric_2 <- metric
  }
  clusters = read.csv(paste0(working_directory,"vsearch_clusters/collapse_",metric_2, ".txt",sep=""), header=F, sep="\t")
  
  #Reformatting table of interest 
  no_of_isolates <- length(row.names(clusters))
  
  #Create new table with the right dimensions
  new_table <- data.frame(matrix(NA, nrow = no_of_isolates, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_isolates){
    new_table[number,1] <- paste("Cluster_", number,sep="")
  }
  
  column_names <- c("cluster", paste(row.names(samples_meta),sep=""))
  colnames(new_table) <- column_names
  
  row.names(new_table) <- new_table$cluster
  new_table <- new_table %>% dplyr::select(-cluster)
  
  #Rename and reorder column names of otu_table
  colnames(otu_table) = gsub("_t3", "", colnames(otu_table))
  colnames(otu_table) = gsub("_t0", "", colnames(otu_table))
  colnames(otu_table) = gsub("_undil", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v5_v7", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v3_v4", "", colnames(otu_table))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_table_2 <- otu_table[column_names_order]
  otu_table <- otu_table_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:nrow(clusters)){
    isolate_vector <- as.character(clusters[row,])
    isolate_vector <- isolate_vector[!is.na(isolate_vector)]
    isolate_vector <- isolate_vector[!isolate_vector == ""]
    
    if (length(isolate_vector) > 1){
      other_list <- row.names(otu_table) 
      the_isolate <- isolate_vector[isolate_vector %in% other_list]
    } else {
      the_isolate <- isolate_vector
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_table[row.names(otu_table) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    new_table[row,] <- otu_table_3
  }
  
  #Create new table with the right dimensions for metagenome table
  meta_table <- data.frame(matrix(NA, nrow = no_of_isolates, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_isolates){
    meta_table[number,1] <- paste("Cluster_", number,sep="")
  }
  
  column_names <- c("cluster", paste(row.names(samples_meta),sep=""))
  colnames(meta_table) <- column_names
  
  row.names(meta_table) <- meta_table$cluster
  meta_table <- meta_table %>% dplyr::select(-cluster)
  
  #Reshuffle columns of the metagenome data
  colnames(otu_meta) = gsub("_undil", "", colnames(otu_meta))
  colnames(otu_meta) = gsub("_HL_", "_HL_orig_", colnames(otu_meta))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_meta_2 <- otu_meta[colnames(otu_meta) %in% column_names_order]
  otu_meta <- otu_meta_2[column_names_order]
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:nrow(clusters)){
    isolate_vector <- as.character(clusters[row,])
    other_list <- row.names(otu_meta) 
    the_isolate <- isolate_vector[isolate_vector %in% other_list]
    the_isolate <- the_isolate[!is.na(the_isolate)]
    the_isolate <- the_isolate[!the_isolate == ""]
    
    #subset otu_table of choice
    otu_table_2 <- otu_meta[row.names(otu_meta) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    meta_table[row,] <- otu_table_3
  }
  
  #Melt tables
  meta_table$Cluster <- row.names(meta_table)
  meta_table_2 <- melt(meta_table)
  new_table$Cluster <- row.names(new_table)
  new_table_2 <- melt(new_table)
  length(new_table[new_table[,2] > 0,2])
  #Make relative abundances
  meta_table_3 <- Hank_the_normalizer(meta_table_2,"variable","value")
  new_table_3 <- Hank_the_normalizer(new_table_2,"variable","value")
  
  meta_table_3$ID <- paste(meta_table_3$variable,meta_table_3$Cluster,sep="_")
  new_table_3$ID <- paste(new_table_3$variable,new_table_3$Cluster,sep="_")
  
  fused_table <- data.frame(matrix(NA, nrow = length(meta_table_3$ID), ncol = 0))
  fused_table$V1 <- meta_table_3$Rel
  row.names(fused_table) <- meta_table_3$ID
  fused_table$V2 <- new_table_3$Rel[match(row.names(fused_table),new_table_3$ID)]
  colnames(fused_table) <- c("metagenome_abundances", "metric_abundances")
  
  fused_table_2 <- na.omit(fused_table[fused_table$metric_abundances != 0,])
  fused_table_3 <- na.omit(fused_table_2[fused_table_2$metagenome_abundances != 0,])
  fused_table_3$Strain <- meta_table_3$Cluster[match(row.names(fused_table_3), meta_table_3$ID)]
  
  data_tog_2 <- data.frame()
  
  for (strain in unique(fused_table_3$Strain)){
    fused_table_3_sub <- fused_table_3[fused_table_3$Strain == paste(strain),]
    fused_table_3_sub$NRMSE <- NA
    
    for (i in 1:length(fused_table_3_sub$metagenome_abundances)){
      sub_table <- fused_table_3_sub[i,]
      t <- sub_table$metagenome_abundances
      w <- sub_table$metric_abundances
      
      value <- (abs(w - t)^2)/((w+t)/2)
      fused_table_3_sub$NRMSE[i] <- value
    }
    
    #value_sub <- sqrt(mean(fused_table_3_sub$NRMSE))
    value_sub <- sqrt(1/length(fused_table_3_sub$NRMSE) * sum(fused_table_3_sub$NRMSE))
    
    value_sub_av <- sum(fused_table_3_sub$metagenome_abundances)/length(fused_table_3_sub$metagenome_abundances)
    
    data_tog <- data.frame(t(data.frame(c(paste(genus), value_sub,value_sub_av))))
    data_tog_2 <- rbind(data_tog_2, data_tog)
  }
  
  data_tog_2$metric <- paste(metric) 
  data_tog_3 <- rbind(data_tog_3, data_tog_2)
  
}

row.names(data_tog_3) <- NULL
colnames(data_tog_3) <- c("Cluster", "NRMSE","Relative_abundance","Metric")

data_tog_3$Amplicon <- sapply(strsplit(data_tog_3$Metric, "_"),"[[",2)
data_tog_3$Method <- sapply(strsplit(data_tog_3$Metric, "_"),"[[",1)

data_tog_3$NRMSE <- as.numeric(data_tog_3$NRMSE)
data_tog_3$Relative_abundance <- as.numeric(data_tog_3$Relative_abundance)
data_tog_3$Group <- data_tog_3$Method
data_tog_3$Group[data_tog_3$Method %in% c("100","99","98","97","95")] <- "Qiime2-VSEARCH"
data_tog_3$Group[grep("Amplicon", data_tog_3$Method)] <- "Kraken2 - Amplicon"
data_tog_3$Group[grep("Fingerprint", data_tog_3$Method)] <- "Kraken2 - Fingerprint"

data_tog_3$Method_2 <- data_tog_3$Method
data_tog_3$Method_2 <- gsub("100","VS-100%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("99","VS-99%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("98","VS-98%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("97","VS-97%",data_tog_3$Method_2)
data_tog_3$Method_2 <- gsub("95","VS-95%",data_tog_3$Method_2)

data_tog_3$Group[grep("Fingerprint", data_tog_3$Method)] <- "Kraken2 - Fingerprint"
data_tog_3$Group[grep("Amplicon", data_tog_3$Method)] <- "Kraken2 - Amplicon"

data_tog_3$Method_2[grep("Kraken2:Fingerprint", data_tog_3$Method_2)] <- "KM-2-F"
data_tog_3$Method_2[grep("Kraken3:Fingerprint", data_tog_3$Method_2)] <- "KM-3-F"
data_tog_3$Method_2[grep("Kraken5:Fingerprint", data_tog_3$Method_2)] <- "KM-5-F"
data_tog_3$Method_2[grep("Kraken2:Amplicon", data_tog_3$Method_2)] <- "KM-2-A"
data_tog_3$Method_2[grep("Kraken3:Amplicon", data_tog_3$Method_2)] <- "KM-3-A"
data_tog_3$Method_2[grep("Kraken5:Amplicon", data_tog_3$Method_2)] <- "KM-5-A"

#Stats
anova_results <- list()
# Loop through each plant type to perform ANOVA, Tukey's test, and get letters
for(metric in unique(data_tog_3$Amplicon)) {
  # Subset data for the current plant
  data_tog_3_sub <- data_tog_3[data_tog_3$Amplicon == paste(metric), ]
  
  # Perform ANOVA
  fitAnova <- aov(NRMSE ~ Method, data=data_tog_3_sub)
  
  # Perform Tukey's post-hoc test
  Tukey <- TukeyHSD(fitAnova)
  
  # Get letters
  letters_anova <- multcompView::multcompLetters4(fitAnova, Tukey)$Method$Letters
  
  # Store results
  anova_results[[metric]] <- letters_anova
}
# Combine results into a data frame for plotting
ltlbl_combined <- do.call(rbind, lapply(names(anova_results), function(metric) {
  data.frame(Metric = metric, Method = names(anova_results[[metric]]), Letters = anova_results[[metric]])
}))
# Order factors based on your original setup
ltlbl_combined$Metric <- factor(ltlbl_combined$Metric, levels = c("V3V4","V5V7"))

ltlbl_combined$Method <- factor(ltlbl_combined$Method, levels = c("95","97","98","99","100","Kraken2:Amplicon","Kraken3:Amplicon",
                                                                  "Kraken5:Amplicon","Kraken2:Fingerprint","Kraken3:Fingerprint",
                                                                  "Kraken5:Fingerprint","Mothur","Salmon","SyFi"))
ltlbl_combined <- ltlbl_combined[order(ltlbl_combined$Metric, ltlbl_combined$Method), ]

data_tog_3$Method_2 <- factor(data_tog_3$Method_2, levels = c("VS-95%","VS-97%","VS-98%","VS-99%","VS-100%","KM-2-A","KM-3-A", "KM-5-A",
                                                              "KM-2-F","KM-3-F","KM-5-F","Mothur","Salmon","SyFi"))

data_tog_3$Amplicon <- gsub("V3V4","V3-V4", data_tog_3$Amplicon)
data_tog_3$Amplicon <- gsub("V5V7","V5-V7", data_tog_3$Amplicon)

plot_violin_strain <- ggplot(data_tog_3, aes(x = Method_2, y = NRMSE, colour = Group)) + 
  geom_violin(scale = "width") +
  geom_jitter(
    aes(size = Relative_abundance), 
    shape = 16, 
    position = position_jitter(0.2), 
    show.legend = TRUE
  ) +
  scale_color_manual(values =c("#d55e00","#cc79a7","#0072b2","#f0e442","#009e73","#8BB5E7") ) +
  ylab("NRMSE") + 
  xlab("Method") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, angle = 25, hjust = 1), 
    axis.title.y = element_text(size = 18), 
    axis.title.x = element_blank(), 
    axis.text.y = element_text(size = 14), 
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 14), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 14)
  ) +
  ggtitle("Strain - Weighted accuracy") +
  guides(colour = guide_legend(title = "Method")) +
  guides(size = guide_legend(title = "Relative abundance")) +
  facet_wrap(~Amplicon, scales = "free_x") +
  stat_summary(geom = 'text', label = ltlbl_combined$Letters, fun.y = max, aes(y = max(NRMSE)*1.05), show.legend=FALSE)

plot_violin_strain

###Figure 4 all together ====
together_plot <- ggarrange(plot_violin_genus, plot_violin_species, plot_violin_strain, ncol = 1, nrow = 3,labels = c("A", "B", "C"), common.legend = T)

pdf(paste(results.dir,"Figure_4.pdf", sep=""), width=15, height=25)
print(together_plot)
dev.off()

###Figure 5 - SyFi performance summary =====
Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

#Taxonomy file
taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

#Vector of each dataset - V3V4 or V5V7 amplicons/fingerprints - Qiime2-VSEARCH clustering at 95%, 97%, 98%, 99%, 100% and SyFi (Salmon pseudoalignment of amplicon reads to SyFi-generated fingerprints)
vector <- c("95_V3V4", "97_V3V4", "98_V3V4", "99_V3V4", "100_V3V4", 
            "Salmon_V3V4", "SyFi_V3V4","Mothur_V3V4",
            "Kraken2:Amplicon_V3V4","Kraken3:Amplicon_V3V4","Kraken5:Amplicon_V3V4",
            "Kraken2:Fingerprint_V3V4","Kraken3:Fingerprint_V3V4","Kraken5:Fingerprint_V3V4",
            "95_V5V7", "97_V5V7", "98_V5V7", "99_V5V7", "100_V5V7", 
            "Salmon_V5V7","SyFi_V5V7", "Mothur_V5V7",
            "Kraken2:Amplicon_V5V7","Kraken3:Amplicon_V5V7","Kraken5:Amplicon_V5V7",
            "Kraken2:Fingerprint_V5V7","Kraken3:Fingerprint_V5V7","Kraken5:Fingerprint_V5V7")
#Create data frame from vector 
correlations <- as.data.frame(vector)
correlations$NRMSE <- NA

#For loop that loads in each dataset of the vector
#It loads the same shotgun-metagenome-sequenced dataset every time
#And the clustering output that says which SynCom isolate amplicon sequences or fingerprints cluster together based on the implemented sequence identity clustering percentage
#For SyFi this percentage is set on 100%

#the loop will cluster the microbiome abundances of the dataset according to the clustering
#Similarly it will use the same clustering for the isolate counts in the shotgun metagenome-sequenced dataset
#Isolates with 0 counts in the loaded amplicon dataset or the shotgun metagenome-sequenced dataset are discarded

for (metric in vector) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory, "comparison_datasets/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert dataset for comparison (SyFi, or Qiime2-VSEARCH-generated table)
  otu_table = read.table(paste0(working_directory, "comparison_datasets/AtSC_", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep = ""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  #Insert clusters file (which isolate is in which cluster based on amplicon sequence or fingerprint and at different sequence identities)
  if (metric %in% c("Salmon_V3V4", "Kraken2:Amplicon_V3V4","Kraken3:Amplicon_V3V4","Kraken5:Amplicon_V3V4")){
    metric_2 <- "100_V3V4"
  } else if (metric %in% c("Salmon_V5V7", "Kraken2:Amplicon_V5V7","Kraken3:Amplicon_V5V7","Kraken5:Amplicon_V5V7")){
    metric_2 <- "100_V5V7"
  } else if (metric %in% c("Kraken2:Fingerprint_V3V4","Kraken3:Fingerprint_V3V4","Kraken5:Fingerprint_V3V4")) {
    metric_2 <- "SyFi_V3V4"
  } else if (metric %in% c("Kraken2:Fingerprint_V5V7","Kraken3:Fingerprint_V5V7","Kraken5:Fingerprint_V5V7")) {
    metric_2 <- "SyFi_V5V7"
  } else {
    metric_2 <- metric
  }
  clusters = read.csv(paste0(working_directory,"vsearch_clusters/collapse_",metric_2, ".txt",sep=""), header=F, sep="\t")
  
  #Reformatting table of interest 
  no_of_isolates <- length(row.names(clusters))
  
  #Create new table with the right dimensions
  new_table <- data.frame(matrix(NA, nrow = no_of_isolates, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_isolates){
    new_table[number,1] <- paste("Cluster_", number,sep="")
  }
  
  column_names <- c("cluster", paste(row.names(samples_meta),sep=""))
  colnames(new_table) <- column_names
  
  row.names(new_table) <- new_table$cluster
  new_table <- new_table %>% dplyr::select(-cluster)
  
  #Rename and reorder column names of otu_table
  colnames(otu_table) = gsub("_t3", "", colnames(otu_table))
  colnames(otu_table) = gsub("_t0", "", colnames(otu_table))
  colnames(otu_table) = gsub("_undil", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v5_v7", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v3_v4", "", colnames(otu_table))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_table_2 <- otu_table[column_names_order]
  otu_table <- otu_table_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:nrow(clusters)){
    isolate_vector <- as.character(clusters[row,])
    isolate_vector <- isolate_vector[!is.na(isolate_vector)]
    isolate_vector <- isolate_vector[!isolate_vector == ""]
    
    if (length(isolate_vector) > 1){
      other_list <- row.names(otu_table) 
      the_isolate <- isolate_vector[isolate_vector %in% other_list]
    } else {
      the_isolate <- isolate_vector
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_table[row.names(otu_table) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    new_table[row,] <- otu_table_3
  }
  
  #Create new table with the right dimensions for metagenome table
  meta_table <- data.frame(matrix(NA, nrow = no_of_isolates, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_isolates){
    meta_table[number,1] <- paste("Cluster_", number,sep="")
  }
  
  column_names <- c("cluster", paste(row.names(samples_meta),sep=""))
  colnames(meta_table) <- column_names
  
  row.names(meta_table) <- meta_table$cluster
  meta_table <- meta_table %>% dplyr::select(-cluster)
  
  #Reshuffle columns of the metagenome data
  colnames(otu_meta) = gsub("_undil", "", colnames(otu_meta))
  colnames(otu_meta) = gsub("_HL_", "_HL_orig_", colnames(otu_meta))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_meta_2 <- otu_meta[colnames(otu_meta) %in% column_names_order]
  otu_meta <- otu_meta_2[column_names_order]
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:nrow(clusters)){
    isolate_vector <- as.character(clusters[row,])
    other_list <- row.names(otu_meta) 
    the_isolate <- isolate_vector[isolate_vector %in% other_list]
    the_isolate <- the_isolate[!is.na(the_isolate)]
    the_isolate <- the_isolate[!the_isolate == ""]
    
    #subset otu_table of choice
    otu_table_2 <- otu_meta[row.names(otu_meta) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    meta_table[row,] <- otu_table_3
  }
  
  #Melt tables
  meta_table$Cluster <- row.names(meta_table)
  meta_table_2 <- melt(meta_table)
  new_table$Cluster <- row.names(new_table)
  new_table_2 <- melt(new_table)
  length(new_table[new_table[,2] > 0,2])
  #Make relative abundances
  meta_table_3 <- Hank_the_normalizer(meta_table_2,"variable","value")
  new_table_3 <- Hank_the_normalizer(new_table_2,"variable","value")
  
  meta_table_3$ID <- paste(meta_table_3$variable,meta_table_3$Cluster,sep="_")
  new_table_3$ID <- paste(new_table_3$variable,new_table_3$Cluster,sep="_")
  
  fused_table <- data.frame(matrix(NA, nrow = length(meta_table_3$ID), ncol = 0))
  fused_table$V1 <- meta_table_3$Rel
  row.names(fused_table) <- meta_table_3$ID
  fused_table$V2 <- new_table_3$Rel[match(row.names(fused_table),new_table_3$ID)]
  colnames(fused_table) <- c("metagenome_abundances", "metric_abundances")
  
  fused_table_2 <- na.omit(fused_table[fused_table$metric_abundances != 0,])
  fused_table_3 <- na.omit(fused_table_2[fused_table_2$metagenome_abundances != 0,])
  fused_table_3$Strain <- meta_table_3$Cluster[match(row.names(fused_table_3), meta_table_3$ID)]
  
  fused_table_3$NRMSE <- NA
  
  for (i in 1:length(fused_table_3$metagenome_abundances)){
    sub_table <- fused_table_3[i,]
    t <- sub_table$metagenome_abundances
    w <- sub_table$metric_abundances
    
    value <- (abs(w - t)^2)/((w+t)/2)
    fused_table_3$NRMSE[i] <- value
  }
  
  value_sub <- sqrt(1/length(fused_table_3$NRMSE) * sum(fused_table_3$NRMSE))
  
  correlations$NRMSE[correlations$vector == paste0(metric)] <- value_sub  
}

correlations$NRMSE <- 1/correlations$NRMSE

#Vsearch clustering tables
correlations$Combination <- correlations$vector
correlations$Amplicon <- sapply(strsplit(correlations$vector, "_"),"[[",2)
No_clusters = read.table(paste(working_directory,"No_clusters.txt", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)
No_clusters$vector <- gsub("V3-V4","V3V4", No_clusters$vector)
No_clusters$vector <- gsub("V5-V7","V5V7", No_clusters$vector)

correlations$vector <- gsub("SyFi_V3V4", "SyFi", correlations$vector)
correlations$vector <- gsub("SyFi_V5V7", "SyFi", correlations$vector)

correlations$vector <- gsub("_V3V4", "%", correlations$vector)
correlations$vector <- gsub("_V5V7", "%", correlations$vector)
correlations$vector <- gsub("9", "VS-9", correlations$vector)
correlations$vector <- gsub("9VS-", "9", correlations$vector)
correlations$vector <- gsub("100", "VS-100", correlations$vector)
correlations$vector <- gsub("Salmon%", "Salmon", correlations$vector)
correlations$vector <- gsub(":Amplicon%", "-A", correlations$vector)
correlations$vector <- gsub(":Fingerprint%", "-F", correlations$vector)
correlations$vector <- gsub("Genome%", "Genome", correlations$vector)
correlations$vector <- gsub("Mothur%", "Mothur", correlations$vector)
correlations$vector <- gsub("Kraken", "KM-", correlations$vector)

correlations$No_of_Clusters <- as.numeric(No_clusters$Clusters[match(correlations$Combination, No_clusters$vector)])

correlations$Amplicon <- gsub("V3V4", "V3-V4", correlations$Amplicon)
correlations$Amplicon <- gsub("V5V7", "V5-V7", correlations$Amplicon)
correlations$Amplicon[correlations$Amplicon == "SyFi"] <- c("V3-V4","V5-V7")

correlations$Group <- correlations$vector 
correlations$Group[correlations$vector %in% c("VS-100%","VS-99%","VS-98%","VS-97%","VS-95%")] <- "Qiime2-VSEARCH"
correlations$Group[grep("-A", correlations$Group)] <- "Kraken2 - Amplicon"
correlations$Group[grep("-F", correlations$Group)] <- "Kraken2 - Fingerprint"

correlations$vector <- factor(correlations$vector, levels = c("VS-95%","VS-97%","VS-98%","VS-99%","VS-100%","KM-2-A","KM-3-A",
                                                              "KM-5-A","KM-2-F","KM-3-F",
                                                              "KM-5-F","Mothur","Salmon","SyFi"))

plot_occ <- ggscatter(correlations, x = "NRMSE", y = "No_of_Clusters", color = "Group") + 
  ggtitle("SyFi weighted accuracy") + 
  theme_classic() +
  geom_text_repel(aes(label = vector), size = 4) + 
  scale_color_manual(values =c("#d55e00","#cc79a7","#0072b2","#f0e442","#009e73","#8BB5E7") ) +
  ylab("Number of clusters") + 
  xlab("1/NRMSE") +
  labs(color = "Metric") +
  theme(axis.text.x = element_text(size = 10), axis.title = element_text(size = 12), 
        axis.text.y = element_text(size=10), legend.title = element_text(size=12), 
        legend.text = element_text(size=10)) +
  facet_wrap(~Amplicon, scales = "free_x") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

plot_occ

pdf(paste(results.dir,"Figure_5.pdf", sep=""), width=14, height=5)
print(plot_occ)
dev.off()


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



###Figure S2 - Genome assembly effect on SyFi output =====

#Figure S2A - Genome completeness plot
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
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size =12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18)) +  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("Genome completeness (%)") + 
  xlab("Fingerprint")
plot_comp_2 <- plot_comp + stat_compare_means(comparisons = list(c("Present", "Absent")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
plot_comp_2

#Figure S2B and C - Genome coverage and GC content
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
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
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
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size =12), axis.title = element_text(size = 14), axis.text.y = element_text(size=12), legend.title = element_text(size=14), legend.text = element_text(size=12), plot.title = element_text(size=18)) +  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("GC content (%)") + 
  xlab("Fingerprint")
plot_GC_2 <- plot_GC + stat_compare_means(comparisons = list(c("Present", "Absent")), method = "t.test", vjust = 1.5) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
plot_GC_2

#Figure S2D - L50
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
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
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

pdf(paste(results.dir,"Figure_S2.pdf", sep=""), width=10, height=8)
print(all)
dev.off()

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



###Figure S4 - 16S copy number and haplotypes comparison - SyFi vs complete genomes =====

#This script generates the figures showing the 16S copy number and haplotype distribution from closed bacterial genomes
table_16S <- read.table(paste(working_directory, "complete_genomes/16S_copy_number.tsv", sep = ""), header =T, sep = "\t")
table_16S_hap <- read.table(paste(working_directory,"complete_genomes/16S_haplotypes.tsv", sep =""), header =T,  sep = "\t")
table_16S$Haplotype <- table_16S_hap$X16S_haplotypes[match(table_16S$isolate, table_16S_hap$isolate)]

taxonomy <- read.table(paste(working_directory, "complete_genomes/taxonomy.tsv", sep =""), header =T,  sep = "\t")

table_16S_2 <- table_16S[table_16S$X16S_copy_number != 0,]

#Figure S4A
plot_16S_CN <- ggplot(table_16S_2, aes(x=X16S_copy_number)) +
  geom_bar() +
  ggtitle("16S copy number - complete genomes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S copy number") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_CN

#Figure S4B
plot_16S_hap <- ggplot(table_16S_2, aes(x=Haplotype)) +
  geom_bar() +
  ggtitle("16S haplotypes - complete genomes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Number of isolates") + 
  xlab("16S haplotypes") +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_16S_hap

#Family plot - Figure S4C
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
  scale_color_manual(values =c("#D81B60") ) +
  guides(colour = guide_legend(title = "16S")) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")
plot_family_2

figure <- ggarrange(plot_16S_CN, plot_16S_hap,
                    labels = c("A", "B" ),font.label = list(size = 24),ncol = 2, nrow = 1, widths = 15, heights = 8)

figure_family <- ggarrange(plot_family_2,
                           labels = c("C" ),font.label = list(size = 24),ncol = 1, nrow = 1, widths = 15, heights = 8)
figure_combined <- ggarrange(figure,figure_family, ncol =1, nrow =2, widths = 15, heights = c(1,2))
figure_combined

#This produces figures S4A, B, and C. 
pdf(paste(results.dir,"Figures_S4ABC.pdf", sep=""), width=13, height=13)
print(figure_combined)
dev.off()

#Figures 4D and E are generated using tidyheatmaps with the commands below
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

#Figure S4D
Heatmap <- tidy_heatmap(table_16_1,rows = "Family", columns = "Database", values = "value", annotation_col = c(Database, variable),scale = "none", gaps_col = variable, fontsize_row = 6, main = "SyFi vs NCBI")

#Figure S4E
Heatmap_2 <- tidy_heatmap(table_16_2,rows = "Family", columns = "Database", values = "value", annotation_col = c(Database, variable),scale = "none", gaps_col = variable, fontsize_row = 6, main = "SyFi vs NCBI")

ggsave("Figure_S4D.png", plot = Heatmap, width = 12, height = 10, units = "cm", device = "png", dpi = 300, path = results.dir)
ggsave("Figure_S4E.png", plot = Heatmap_2, width = 12, height = 10, units = "cm", device = "png", dpi = 300, path = results.dir)

###Figure S6 - SyFi Accuracy with different Salmon minscorefraction settings=====

Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}
vector_2 <- c("65_V3V4", "75_V3V4", "85_V3V4", "95_V3V4","65_V5V7", "75_V5V7", "85_V5V7", "95_V5V7")

#Create data frame from vector 
correlations <- as.data.frame(vector_2)
correlations$NRMSE <- NA

for (metric in vector_2) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory,"minscorefraction/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert table for comparison
  otu_table = read.table(paste0(working_directory,"minscorefraction/", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep = ""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  #Insert clusters file
  metric_2 <- strsplit(paste(metric), "_")[[1]][2]
  clusters = read.csv(paste0(working_directory,"vsearch_clusters/collapse_SyFi_",metric_2, ".txt",sep=""), header=F, sep="\t")
  
  #Reformatting table of interest 
  no_of_isolates <- length(row.names(clusters))
  
  #Create new table with the right dimensions
  new_table <- data.frame(matrix(NA, nrow = no_of_isolates, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_isolates){
    new_table[number,1] <- paste("Cluster_", number,sep="")
  }
  
  column_names <- c("cluster", paste(row.names(samples_meta),sep=""))
  colnames(new_table) <- column_names
  
  row.names(new_table) <- new_table$cluster
  new_table <- new_table %>% dplyr::select(-cluster)
  
  #Rename and reorder column names of otu_table
  colnames(otu_table) = gsub("_t3", "", colnames(otu_table))
  colnames(otu_table) = gsub("_t0", "", colnames(otu_table))
  colnames(otu_table) = gsub("_undil", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v5_v7", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v3_v4", "", colnames(otu_table))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_table_2 <- otu_table[column_names_order]
  otu_table <- otu_table_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:nrow(clusters)){
    isolate_vector <- as.character(clusters[row,])
    isolate_vector <- isolate_vector[!is.na(isolate_vector)]
    isolate_vector <- isolate_vector[!isolate_vector == ""]
    
    if (length(isolate_vector) > 1){
      other_list <- row.names(otu_table) 
      the_isolate <- isolate_vector[isolate_vector %in% other_list]
    } else {
      the_isolate <- isolate_vector
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_table[row.names(otu_table) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    new_table[row,] <- otu_table_3
  }
  
  #Create new table with the right dimensions for metagenome table
  meta_table <- data.frame(matrix(NA, nrow = no_of_isolates, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_isolates){
    meta_table[number,1] <- paste("Cluster_", number,sep="")
  }
  
  column_names <- c("cluster", paste(row.names(samples_meta),sep=""))
  colnames(meta_table) <- column_names
  
  row.names(meta_table) <- meta_table$cluster
  meta_table <- meta_table %>% dplyr::select(-cluster)
  
  #Reshuffle columns of the metagenome data
  colnames(otu_meta) = gsub("_undil", "", colnames(otu_meta))
  colnames(otu_meta) = gsub("_t0", "", colnames(otu_meta))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_meta_2 <- otu_meta[colnames(otu_meta) %in% column_names_order]
  otu_meta <- otu_meta_2[column_names_order]
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:nrow(clusters)){
    isolate_vector <- as.character(clusters[row,])
    other_list <- row.names(otu_meta) 
    the_isolate <- isolate_vector[isolate_vector %in% other_list]
    the_isolate <- the_isolate[!is.na(the_isolate)]
    the_isolate <- the_isolate[!the_isolate == ""]
    
    #subset otu_table of choice
    otu_table_2 <- otu_meta[row.names(otu_meta) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    meta_table[row,] <- otu_table_3
  }
  
  #Melt tables
  meta_table$Cluster <- row.names(meta_table)
  meta_table_2 <- melt(meta_table)
  new_table$Cluster <- row.names(new_table)
  new_table_2 <- melt(new_table)
  length(new_table[new_table[,2] > 0,2])
  #Make relative abundances
  meta_table_3 <- Hank_the_normalizer(meta_table_2,"variable","value")
  new_table_3 <- Hank_the_normalizer(new_table_2,"variable","value")
  
  meta_table_3$ID <- paste(meta_table_3$variable,meta_table_3$Cluster,sep="_")
  new_table_3$ID <- paste(new_table_3$variable,new_table_3$Cluster,sep="_")
  
  fused_table <- data.frame(matrix(NA, nrow = length(meta_table_3$ID), ncol = 0))
  fused_table$V1 <- meta_table_3$Rel
  row.names(fused_table) <- meta_table_3$ID
  fused_table$V2 <- new_table_3$Rel[match(row.names(fused_table),new_table_3$ID)]
  colnames(fused_table) <- c("metagenome_abundances", "metric_abundances")
  
  fused_table_2 <- na.omit(fused_table[fused_table$metric_abundances != 0,])
  fused_table_3 <- na.omit(fused_table_2[fused_table_2$metagenome_abundances != 0,])
  fused_table_3$Strain <- meta_table_3$Cluster[match(row.names(fused_table_3), meta_table_3$ID)]
  
  fused_table_3$NRMSE <- NA
  
  for (i in 1:length(fused_table_3$metagenome_abundances)){
    sub_table <- fused_table_3[i,]
    t <- sub_table$metagenome_abundances
    w <- sub_table$metric_abundances
    
    value <- (abs(w - t)^2)/((w+t)/2)
    fused_table_3$NRMSE[i] <- value
  }
  
  value_sub <- sqrt(1/length(fused_table_3$NRMSE) * sum(fused_table_3$NRMSE))
  
  correlations$NRMSE[correlations$vector == paste0(metric)] <- value_sub  
}

correlations$Amplicon <- c("V3-V4", "V3-V4", "V3-V4", "V3-V4", "V5-V7", "V5-V7", "V5-V7", "V5-V7")
correlations$vector_2<-  gsub("_....*","",correlations$vector_2)

plot_line <- ggplot(correlations, aes(x = vector_2, y = NRMSE, group = Amplicon)) + 
  geom_line(aes(color = Amplicon),size=1.2) +
  geom_point(aes(color = Amplicon)) + 
  ggtitle("SyFi with different pseudoalignment parameters") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("NRMSE") + 
  xlab("Percentage identity pseudoalignment") +
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24)) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
plot_line

pdf(paste(results.dir,"Figure_S6.pdf", sep=""), width=8, height=6)
print(plot_line)
dev.off()

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

###Table S5 - presence/absence - Genus =====
taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

vector <- c("95_V3V4", "97_V3V4", "98_V3V4", "99_V3V4", "100_V3V4", 
            "Salmon_V3V4", "SyFi_V3V4","Mothur_V3V4",
            "Kraken2:Amplicon_V3V4","Kraken3:Amplicon_V3V4","Kraken5:Amplicon_V3V4",
            "Kraken2:Fingerprint_V3V4","Kraken3:Fingerprint_V3V4","Kraken5:Fingerprint_V3V4",
            "95_V5V7", "97_V5V7", "98_V5V7", "99_V5V7", "100_V5V7", 
            "Salmon_V5V7","SyFi_V5V7", "Mothur_V5V7",
            "Kraken2:Amplicon_V5V7","Kraken3:Amplicon_V5V7","Kraken5:Amplicon_V5V7",
            "Kraken2:Fingerprint_V5V7","Kraken3:Fingerprint_V5V7","Kraken5:Fingerprint_V5V7")

#Function Relative abundance
Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

fused_table_4_occ <- data.frame()

for (metric in vector) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory, "comparison_datasets_presence_absence/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert table for comparison
  otu_table = read.table(paste0(working_directory,"comparison_datasets_presence_absence/AtSC_", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep =""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  #Reformatting table of interest 
  no_of_genera <- length(na.omit(unique(taxonomy$genus)))
  
  taxonomy <- read.table(paste(working_directory, "AtSC_taxonomy_GTDB.tsv",sep= ""), sep = "\t", row.names =1, header =T)
  
  #Create new table with the right dimensions
  new_table <- data.frame(matrix(NA, nrow = no_of_genera, ncol = length(row.names(samples_meta))+1))
  
  genera_all <- na.omit(unique(taxonomy$genus))
  
  for (number in 1:no_of_genera){
    number_2 <- genera_all[number]
    new_table[number,1] <- paste(number_2)
  }
  
  column_names <- c("Genus", paste(row.names(samples_meta),sep=""))
  colnames(new_table) <- column_names
  
  row.names(new_table) <- new_table$Genus
  new_table <- new_table %>% dplyr::select(-Genus)
  
  #Rename and reorder column names of otu_table
  colnames(otu_table) = gsub("_t3", "", colnames(otu_table))
  colnames(otu_table) = gsub("_t0", "", colnames(otu_table))
  colnames(otu_table) = gsub("_undil", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v5_v7", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v3_v4", "", colnames(otu_table))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_table_2 <- otu_table[column_names_order]
  otu_table <- otu_table_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:length(genera_all)){
    genus_pick <- genera_all[row]
    
    isolates <- na.omit(row.names(taxonomy)[taxonomy$genus == paste(genus_pick)])
    
    if (length(isolates) > 1){
      other_list <- row.names(otu_table) 
      the_isolate <- isolates[isolates %in% other_list]
    } else {
      the_isolate <- isolates
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_table[row.names(otu_table) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    new_table[row,] <- otu_table_3
  }
  
  #Create new table with the right dimensions for metagenome table
  meta_table <- data.frame(matrix(NA, nrow = no_of_genera, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_genera){
    number_2 <- genera_all[number]
    meta_table[number,1] <- paste(number_2)
  }
  
  column_names <- c("Genus", paste(row.names(samples_meta),sep=""))
  colnames(meta_table) <- column_names
  
  row.names(meta_table) <- meta_table$Genus
  meta_table <- meta_table %>% dplyr::select(-Genus)
  
  #Reshuffle columns of the metagenome data
  colnames(otu_meta) = gsub("_undil", "", colnames(otu_meta))
  #colnames(otu_meta) = gsub("_HL_", "_HL_orig_", colnames(otu_meta))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_meta_2 <- otu_meta[colnames(otu_meta) %in% column_names_order]
  otu_meta <- otu_meta_2[column_names_order]
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:length(genera_all)){
    genus_pick <- genera_all[row]
    isolates <- na.omit(row.names(taxonomy)[taxonomy$genus == paste(genus_pick)])
    
    if (length(isolates) > 1){
      other_list <- row.names(otu_meta) 
      the_isolate <- isolates[isolates %in% other_list]
    } else {
      the_isolate <- isolates
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_meta[row.names(otu_meta) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    meta_table[row,] <- otu_table_3
  }
  
  #Melt tables
  meta_table$Genus <- row.names(meta_table)
  meta_table_2 <- melt(meta_table)
  new_table$Genus <- row.names(new_table)
  new_table_2 <- melt(new_table)
  
  #Make relative abundances
  meta_table_3 <- Hank_the_normalizer(meta_table_2,"variable","value")
  new_table_3 <- Hank_the_normalizer(new_table_2,"variable","value")
  
  #Make relative abundances
  meta_table_3$Rel[meta_table_3$Rel < 0.005] <- 0
  meta_table_3$Rel[meta_table_3$Rel >= 0.005] <- 1 
  new_table_3$Rel[new_table_3$Rel < 0.005] <- 0
  new_table_3$Rel[new_table_3$Rel >= 0.005] <- 1 
  
  meta_table_3$ID <- paste(meta_table_3$variable,meta_table_3$Genus,sep="_")
  new_table_3$ID <- paste(new_table_3$variable,new_table_3$Genus,sep="_")
  
  fused_table <- data.frame(matrix(NA, nrow = length(meta_table_3$ID), ncol = 0))
  fused_table$V1 <- meta_table_3$Rel
  row.names(fused_table) <- meta_table_3$ID
  fused_table$V2 <- new_table_3$Rel[match(row.names(fused_table),new_table_3$ID)]
  colnames(fused_table) <- c("metagenome_abundances", "metric_abundances")
  
  fused_table$metric <- paste(metric) 
  fused_table_4_occ <- rbind(fused_table_4_occ,fused_table)
}

#Presence/Absence
fused_table_4_occ$all <- paste(fused_table_4_occ$metric, fused_table_4_occ$metagenome_abundances, fused_table_4_occ$metric_abundances, sep = "_")

all <- data.frame(table(fused_table_4_occ$all))

all_data <- data.frame()

for (group in unique(fused_table_4_occ$metric)){
  all_new <- data.frame(t(data.frame(c(paste(group),"Correctly present", paste(group,1,1, sep = "_")), c(paste(group), "Incorrectly present", paste(group,0,1, sep = "_")), c(paste(group), "Correctly absent", paste(group,0,0, sep = "_")), c(paste(group), "Incorrectly absent", paste(group,1,0, sep = "_")))))
  all_data <- rbind(all_data, all_new)
}

row.names(all_data) <- NULL
colnames(all_data) <- c("Metric", "Group", "all")

all_data$No_of_SynCom_isolates <- all$Freq[match(all_data$all, all$Var1)]

all_data <- Hank_the_normalizer(all_data, "Metric", "No_of_SynCom_isolates")
all_data$Rel <- as.numeric(all_data$Rel)
all_data$Rel_2 <- paste(round(all_data$Rel*100,1), "%", sep = "")

all_data

write.table(all_data, paste(results.dir, "Table_S5_genus.tsv", sep = ""),sep = "\t", quote =F, col.names =T, row.names =F)

###Table S5 - presence/absence - Species =====
taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

vector <- c("95_V3V4", "97_V3V4", "98_V3V4", "99_V3V4", "100_V3V4", 
            "Salmon_V3V4", "SyFi_V3V4","Mothur_V3V4",
            "Kraken2:Amplicon_V3V4","Kraken3:Amplicon_V3V4","Kraken5:Amplicon_V3V4",
            "Kraken2:Fingerprint_V3V4","Kraken3:Fingerprint_V3V4","Kraken5:Fingerprint_V3V4",
            "95_V5V7", "97_V5V7", "98_V5V7", "99_V5V7", "100_V5V7", 
            "Salmon_V5V7","SyFi_V5V7", "Mothur_V5V7",
            "Kraken2:Amplicon_V5V7","Kraken3:Amplicon_V5V7","Kraken5:Amplicon_V5V7",
            "Kraken2:Fingerprint_V5V7","Kraken3:Fingerprint_V5V7","Kraken5:Fingerprint_V5V7")

#Function Relative abundance
Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

fused_table_4_occ <- data.frame()

for (metric in vector) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory, "comparison_datasets_presence_absence/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert table for comparison
  otu_table = read.table(paste0(working_directory,"comparison_datasets_presence_absence/AtSC_", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep =""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  taxonomy <- read.table(paste(working_directory, "species_clustering.txt",sep= ""), sep = "\t", row.names =1, header =T)
  
  #Reformatting table of interest 
  no_of_species <- length(na.omit(unique(taxonomy$primary_cluster)))
  
  #Create new table with the right dimensions
  new_table <- data.frame(matrix(NA, nrow = no_of_genera, ncol = length(row.names(samples_meta))+1))
  
  species_all <- na.omit(unique(taxonomy$primary_cluster))
  
  for (number in 1:no_of_species){
    number_2 <- species_all[number]
    new_table[number,1] <- paste(number_2)
  }
  
  column_names <- c("Species", paste(row.names(samples_meta),sep=""))
  colnames(new_table) <- column_names
  
  row.names(new_table) <- new_table$Species
  new_table <- new_table %>% dplyr::select(-Species)
  
  #Rename and reorder column names of otu_table
  colnames(otu_table) = gsub("_t3", "", colnames(otu_table))
  colnames(otu_table) = gsub("_t0", "", colnames(otu_table))
  colnames(otu_table) = gsub("_undil", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v5_v7", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v3_v4", "", colnames(otu_table))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_table_2 <- otu_table[column_names_order]
  otu_table <- otu_table_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:length(species_all)){
    species_pick <- species_all[row]
    
    isolates <- na.omit(row.names(taxonomy)[taxonomy$primary_cluster == paste(species_pick)])
    
    if (length(isolates) > 1){
      other_list <- row.names(otu_table) 
      the_isolate <- isolates[isolates %in% other_list]
    } else {
      the_isolate <- isolates
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_table[row.names(otu_table) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    new_table[row,] <- otu_table_3
  }
  
  #Create new table with the right dimensions for metagenome table
  meta_table <- data.frame(matrix(NA, nrow = no_of_species, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_species){
    number_2 <- species_all[number]
    meta_table[number,1] <- paste(number_2)
  }
  
  column_names <- c("Genus", paste(row.names(samples_meta),sep=""))
  colnames(meta_table) <- column_names
  
  row.names(meta_table) <- meta_table$Genus
  meta_table <- meta_table %>% dplyr::select(-Genus)
  
  #Reshuffle columns of the metagenome data
  colnames(otu_meta) = gsub("_undil", "", colnames(otu_meta))
  #colnames(otu_meta) = gsub("_HL_", "_HL_orig_", colnames(otu_meta))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_meta_2 <- otu_meta[colnames(otu_meta) %in% column_names_order]
  otu_meta <- otu_meta_2[column_names_order]
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:length(species_all)){
    species_pick <- species_all[row]
    isolates <- na.omit(row.names(taxonomy)[taxonomy$primary_cluster == paste(species_pick)])
    
    if (length(isolates) > 1){
      other_list <- row.names(otu_meta) 
      the_isolate <- isolates[isolates %in% other_list]
    } else {
      the_isolate <- isolates
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_meta[row.names(otu_meta) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    meta_table[row,] <- otu_table_3
  }
  
  #Melt tables
  meta_table$Species <- row.names(meta_table)
  meta_table_2 <- melt(meta_table)
  new_table$Species <- row.names(new_table)
  new_table_2 <- melt(new_table)
  
  #Make relative abundances
  meta_table_3 <- Hank_the_normalizer(meta_table_2,"variable","value")
  new_table_3 <- Hank_the_normalizer(new_table_2,"variable","value")
  
  #Make relative abundances
  meta_table_3$Rel[meta_table_3$Rel < 0.005] <- 0
  meta_table_3$Rel[meta_table_3$Rel >= 0.005] <- 1 
  new_table_3$Rel[new_table_3$Rel < 0.005] <- 0
  new_table_3$Rel[new_table_3$Rel >= 0.005] <- 1 
  
  meta_table_3$ID <- paste(meta_table_3$variable,meta_table_3$Species,sep="_")
  new_table_3$ID <- paste(new_table_3$variable,new_table_3$Species,sep="_")
  
  fused_table <- data.frame(matrix(NA, nrow = length(meta_table_3$ID), ncol = 0))
  fused_table$V1 <- meta_table_3$Rel
  row.names(fused_table) <- meta_table_3$ID
  fused_table$V2 <- new_table_3$Rel[match(row.names(fused_table),new_table_3$ID)]
  colnames(fused_table) <- c("metagenome_abundances", "metric_abundances")
  
  fused_table$metric <- paste(metric) 
  fused_table_4_occ <- rbind(fused_table_4_occ,fused_table)
}

#Presence/Absence
fused_table_4_occ$all <- paste(fused_table_4_occ$metric, fused_table_4_occ$metagenome_abundances, fused_table_4_occ$metric_abundances, sep = "_")

all <- data.frame(table(fused_table_4_occ$all))

all_data <- data.frame()

for (group in unique(fused_table_4_occ$metric)){
  all_new <- data.frame(t(data.frame(c(paste(group),"Correctly present", paste(group,1,1, sep = "_")), c(paste(group), "Incorrectly present", paste(group,0,1, sep = "_")), c(paste(group), "Correctly absent", paste(group,0,0, sep = "_")), c(paste(group), "Incorrectly absent", paste(group,1,0, sep = "_")))))
  all_data <- rbind(all_data, all_new)
}

row.names(all_data) <- NULL
colnames(all_data) <- c("Metric", "Group", "all")

all_data$No_of_SynCom_isolates <- all$Freq[match(all_data$all, all$Var1)]

all_data <- Hank_the_normalizer(all_data, "Metric", "No_of_SynCom_isolates")
all_data$Rel <- as.numeric(all_data$Rel)
all_data$Rel_2 <- paste(round(all_data$Rel*100,1), "%", sep = "")

all_data

write.table(all_data, paste(results.dir, "Table_S5_species.tsv", sep = ""),sep = "\t", quote =F, col.names =T, row.names =F)


###Table S5 - presence/absence - Strain =====
taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

vector <- c("95_V3V4", "97_V3V4", "98_V3V4", "99_V3V4", "100_V3V4", 
            "Salmon_V3V4", "SyFi_V3V4","Mothur_V3V4",
            "Kraken2:Amplicon_V3V4","Kraken3:Amplicon_V3V4","Kraken5:Amplicon_V3V4",
            "Kraken2:Fingerprint_V3V4","Kraken3:Fingerprint_V3V4","Kraken5:Fingerprint_V3V4",
            "95_V5V7", "97_V5V7", "98_V5V7", "99_V5V7", "100_V5V7", 
            "Salmon_V5V7","SyFi_V5V7", "Mothur_V5V7",
            "Kraken2:Amplicon_V5V7","Kraken3:Amplicon_V5V7","Kraken5:Amplicon_V5V7",
            "Kraken2:Fingerprint_V5V7","Kraken3:Fingerprint_V5V7","Kraken5:Fingerprint_V5V7")

#Function Relative abundance
Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

fused_table_4_occ <- data.frame()

for (metric in vector) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory, "comparison_datasets_presence_absence/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert table for comparison
  otu_table = read.table(paste0(working_directory,"comparison_datasets_presence_absence/AtSC_", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep =""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  #Insert clusters file (which isolate is in which cluster based on amplicon sequence or fingerprint and at different sequence identities)
  if (metric %in% c("Salmon_V3V4", "Kraken2:Amplicon_V3V4","Kraken3:Amplicon_V3V4","Kraken5:Amplicon_V3V4")){
    metric_2 <- "100_V3V4"
  } else if (metric %in% c("Salmon_V5V7", "Kraken2:Amplicon_V5V7","Kraken3:Amplicon_V5V7","Kraken5:Amplicon_V5V7")){
    metric_2 <- "100_V5V7"
  } else if (metric %in% c("Kraken2:Fingerprint_V3V4","Kraken3:Fingerprint_V3V4","Kraken5:Fingerprint_V3V4")) {
    metric_2 <- "SyFi_V3V4"
  } else if (metric %in% c("Kraken2:Fingerprint_V5V7","Kraken3:Fingerprint_V5V7","Kraken5:Fingerprint_V5V7")) {
    metric_2 <- "SyFi_V5V7"
  } else {
    metric_2 <- metric
  }
  
  clusters = read.csv(paste0(working_directory,"vsearch_clusters_presence_absence/collapse_",metric_2,".txt",sep=""), header=F, sep="\t")
  
  #Reformatting table of interest 
  no_of_isolates <- length(row.names(clusters))
  
  #Create new table with the right dimensions
  new_table <- data.frame(matrix(NA, nrow = no_of_isolates, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_isolates){
    new_table[number,1] <- paste("Cluster_", number,sep="")
  }
  
  column_names <- c("cluster", paste(row.names(samples_meta),sep=""))
  colnames(new_table) <- column_names
  
  row.names(new_table) <- new_table$cluster
  new_table <- new_table %>% dplyr::select(-cluster)
  
  #Reorder column names of otu_table
  colnames(otu_table) = gsub("_t3", "", colnames(otu_table))
  colnames(otu_table) = gsub("_t0", "", colnames(otu_table))
  colnames(otu_table) = gsub("_undil", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v5_v7", "", colnames(otu_table))
  colnames(otu_table) = gsub("_v3_v4", "", colnames(otu_table))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_table_2 <- otu_table[column_names_order]
  otu_table <- otu_table_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:nrow(clusters)){
    isolate_vector <- as.character(clusters[row,])
    isolate_vector <- isolate_vector[!is.na(isolate_vector)]
    isolate_vector <- isolate_vector[!isolate_vector == ""]
    
    if (length(isolate_vector) > 1){
      other_list <- row.names(otu_table) 
      the_isolate <- isolate_vector[isolate_vector %in% other_list]
    } else {
      the_isolate <- isolate_vector
    }
    
    #subset otu_table of choice
    otu_table_2 <- otu_table[row.names(otu_table) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    new_table[row,] <- otu_table_3
  }
  
  #Create new table with the right dimensions for metagenome table
  meta_table <- data.frame(matrix(NA, nrow = no_of_isolates, ncol = length(row.names(samples_meta))+1))
  
  for (number in 1:no_of_isolates){
    meta_table[number,1] <- paste("Cluster_", number,sep="")
  }
  
  column_names <- c("cluster", paste(row.names(samples_meta),sep=""))
  colnames(meta_table) <- column_names
  
  row.names(meta_table) <- meta_table$cluster
  meta_table <- meta_table %>% dplyr::select(-cluster)
  
  colnames(otu_meta) = gsub("_undil", "", colnames(otu_meta))
  column_names_order <- c(paste(row.names(samples_meta),sep=""))
  otu_meta_2 <- otu_meta[column_names_order]
  otu_meta <- otu_meta_2
  
  #Merge values from isolates of the same cluster - table for comparison
  for (row in 1:nrow(clusters)){
    isolate_vector <- as.character(clusters[row,])
    other_list <- row.names(otu_meta) 
    the_isolate <- isolate_vector[isolate_vector %in% other_list]
    the_isolate <- the_isolate[!is.na(the_isolate)]
    the_isolate <- the_isolate[!the_isolate == ""]
    
    #subset otu_table of choice
    otu_table_2 <- otu_meta[row.names(otu_meta) %in% paste0(the_isolate),]
    otu_table_3 <- colSums(otu_table_2)
    
    meta_table[row,] <- otu_table_3
  }
  
  #Melt tables
  meta_table$Cluster <- row.names(meta_table)
  meta_table_2 <- melt(meta_table)
  new_table$Cluster <- row.names(new_table)
  new_table_2 <- melt(new_table)
  
  #Make relative abundances
  meta_table_3 <- Hank_the_normalizer(meta_table_2,"variable","value")
  new_table_3 <- Hank_the_normalizer(new_table_2,"variable","value")
  
  #Make relative abundances
  meta_table_3$Rel[meta_table_3$Rel < 0.005] <- 0
  meta_table_3$Rel[meta_table_3$Rel >= 0.005] <- 1 
  new_table_3$Rel[new_table_3$Rel < 0.005] <- 0
  new_table_3$Rel[new_table_3$Rel >= 0.005] <- 1 
  
  meta_table_3$ID <- paste(meta_table_3$variable,meta_table_3$Cluster,sep="_")
  new_table_3$ID <- paste(new_table_3$variable,new_table_3$Cluster,sep="_")
  
  fused_table <- data.frame(matrix(NA, nrow = length(meta_table_3$ID), ncol = 0))
  fused_table$V1 <- meta_table_3$Rel
  row.names(fused_table) <- meta_table_3$ID
  fused_table$V2 <- new_table_3$Rel[match(row.names(fused_table),new_table_3$ID)]
  colnames(fused_table) <- c("metagenome_abundances", "metric_abundances")
  
  fused_table$metric <- paste(metric) 
  fused_table_4_occ <- rbind(fused_table_4_occ,fused_table)
}

#Presence/Absence
fused_table_4_occ$all <- paste(fused_table_4_occ$metric, fused_table_4_occ$metagenome_abundances, fused_table_4_occ$metric_abundances, sep = "_")

all <- data.frame(table(fused_table_4_occ$all))

all_data <- data.frame()

for (group in unique(fused_table_4_occ$metric)){
  all_new <- data.frame(t(data.frame(c(paste(group),"Correctly present", paste(group,1,1, sep = "_")), c(paste(group), "Incorrectly present", paste(group,0,1, sep = "_")), c(paste(group), "Correctly absent", paste(group,0,0, sep = "_")), c(paste(group), "Incorrectly absent", paste(group,1,0, sep = "_")))))
  all_data <- rbind(all_data, all_new)
}

row.names(all_data) <- NULL
colnames(all_data) <- c("Metric", "Group", "all")

all_data$No_of_SynCom_isolates <- all$Freq[match(all_data$all, all$Var1)]

all_data <- Hank_the_normalizer(all_data, "Metric", "No_of_SynCom_isolates")
all_data$Rel <- as.numeric(all_data$Rel)
all_data$Rel_2 <- paste(round(all_data$Rel*100,1), "%", sep = "")

all_data

write.table(all_data, paste(results.dir, "Table_S5_strain.tsv", sep = ""),sep = "\t", quote =F, col.names =T, row.names =F)


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
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
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
  scale_color_manual(values =c("#D81B60","#1E88E5") ) +
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
