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

