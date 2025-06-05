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
