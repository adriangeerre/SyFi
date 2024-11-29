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

###Figure S7 - SyFi accuracy - full =====

taxonomy = read.table(paste(working_directory,"AtSC_taxonomy_GTDB.tsv", sep =""), header=T,sep="\t",quote="\"", fill = FALSE)

vector <- c("95_V3V4", "97_V3V4", "98_V3V4", "99_V3V4", "100_V3V4", "Salmon_V3V4","SyFi_V3V4","95_V5V7", "97_V5V7", "98_V5V7", "99_V5V7", "100_V5V7", "Salmon_V5V7","SyFi_V5V7")

#Function Relative abundance
Hank_the_normalizer <- function(df,group,amount){
  df_2 <- df %>% dplyr::group_by_at(group) %>% dplyr::summarise(total=sum(.data[[amount]]))
  df_3 <- df_2$total
  names(df_3) <- df_2[[group]]
  df$total <- df_3[as.character(df[[group]])]
  df$Rel <- df[[amount]] / df$total
  return(df)
}

fused_table_4 <- data.frame()

correlations <- as.data.frame(vector)
correlations$Correlation <- NA

for (metric in vector) {
  #Insert metagenome table
  otu_meta = read.table(paste(working_directory, "comparison_datasets/AtSC_meta_norm.tsv", sep =""), header=T, sep="\t")
  
  row.names(otu_meta) <- otu_meta$isolate
  otu_meta <- otu_meta %>% dplyr::select(-isolate)
  
  #Insert table for comparison
  otu_table = read.table(paste0(working_directory, "comparison_datasets/AtSC_", metric, ".tsv",sep =""), header=T, sep="\t")
  
  row.names(otu_table) <- otu_table$isolate
  otu_table <- otu_table %>% dplyr::select(-isolate)
  
  #Insert metadata file
  samples_meta = read.table(paste(working_directory,"metadata.txt", sep =""), header=TRUE,sep="\t")
  
  row.names(samples_meta) <- samples_meta$sample_id
  samples_meta <- samples_meta %>% dplyr::select(-sample_id)
  
  #Insert clusters file (which isolate is in which cluster based on amplicon sequence or fingerprint and at different sequence identities)
  if (metric == "Salmon_V3V4"){
    metric_2 <- "100_V3V4"
  } else if (metric == "Salmon_V5V7"){
    metric_2 <- "100_V5V7"
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
  
  fused_table_3$metric <- paste(metric) 
  fused_table_4 <- rbind(fused_table_4,fused_table_3)
  
  value <- cor(fused_table_3$metagenome_abundances,fused_table_3$metric_abundances)
  
  correlations$Correlation[correlations$vector == paste0(metric)] <- value
}


#Vsearch clustering tables - Abundances

fused_table_V3V4 <- fused_table_4[grepl("V3V4",fused_table_4$metric),]
fused_table_V3V4$test <- gsub("_....", "", fused_table_V3V4$metric)
fused_table_V3V4$test_2 <- "%"
fused_table_V3V4$Metric <- paste(fused_table_V3V4$test,fused_table_V3V4$test_2,sep="")
fused_table_V3V4$Metric[fused_table_V3V4$Metric == "SyFi%"] <- "SyFi"
fused_table_V3V4$Metric[fused_table_V3V4$Metric == "Salmon%"] <- "Salmon"

fused_table_V3V4$Metric <- factor(fused_table_V3V4$Metric, levels = c("95%", "97%", "98%", "99%", "100%", "Salmon", "SyFi"))

V3V4 <- ggscatter(fused_table_V3V4, x="metagenome_abundances", y="metric_abundances", color = "Metric",conf.int = TRUE,alpha = 0.2,
                  palette = "jco", add = "reg.line", xlim = c(0,0.8), ylim = c(0,0.8)) +
  stat_cor(aes(color = Metric), label.x = 0.55) +
  ggtitle("V3-V4") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("Metric - relative abundances") + 
  xlab("Metagenome - relative abundances") + 
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
V3V4

fused_table_V5V7 <- fused_table_4[grepl("V5V7",fused_table_4$metric),]
fused_table_V5V7$test <- gsub("_....", "", fused_table_V5V7$metric)
fused_table_V5V7$test_2 <- "%"
fused_table_V5V7$Metric <- paste(fused_table_V5V7$test,fused_table_V5V7$test_2,sep="")
fused_table_V5V7$Metric[fused_table_V5V7$Metric == "SyFi%"] <- "SyFi"
fused_table_V5V7$Metric[fused_table_V5V7$Metric == "Salmon%"] <- "Salmon"
fused_table_V5V7$Metric <- factor(fused_table_V5V7$Metric, levels = c("95%", "97%", "98%", "99%", "100%", "Salmon", "SyFi"))

V5V7 <- ggscatter(fused_table_V5V7, x="metagenome_abundances", y="metric_abundances", color = "Metric",conf.int = TRUE,alpha = 0.2,
                  palette = "jco", add = "reg.line", xlim = c(0,0.8), ylim = c(0,0.8)) +
  stat_cor(aes(color = Metric), label.x = 0.55) +
  ggtitle("V5-V7") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ylab("Metric - relative abundances") + 
  xlab("Metagenome - relative abundances") + 
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18), axis.text.y = element_text(size=14), legend.title = element_text(size=18), legend.text = element_text(size=14), plot.title = element_text(size=24))
V5V7

figure <- ggarrange(V3V4, V5V7,
                    labels = c("A", "B" ),font.label = list(size = 28),ncol = 2, nrow = 1, widths = 15, heights = 8, common.legend = T)

pdf(paste(results.dir,"Figure_S7.pdf", sep=""), width=15, height=8)
print(figure)
dev.off()


