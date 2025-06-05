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


