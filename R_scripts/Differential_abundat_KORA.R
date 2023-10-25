library("ggpubr")
library(dplyr)
#select the column number which has the nature of the samples
coln<-4
setwd("F:/Cohort/github/dataset/")
map_data<-read.table("KORA/KORA_METADATA.tsv",sep="\t",header=TRUE,check.names = FALSE)
class(map_data$Sample_ID)<-"character"
map_data<-dplyr::filter(map_data,DIABETES=="Healthy" | DIABETES=="Diabetes")
map_data<-dplyr::filter(map_data,DIABETES=="Healthy" | DIABETES=="Diabetes")
KO_WGS<-read.table("KORA/KORA_MGS_KO.tsv",sep="\t",header = TRUE,check.names = FALSE,row.names=1)
KO_WGS<-KO_WGS %>% dplyr::select(one_of(dput(as.character(map_data$Sample_ID))))
KO_WGS<-OTUtable::filter_taxa(KO_WGS, abundance=0,persistence=5)
#KO_WGS<-microbiome::transform(KO_WGS, 'compositional')
KO_WGS<-as.data.frame(KO_WGS)
  #common<-intersect(rownames(KO_table),rownames(KO_WGS))
  #KO_table<-KO_table[rownames(KO_table) %in% common, ]
  #KO_WGS<-KO_WGS[rownames(KO_WGS) %in% common, ]
KO_WGS<-as.data.frame(t(KO_WGS))
KO_WGS$Sample_ID<-rownames(KO_WGS)
KO_WGS<-dplyr::left_join(KO_WGS,map_data[,c(1,coln)])
rownames(KO_WGS)<-KO_WGS$Sample_ID
KO_WGS$Sample_ID<-NULL
KO_WGS<-reshape2::melt(KO_WGS)
colnames(KO_WGS)<-c("CASE","KO","Abundance")
  
Wilcox_wgs<-compare_means(Abundance ~ CASE,  data=KO_WGS, group.by = 'KO',
                            paired = FALSE, p.adjust.method = "BH",method = "wilcox.test") 

# List files with the specified prefix and suffix
file_list <- list.files(c("KORA/"),pattern = "_REL_KO.tsv",full.names = TRUE)
print(file_list)

# Create a list to store shuffled matrices
combined_matrices <- vector("list", length(file_list))
names(combined_matrices) <- basename(file_list)

# Generate and store shuffled matrices
for (i in 1:length(file_list)){
  # Print the current file path
  print(file_list[[i]])
  KO_table<-read.table(file_list[[i]],sep="\t",header=TRUE,check.names = FALSE,row.names=1)
  KO_table<-KO_table%>% dplyr::select(one_of(dput(as.character(map_data$Sample_ID))))
  #(at any abundance) in at least 5% of samples:
  KO_table<-OTUtable::filter_taxa(KO_table, abundance=0, persistence=5)
  combined_matrices[[i]] <- KO_table
  #names(combined_matrices[[i]]) <- basename(file_list[[i]])
}
result_list <- list()


# Define a function to process each element of shuffled_matrices
process_element <- function(j) {
  KO_table <- as.data.frame(t(combined_matrices[[j]]))
  KO_table$Sample_ID <- rownames(KO_table)
  KO_table <- dplyr::left_join(KO_table, map_data[, c(1, coln)])
  rownames(KO_table) <- KO_table$Sample_ID
  KO_table$Sample_ID <- NULL
  KO_table <- reshape2::melt(KO_table)
  colnames(KO_table) <- c("CASE", "KO", "Abundance")
  Wilcox_p1 <- ggpubr::compare_means(Abundance ~ CASE, data = KO_table, group.by = 'KO',
                                     paired = FALSE, p.adjust.method = "BH", method = "wilcox.test")
  Wilcox_p1_sig <- dplyr::filter(Wilcox_p1, p < 0.05)
  colnames(Wilcox_p1) <- paste(colnames(Wilcox_p1), "p1", sep = "_")
  colnames(Wilcox_p1)[1] <- "KO"
  colnames(Wilcox_wgs)[1] <- "KO"
  Wilcox_p1_1 <- Wilcox_p1[, c(1, 6, 7)]
  Wilcox_wgs_1 <- Wilcox_wgs[, c(1, 6, 7)]
  KORA_p1 <- merge(Wilcox_p1_1, Wilcox_wgs_1, all = TRUE)
  KORA_p1_adj <- KORA_p1[, c(1, 2, 4)]
  colnames(KORA_p1_adj) <- c("KO", "prediction", "metagenom")
  KORA_p1_adj <- KORA_p1_adj[, c(1, 3, 2)]
  KORA_p1_p <- KORA_p1[, c(1, 3, 5)]
  colnames(KORA_p1_p) <- c("KO", "prediction", "metagenom")
  KORA_p1_p <- KORA_p1_p[, c(1, 3, 2)]
  common_genes <- intersect(Wilcox_wgs$KO, Wilcox_p1$KO)
  Wilcox_wgs_1 <- subset(Wilcox_wgs, KO %in% common_genes)
  Wilcox_p1_1 <- subset(Wilcox_p1_1, KO %in% common_genes)
  threshold <- 0.05
  tp <- sum(Wilcox_wgs_1$p <= threshold & Wilcox_p1_1$p.format_p1 <= threshold)
  tn <- sum(Wilcox_wgs_1$p > threshold & Wilcox_p1_1$p.format_p1 > threshold)
  fp <- sum(Wilcox_wgs_1$p > threshold & Wilcox_p1_1$p.format_p1 <= threshold)
  fn <- sum(Wilcox_wgs_1$p <= threshold & Wilcox_p1_1$p.format_p1 > threshold)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * precision * recall / (precision + recall)
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  result <- data.frame(
   Data = sub("\\..*", "", names(combined_matrices)[j]),
    precision = precision,
    recall = recall,
    f1 = f1,
    accuracy = accuracy
  )
  return(result)
}

# Apply the function to each element of shuffled_matrices
result_list <- lapply(1:length(combined_matrices), process_element)
# Combine the results into a single data frame
result_df <- do.call(rbind, result_list)
result_df$Cohort<-gsub("(.+?)(\\_.*)", "\\1", result_df$Data)
result_df$Methods<-gsub("^[^_]*_([^_]+)_.*$", "\\1",result_df$Data)
result_df$Pattern <- ifelse(grepl("CUSTOM", result_df$Data), "CUSTOM", "DEFAULT")
write.table(result_df,file="differential_abundance_accuracy_matrix_KORA.tsv",sep="\t",quote=FALSE,row.names = FALSE)
