library("ggpubr")
library(dplyr)
#select the column number which has the nature of the samples
coln<-4
source("F:/Cohort/github/Scripts/R_scripts/Differential_abundat_core.R")
setwd("F:/Cohort/github/dataset/")
map_data<-read.table("POPGEN/POPGEN_METADATA.tsv",sep="\t",header=TRUE,check.names = FALSE)
class(map_data$Sample_ID)<-"character"
#map_data<-dplyr::filter(map_data,DIABETES=="Healthy" | DIABETES=="Diabetes")
#map_data<-dplyr::filter(map_data,DIABETES=="Healthy" | DIABETES=="Diabetes")
KO_WGS<-read.table("POPGEN/POPGEN_MGS_KO.tsv",sep="\t",header = TRUE,check.names = FALSE,row.names=1)
KO_WGS<-KO_WGS %>% dplyr::select(one_of(dput(as.character(map_data$Sample_ID))))
KO_WGS<-OTUtable::filter_taxa(KO_WGS, abundance=0,persistence=5)
#KO_WGS<-microbiome::transform(KO_WGS, 'compositional')
KO_WGS<-as.data.frame(KO_WGS)
#common<-intersect(rownames(KO_table),rownames(KO_WGS))
#KO_table<-KO_table[rownames(KO_table) %in% common, ]
#KO_WGS<-KO_WGS[rownames(KO_WGS) %in% common, ]
KO_WGS<-as.data.frame(t(KO_WGS))
KO_WGS$Sample_ID<-rownames(KO_WGS)
#change the column number according to map data
#POPGEN: 1,4
#POPGEN: 1,5
KO_WGS<-dplyr::left_join(KO_WGS,map_data[,c(1,coln)])
rownames(KO_WGS)<-KO_WGS$Sample_ID
KO_WGS$Sample_ID<-NULL
KO_WGS<-reshape2::melt(KO_WGS)
colnames(KO_WGS)<-c("CASE","KO","Abundance")

Wilcox_wgs<-compare_means(Abundance ~ CASE,  data=KO_WGS, group.by = 'KO',
                          paired = FALSE, p.adjust.method = "BH",method = "wilcox.test") 
#Wilcox_wgs_BH<-compare_means(Abundance ~ CASE,  data=KO_WGS, group.by = 'KO',
#write.table(Wilcox_wgs,file="Wilcox_results_mgs_new.tsv",sep="\t",quote=FALSE,row.names = FALSE)
#paired = FALSE, p.adjust.method = "BH",method = "wilcox.test")

Wilcox_wgs_sig<-dplyr::filter(Wilcox_wgs,p<0.05)


# List files with the specified prefix and suffix
file_list <- list.files(c("POPGEN/"),pattern = "_REL_KO.tsv",full.names = TRUE)
print(file_list)


# Create a list to store shuffled matrices
combined_matrices <- list()
for (file_path in file_list) {
  combined_matrices <- process_file(file_path, map_data, combined_matrices)
}

# Generate and store shuffled matrices

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
  common_genes <- intersect(Wilcox_wgs$KO, Wilcox_p1$KO)
  Wilcox_wgs_1 <- subset(Wilcox_wgs, KO %in% common_genes)
  Wilcox_p1_1 <- subset(Wilcox_p1, KO %in% common_genes)
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
    Method = sub("\\..*", "", names(combined_matrices)[j]),
    precision = precision,
    recall = recall,
    f1 = f1,
    accuracy = accuracy
  )
  return(result)
}

# Apply the function to each element of shuffled_matrices
result_list <- lapply(1:length(combined_matrices), process_element)
result_df <- do.call(rbind, result_list)
result_df$Cohort<-gsub("(.+?)(\\_.*)", "\\1", result_df$Data)
result_df$Methods<-gsub("^[^_]*_([^_]+)_.*$", "\\1",result_df$Data)
result_df$Pattern <- ifelse(grepl("CUSTOM", result_df$Data), "CUSTOM", "DEFAULT")
# Combine the results into a single data frame
result_df <- do.call(rbind, result_list)
write.table(result_df,file="differential_abundance_accuracy_matrix_POPGEN.tsv",sep="\t",quote=FALSE,row.names = FALSE)
