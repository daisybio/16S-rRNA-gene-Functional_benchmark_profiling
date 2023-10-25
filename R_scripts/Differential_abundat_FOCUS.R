library("ggpubr")
library(dplyr)
#select the column number which has the nature of the samples
coln<-4
source("F:/Cohort/github/Scripts/R_scripts/Differential_abundat_core.R")
setwd("F:/Cohort/github/dataset/")
map_data<-read.table("FOCUS/FOCUS_METADATA.tsv",sep="\t",header=TRUE,check.names = FALSE)
class(map_data$Sample_ID)<-"character"
#map_data<-dplyr::filter(map_data,DIABETES=="Healthy" | DIABETES=="Diabetes")
#map_data<-dplyr::filter(map_data,DIABETES=="Healthy" | DIABETES=="Diabetes")
KO_WGS<-read.table("FOCUS/FOCUS_MGS_KO.tsv",sep="\t",header = TRUE,check.names = FALSE,row.names=1)
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
#FOCUS: 1,4
#FOCUS: 1,5
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
file_list <- list.files(c("FOCUS/"),pattern = "_REL_KO.tsv",full.names = TRUE)
print(file_list)


# Create a list to store shuffled matrices
combined_matrices <- list()
for (file_path in file_list) {
  combined_matrices <- process_file(file_path, map_data, combined_matrices)
}

result_list <- list()
# Apply the function to each element of shuffled_matrices
result_list <- lapply(1:length(combined_matrices), process_element)
# Combine the results into a single data frame
result_df <- do.call(rbind, result_list)
result_df <- do.call(rbind, result_list)
result_df$Cohort<-gsub("(.+?)(\\_.*)", "\\1", result_df$Data)
result_df$Methods<-gsub("^[^_]*_([^_]+)_.*$", "\\1",result_df$Data)
result_df$Pattern <- ifelse(grepl("CUSTOM", result_df$Data), "CUSTOM", "DEFAULT")

write.table(result_df,file="differential_abundance_accuracy_matrix_FOCUS.tsv",sep="\t",quote=FALSE,row.names = FALSE)
