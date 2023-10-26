library("ggpubr")
library(dplyr)
coln<-2
setwd("../simulation/")
source("../R_scripts/Differential_abundat.R")
map_data<-read.table("differential_analysis/SIMUATED_METADATA.tsv",sep="\t",header=TRUE,check.names = FALSE)
map_data<-dplyr::filter(map_data,Nature=="GI" | Nature=="oral")
class(map_data$Sample_ID)<-"character"
KO_WGS<-read.table("differential_analysis/MGS_sim.tsv",sep="\t",header = TRUE,check.names = FALSE,row.names=1)
KO_WGS<-KO_WGS %>% dplyr::select(one_of(dput(as.character(map_data$Sample_ID))))
KO_WGS<-OTUtable::filter_taxa(KO_WGS, abundance=0,persistence=5)
KO_WGS<-as.data.frame(t(KO_WGS))
KO_WGS$Sample_ID<-rownames(KO_WGS)
KO_WGS<-dplyr::left_join(KO_WGS,map_data[,c(1,coln)])
rownames(KO_WGS)<-KO_WGS$Sample_ID
KO_WGS$Sample_ID<-NULL
KO_WGS<-reshape2::melt(KO_WGS)
colnames(KO_WGS)<-c("CASE","KO","Abundance")
Wilcox_wgs<-compare_means(Abundance ~ CASE,  data=KO_WGS, group.by = 'KO',
                          paired = FALSE, p.adjust.method = "BH",method = "wilcox.test") 

Wilcox_wgs_sig<-dplyr::filter(Wilcox_wgs,p<0.05)
file_list <- list.files(c("F:/Cohort/meta-apo/dataset1/simulation/differential_analysis/"),pattern = "*_simulation.tsv",full.names = TRUE)
print(file_list)
# Create a list to store shuffled matrices
combined_matrices <- list()
for (file_path in file_list) {
  combined_matrices <- process_file(file_path, map_data, combined_matrices)
}
result_list <- list()
# Define a function to process each element of shuffled_matrices
# Apply the function to each element of shuffled_matrices
result_list <- lapply(1:length(combined_matrices), process_element)
# Combine the results into a single data frame
result_df <- do.call(rbind, result_list)
result_df$Methods<-gsub("(.+?)(\\_.*)", "\\1", result_df$Method)
result_df$Pattern <- ifelse(grepl("CUS|cus", result_df$Method), "CUSTOM", "DEFAULT")
result_df$Nature <- "GI_vs_Oral"
write.table(result_df,file="differential_abundance_accuracy_matrix_simulation_oral_GI.tsv",sep="\t",quote=FALSE,row.names = FALSE)


################Airways vs GI###############################################
map_data<-read.table("differential_analysis/SIMUATED_METADATA.tsv",sep="\t",header=TRUE,check.names = FALSE)
map_data<-dplyr::filter(map_data,Nature=="GI" | Nature=="air")
class(map_data$Sample_ID)<-"character"
KO_WGS<-read.table("differential_analysis/MGS_sim.tsv",sep="\t",header = TRUE,check.names = FALSE,row.names=1)
KO_WGS<-KO_WGS %>% dplyr::select(one_of(dput(as.character(map_data$Sample_ID))))
KO_WGS<-OTUtable::filter_taxa(KO_WGS, abundance=0,persistence=5)
KO_WGS<-as.data.frame(t(KO_WGS))
KO_WGS$Sample_ID<-rownames(KO_WGS)
KO_WGS<-dplyr::left_join(KO_WGS,map_data[,c(1,coln)])
rownames(KO_WGS)<-KO_WGS$Sample_ID
KO_WGS$Sample_ID<-NULL
KO_WGS<-reshape2::melt(KO_WGS)
colnames(KO_WGS)<-c("CASE","KO","Abundance")
Wilcox_wgs<-compare_means(Abundance ~ CASE,  data=KO_WGS, group.by = 'KO',
                          paired = FALSE, p.adjust.method = "BH",method = "wilcox.test") 

Wilcox_wgs_sig<-dplyr::filter(Wilcox_wgs,p<0.05)
file_list <- list.files(c("F:/Cohort/meta-apo/dataset1/simulation/differential_analysis/"),pattern = "*_simulation.tsv",full.names = TRUE)
print(file_list)
# Create a list to store shuffled matrices
combined_matrices <- list()
for (file_path in file_list) {
  combined_matrices <- process_file(file_path, map_data, combined_matrices)
}
result_list <- list()
# Define a function to process each element of shuffled_matrices
# Apply the function to each element of shuffled_matrices
result_list <- lapply(1:length(combined_matrices), process_element)
# Combine the results into a single data frame
result_df <- do.call(rbind, result_list)
result_df$Methods<-gsub("(.+?)(\\_.*)", "\\1", result_df$Method)
result_df$Pattern <- ifelse(grepl("CUS|cus", result_df$Method), "CUSTOM", "DEFAULT")
result_df$Nature <- "GI_vs_Airways"
write.table(result_df,file="differential_abundance_accuracy_matrix_simulation_airways_GI.tsv",sep="\t",quote=FALSE,row.names = FALSE)
