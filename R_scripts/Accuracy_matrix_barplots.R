library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)

file_list <- list.files(c("F:/Cohort/github/dataset/"),pattern = "differential_abundance_accuracy_matrix_*",full.names = TRUE)
print(file_list)
accuracy <- do.call("rbind", lapply(file_list, FUN = function(file) {
  read.table(file, header=TRUE, sep="\t")
}))

# Create a list to store shuffled matrices
accuracy_1<-reshape2::melt(accuracy)
colnames(accuracy_1)<-c("Dataset","Cohort","Methods","Pattern","Metrices","Value")
accuracy_1$Method_1<-paste0(accuracy_1$Methods,"_",accuracy_1$Pattern)
accuracy_1<-dplyr::filter(accuracy_1,Metrices!="accuracy")
accuracy_1$Metrices<-str_to_title(accuracy_1$Metrices)
levels(accuracy_1$Metrices)<-c("F1","Recall","Precision")
accuracy_1$Pattern<-as.factor(accuracy_1$Pattern)
levels(accuracy_1$Pattern)<-c("DEFAULT","CUSTOM")

colnames(accuracy_1)

pdf("F:/Cohort/github/Figures/accuracy_differential_barplt.pdf",width=120, height=50)
ggplot(accuracy_1,aes(x=Method_1, y=Value, fill=Methods,pattern=Pattern)) +
  geom_bar(stat="identity") +
  geom_bar_pattern(stat = "identity",
                   position = "stack",
                   pattern_color = "white",
                   color = "black",
                   aes(pattern_color = white)) +
  scale_pattern_manual(values = c(DEFAULT="none",CUSTOM = "stripe"))+
  facet_grid(Metrices~ Cohort, scales = "free", space = "free", switch="y") +
  theme(panel.grid.major.x = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=100,face="bold", hjust=1,angle=90),
        axis.text.y= element_text(face="bold",size = 100),
        strip.text.x = element_text(face="bold",size = 100),
        legend.position = "none",
        panel.spacing.x=unit(0.3, "lines") , panel.spacing.y=unit(7,"lines"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.text = element_text(size=100, face="bold"),
        legend.text=element_text(size=100,face="bold"),
        legend.key.size = unit(1,"line"), 
        legend.title=element_text(size=100, face="bold"),strip.background =element_rect(fill="#FFEAD2")
  )+ scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  scale_fill_manual(values=c("#554994","#FF8DC7","#F29393","#FFCCB3"))+
  scale_x_discrete(labels=c("METGEM_DEFAULT"= "Metgem","METGEM_CUSTOM" = "Metgem", 
                            "PANFP_DEFAULT"= "PanFP","PANFP_CUSTOM" = "PanFP",
                            "PICRUST2_DEFAULT"="PICRUSt2","PICRUST2_CUSTOM"="PICRUSt2",
                            "TAX4FUN2_DEFAULT"="Tax4Fun2","TAX4FUN2_CUSTOM"="Tax4Fun2"))
dev.off()
