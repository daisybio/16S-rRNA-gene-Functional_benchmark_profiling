##The following scripts were adapted and modified to generate distance matrix and PcoA plot from Parallel-Meta Suite: Interactive and rapid microbiome data analysis on multiple platforms  

## install necessary libraries
p <- c("optparse","vegan", "ade4","ggplot2","grid")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('F:/Cohort/meta-apo/dataset1/simulation/')



## load data
# import meta file
meta_orig <- read.table("F:/Cohort/meta-apo/dataset1/simulation//all_meta.txt",header=TRUE,row.names = 1)
#meta_orig <- read.table(file=opts$meta_data, header=TRUE, row.names=1)
# import dist file
dst_orig<-read.table("F:/Cohort/meta-apo/dataset1/simulation//all.dist",header=TRUE,row.names = 1)
#dst_orig <- read.table(file=opts$dist_file, header=TRUE, row.names=1)
da=dst_orig 
md=meta_orig
meta_orig$Method=as.factor(meta_orig$Method)
meta_orig$Platform=as.factor(meta_orig$Platform)
meta_orig$Nature=as.factor(meta_orig$Nature)

## main calc & draw function definition
  rn <- rownames(da)
  cn <- colnames(da)
  
  meta <- meta_orig
  n <- ncol(meta)
  sampleNumber <- length(rn)

  dst <- as.dist(da)
  pcoa <- dudi.pco(dst, scan=FALSE, nf=3)
  #print(pcoa$li)
  loadings <- signif((pcoa$eig)[1:3] / sum(pcoa$eig)*100, digits=3)
  
  colnames(pcoa$li)<-c("PC1","PC2","PC3")
  data<-as.data.frame(pcoa$li)
  
  # write axes file
 # write.table(data,file=axesfile,sep="\t",quote=FALSE,col.names=NA)
  
  data$group<-rownames(data)
  data$group2<-data$group
  meta$SampleID_1<-rownames(meta)
  data$group2<-data$group
  data$SampleID_1<-data$group
  data<-plyr::join(data,meta)
  
 data_1<-dplyr::filter(data,Method!="METGEM"& Method!="METGEM-CUS")


pdf("F:/Cohort/github/Figures/Simulation_PCA_P1P2_1.pdf",width=80, height=60)
  ggplot() +
    geom_point(data=data,aes(x=PC1,y=PC2,color= Method,shape=Nature),stat='identity',size=30) +
    xlab(paste0("PC2: ",loadings[2],"% variance")) +
    ylab(paste0("PC3: ",loadings[3],"% variance")) + 
    theme(axis.text.x=element_text(size=100,colour="black"),
          axis.text.y=element_text(size=100,colour="black"),
          axis.title.x=element_text(size=100),
          axis.title.y=element_text(size=100),
          legend.key=element_rect(colour="black",size=1),
          legend.text = element_text(size=100),
          legend.title =  element_text(size=100),
          panel.border=element_rect(fill=NA),
          panel.grid.minor=element_blank(),
          panel.background=element_blank())
  dev.off()
  
  ggplot() +
    geom_point(data=data,aes(x=PC1,y=PC2,color= Method),stat='identity',size=3) +
    theme(axis.text.x=element_text(size=12,colour="black"),
          axis.text.y=element_text(size=12,colour="black"),
          axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),
          legend.key=element_rect(colour="black",size=0.2),
          panel.border=element_rect(fill=NA),
          panel.grid.minor=element_blank(),
          panel.background=element_blank())+
    scale_x_continuous(limits = c(-0.19, -0.14))
  
  
 
