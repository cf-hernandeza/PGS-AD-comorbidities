#!/usr/bin/env Rscript
####USAGE######
# ./4_PCA_onlt_ukb.r [Comorbidity] {PCA.tsv} 
#
#
#
#
args <- commandArgs(trailingOnly = TRUE)

### loa libraries
library(data.table)
library(ggplot2)
library(dplyr)


CMRB=args[1]
PCA=args[2]

# load PCA 
data <- fread(PCA)
colnames(data) <- c("IID","PCA1","PCA2","PCA3","PCA4") #renombrar columnas

# Graficar PC1 vs PC2
jpeg("PRS_EPI_all_2022/4_Thresholding/Sample_QC/PCA_results/plot_p1vsp2.jpeg")
    ggplot(data)+
        geom_point(aes(PCA1,PCA2))
dev.off()

# Graficar PC2 vs PC3
jpeg("PRS_EPI_all_2022/4_Thresholding/Sample_QC/PCA_results/plot_p2vsP3.jpeg")
    ggplot(data)+
        geom_point(aes(PCA2,PCA3))
dev.off()

# Graficar PC3 vs PC4 
jpeg("PRS_EPI_all_2022/4_Thresholding/Sample_QC/PCA_results/plot_p3vsP4.jpeg")
    ggplot(data)+
        geom_point(aes(PCA3,PCA4))
dev.off()



jpeg("PRS_EPI_all_2022/4_Thresholding/Sample_QC/PCA_results/eur_group.jpeg")
     ggplot(data,aes(PCA1,PCA2)) +
        geom_point() +
        xlab("PC1") +
        ylab("PC2") +
        ggtitle("EUR group")+
        geom_segment(aes(x = -3.7, y = -3.7, xend = -3.7, yend = 3.7))+
        geom_segment(aes(x = -3.7, y = -3.7, xend = 3.7, yend = -3.7))+
        geom_segment(aes(x = 3.7, y = 3.7, xend = -3.7, yend = 3.7))+
        geom_segment(aes(x = 3.7, y = 3.7, xend = 3.7, yend = -3.7))
dev.off()

temp1 <- data[(data$PCA1 >= -3.7),]
temp2 <- temp1[(temp1$PCA1 <= 3.7),]
temp3 <- temp2[(temp2$PCA2 >= -3.7),]
temp4 <- temp3[(temp3$PCA2 <= 3.7),]

jpeg("PRS_EPI_all_2022/4_Thresholding/Sample_QC/PCA_results/eur_group_filter.jpeg")
     ggplot(temp4,aes(PCA1,PCA2)) +
        geom_point() +
        xlab("PC1") +
        ylab("PC2") +
        ggtitle("EUR group")+
        geom_segment(aes(x = -3.7, y = -3.7, xend = -3.7, yend = 3.7))+
        geom_segment(aes(x = -3.7, y = -3.7, xend = 3.7, yend = -3.7))+
        geom_segment(aes(x = 3.7, y = 3.7, xend = -3.7, yend = 3.7))+
        geom_segment(aes(x = 3.7, y = 3.7, xend = 3.7, yend = -3.7))
dev.off()

filtered_ukb <- temp4$IID

write.table(filtered_ukb,file = "PRS_EPI_all_2022/4_Thresholding/Sample_QC/keep_ancestry_eur.txt", quote = FALSE,row.names=F,col.names = F)

