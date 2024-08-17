#!/home/carlos/miniconda3/envs/r4.1/bin/Rscript --vanilla
#
#
#
#
library(data.table)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(magrittr)
library(GGally)
library(cowplot)
library(magrittr)
############################### Figure sup 2 , Relationship between comorbidity-PGSs

##define functions

#Load data
load_data <- function (PRS) {
  data <- fread(paste0(PRS,"_PRS_pheno_covar.tsv"))
  names(data)[3] <- "zSCORE"
  data
}
#Label Plot
label_plot <- function (plot,label,size) {
    plot + labs(tags=label) + theme(plot.tag = element_text(size=size, face= "bold") )
}
#plot PRS VS PRS
PRS_PRS <- function(){
    load_PRS_data <- function (PRS) {
        data <- fread(paste0(PRS,"_PRS_pheno_covar.tsv"))[,c(1,3)]
        data
    }
    data <- merge(merge(load_PRS_data("EPI"),load_PRS_data("T2D"),by="FID"),
                  fread("EPI_PRS_pheno_covar.tsv")[,c(1,14,24)], by="FID" )
    data <- data[data$AGE >= 65]
    
    data_AD <- data[data$AD == 1 ]
    colnames(data_AD) <- c("FID","PGS EPI (zSCORE)","PGS T2D (zSCORE)","AGE","AD")
    
    
    plot2 <- ggpairs(data_AD, columns= c("PGS T2D (zSCORE)","PGS EPI (zSCORE)"), 
                     title = "",
                     upper = list(continuous = wrap("cor",color="black", size = 5)),
                     lower = list(continuous = wrap("smooth",color="gray", size = 2.5))) +
        theme_bw() 
    p2 <- label_plot(plot2,"B",25)
    p2
}
#Make_plot
p1 <- PRS_PRS()


ggsave(file="All_figures/FigureSup2.pdf",p1,width = 5,height = 4)



