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
 
############################### Figure 1 , Distribution Comorbidity-PGS in Alzheimer’s disease Patients of UK Biobank

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

###### Plot PRS Histogram ####
PRS_histogram <- function(PRS,Disease,subset=FALSE,age=0,test="ks.test"){
    DATA <- load_data(PRS)
    subset <- as.character(subset)
    Disease <- as.character(Disease)
    
    ### Filter  
    data_subset <- DATA[DATA[[subset]] ==1 ]
    
    
    
    ##Filter by Age
    data <- data_subset[data_subset$AGE >= age]
    data[[Disease]] <- as.character(data[[Disease]])
    data[[subset]] <-as.character(data[[subset]])
    
    TEST <-  ks.test(x= data$zSCORE[data[[Disease]] == 0]  ,
                     y= data$zSCORE[data[[Disease]] == 1] )
    stat.test <- c(pvalue=signif(TEST[[2]],3))
    histogram1 <- density(data$zSCORE[data[[Disease]] == 0])
    histogram2 <- density(data$zSCORE[data[[Disease]] == 1])
    
    max_height1 <- max(histogram1$y)
    max_height2 <- max(histogram2$y)
    maxvalue <- ifelse( max_height1 > max_height2,
                        max_height1,
                        max_height2)
    
    
    
    #plot
    p <- ggplot(data,aes_string("zSCORE", fill=Disease  )) + 
        geom_density(weight=10 ,alpha= 0.5 ) +
        geom_vline(xintercept = mean(data$zSCORE[data[[Disease]]==0]),
                   alpha=1, linewidth=1,
                   linetype="longdash",
                   color=ifelse(subset=="AD",
                                "#1b9e77",
                                "orange") ) +
        geom_vline(xintercept = mean(data$zSCORE[data[[Disease]]==1]),
                   alpha=1, linewidth=1,
                   linetype="longdash",color="#d95f02" ) +
        annotate("text", x = 2, y = maxvalue, label = stat.test, size=6)
    p2 <-  p + 
        ylab(ifelse(subset==Disease, "General Population (%)",paste0(subset," Subset (%)")))  +
        xlab(paste0(PRS," PGS (zSCORE)")) +
        scale_fill_manual(name = "Status",
                          values=c(ifelse(subset=="AD",
                                          "#1b9e77",
                                          "orange"),
                                   "#d95f02"),
                          labels= c(ifelse(subset=="AD",
                                           subset,
                                           paste0("Only\n",subset)),
                                    paste0("AD with\n",Disease))
        )  +
        labs(fill="Status") +
        theme_bw() + 
        theme(legend.position = c(0.20,0.80),
              legend.title = element_text(color = "black", size = 20,),
              legend.text = element_text(color = "black", size = 15),
              legend.background = element_rect(fill = "white", colour = "gray30"),
              legend.key.height=unit(1, "cm"),
              legend.key.width = unit(1,"cm"),
              axis.text = element_text(size=20),
              axis.title = element_text(size=20),) 
    
    p2
}

#make Plots
p1 <- PRS_histogram("T2D","T2D","AD",age=65)  %>% label_plot(.,"A",25) 
p2 <- PRS_histogram("MDD","MDD","AD",age=65)  %>% label_plot(.,"B",25)
p3 <- PRS_histogram("MH","MH","AD",age=65)  %>% label_plot(.,"C",25)
p4 <- PRS_histogram("EPI","EPI","AD",age=65) %>% label_plot(.,"D",25)


##Arrenge in grid
lay_vertical <- rbind(c(1,2),
                      c(3,4))

grid_vertical <- arrangeGrob(grobs = list(p1,p2,p3,p4), layout_matrix = lay_vertical)

### PDF Horizontal 20*11
ggsave(file="All_figures/Figure1.pdf",grid_vertical,width = 11,height = 11)
