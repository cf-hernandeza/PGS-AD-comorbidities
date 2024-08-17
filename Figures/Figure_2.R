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
 
############################### Figure 2 , Comorbidity-PGS in Alzheimer’s disease Patients and general population of UK Biobank.

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
###### Plot PRS Violin Plots ####
Box_plot <- function (PRS,Disease,subset,age=0,test="t.test") {
    DATA <- load_data(PRS)
    subset_carga <- as.character(subset)
    Disease <- as.character(Disease)
    
    if (subset == "AD") {
        ### Filter AD 
        data_AD <- DATA[DATA[[subset_carga]] == 1 ]
        data_CMRB <- DATA[DATA[[Disease]] == 0 ] 
        data_subset <- merge(data_AD,data_CMRB,all = TRUE)
    } else {
        ### Filter CMRB 
        data_subset <- DATA[DATA[[Disease]] == 1 ]
        
    }
    
    
    
    
    
    
    
    ##Filter by Age
    data_test <- data_subset[data_subset$AGE >= age]
    ## Generate 3 cohorts and count
    if (subset == "AD") {
        ### with and withou cmrb 
        data_test$PHENO_PLOT <- as.factor(ifelse(data_test[[Disease]] == 1 & data_test$AD == 1  , 2, 
                                                 ifelse(data_test[[Disease]] == 0 & data_test$AD == 1 , 1, 0       )) )
    } else {
        ### without AD 
        data_test$PHENO_PLOT <- as.factor(ifelse(data_test$AD == 1 , 0, 1))
        
    }
    N_records <- data_test %>% 
        group_by(PHENO_PLOT) %>% 
        summarise(count = n(),MEAN=mean(zSCORE), MAX=max(zSCORE),MIN=min(zSCORE),SD= sd(zSCORE) ) %>% 
        mutate(UP= MEAN + 1.96*SD   ,LOW=MEAN - 1.96*SD)
    
    
    
    ## pvalues 
    if (test == "t.test") {
        stat.test <- data_test %>% t_test(zSCORE ~ PHENO_PLOT)
        stat.test <- stat.test %>% add_xy_position(x = "PHENO_PLOT")
        stat.test$p.adj <- stat.test$p
        
    } else if (test == "anova"){
        stat.test <- data_test %>% tukey_hsd(zSCORE ~ PHENO_PLOT)
        stat.test <- stat.test %>% add_xy_position(x = "PHENO_PLOT")
    } 
    
    
    ##plot Violin Plot
    mypal <- brewer.pal(3, "Dark2")[1]
    mypal2 <- brewer.pal(3, "Dark2")[2]
    mypal3 <- brewer.pal(3, "Dark2")[3]
    
    
    p <- ggplot(N_records,aes(x = PHENO_PLOT ,
                              y = MEAN  ) ) +
        scale_x_discrete(
            limits=c("0","1","2"),
            breaks=c("0","1","2"),
            labels= c("GP" ,subset, paste0("AD with\n",Disease) )
        )
    
    p
    p2 <- p+
        geom_text(aes(label=paste0("N=",count)),
                  y=min(N_records$LOW),vjust=2 ,hjust=0.5 ,size=6 ) +
        
        theme_minimal() +
        stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.05, 
                           step.increase = 0.07,
                           size = 6) +
        ylim(min(N_records$LOW) -4,
             max(N_records$UP) + 8)
    p2
    p3 <- p2 + scale_fill_manual(values=c(mypal3,mypal,mypal2) ) +
        theme_bw() + 
        ylab(paste0(PRS," PGS (zSCORE)"))  +
        #xlab(paste0(ifelse(age != 0 , paste0("Patients Over ",age," Years") ,""))) + 
        xlab("") +
        theme(legend.position = "none") + 
        theme(axis.text.x = element_text(size=20),
              axis.title = element_text(size=20),
              axis.text.y= element_text(size=15))
    
    p3
    p3 + geom_boxplot(data = data_test,aes(x= PHENO_PLOT,
                                           y= zSCORE,
                                           fill= PHENO_PLOT),
                      outlier.shape = NA,
                      width=0.4)
}

#make Plots
p1 <- Box_plot("T2D","T2D","AD",test ="anova",age=65) %>% label_plot(.,"A",25)
p2 <- Box_plot("MDD","MDD","AD",test ="anova",age=65) %>% label_plot(.,"B",25)
p3 <- Box_plot("MH","MH","AD",test ="anova",age=65)  %>% label_plot(.,"C",25)
p4 <- Box_plot("EPI","EPI","AD",test = "anova",age=65)  %>% label_plot(.,"D",25)


##Arrenge in grid
lay_vertical <- rbind(c(1,2),
             c(3,4))

grid_vertical <- arrangeGrob(grobs = list(p1,p2,p3,p4), layout_matrix = lay_vertical)

### PDF Horizontal 20*11
ggsave(file="All_figures/Figure2.pdf",grid_vertical,width = 11,height = 11)
