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
library(GGally)
library(cowplot)
library(magrittr)

############################### Figure 3 , Association of comorbidity-PGSs with AD onset in individuals from the UK Biobank

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
###### Plot PRS Violin Plots
PRS_OR <- function(PRS,Disease,subset=FALSE,age=50){
  DATA <- load_data(PRS)
  subset_carga <- as.character(subset)
  if (subset == FALSE) {
    data_test <- DATA[DATA$PD == 0 ]
  } else {
      data_test <- DATA[DATA[[subset_carga]] == 1 ]
  } 
  ##Filter by Age
  data <- data_test[data_test$AGE >= age]
  
  if (subset == "PD") {
      
  ## cargar AAO PD
  AAO <-  load_data(PRS)[,c(1,25,26)]  
  colnames(AAO) <- c("FID","Birth","Year_PD")
  AAO$Year_PD <- as.numeric(format(AAO$Year_PD,'%Y'))
  AAO$AAO_PD <- AAO$Year_PD - AAO$Birth
  
  ## merge data
  data_PD_AAO_all <- merge(data,AAO, by="FID")
  data_PD_AAO_all$AGE_plot <-  data_PD_AAO_all$AAO_PD - data_PD_AAO_all$AGE
  min_age <- abs(min( data_PD_AAO_all$AGE_plot,na.rm=TRUE))
  data_PD_AAO_all$AGE_plot_Positive <-  data_PD_AAO_all$AGE_plot + min_age
  data_PD_AAO_all <- data_PD_AAO_all[!is.na(data_PD_AAO_all$AGE_plot_Positive)]
  data_PD_AAO <- subset(data_PD_AAO_all, data_PD_AAO_all$AAO_PD >= 50 )
  
  data <- data_PD_AAO
  }
  data$SCORE_20_80 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.8), 0, 1) 
  data$SCORE_10_90 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.9), 0, 1) 
  data$SCORE_5_95 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.95), 0, 1)
  data$SCORE_1_99 <- ifelse(data$zSCORE <= quantile(data$zSCORE,0.99), 0, 1)
  
  values_tresh <- data.frame(NULL)
  values_tresh <- cbind(Percentil=c(20,10,5,1))
  tresh <- c()
  tresh <-c(tresh, as.numeric(min(data[data$SCORE_20_80 == 1 ]$zSCORE )))
  tresh <-c(tresh, as.numeric(min(data[data$SCORE_10_90 == 1 ]$zSCORE )))
  tresh <-c(tresh, as.numeric(min(data[data$SCORE_5_95 == 1 ]$zSCORE )))
  tresh <-c(tresh, as.numeric(min(data[data$SCORE_1_99 == 1 ]$zSCORE )))
  values_tresh <- cbind(values_tresh,tresh)
  COEFFS <- data.frame(NULL)
  for (i in c(20,10,5,1)) {
    # Data pvalue
    formula <- as.formula(paste0(Disease, " ~ ", paste0("SCORE_",i,"_",100-i), " + SEX+ AGE +PC1+PC2+PC3+PC4 "))
    glm_tresh <- glm(formula,family = "binomial",data=data )
    BETA <-  summary(glm_tresh)$coeff[2,1]
    SE <- summary(glm_tresh)$coeff[2,2]
    p.value <- summary(glm_tresh)$coeff[2,4]
    valores <- data.frame(BETA,SE,p.value)
    
    #numero casos controles
    var <- paste0("SCORE_",i,"_",100-i)
    case_control <- data %>% 
    group_by(.data[[var]]) %>% 
        summarise(Total=n(),Cases=sum( .data[[Disease]] ),Controls= (Total - Cases))

    Table <- data.frame(NULL)
    Table <- rbind(Table,c(ifelse(subset==FALSE,"GP",subset),
                               Disease,
                               paste0(case_control[1,3],"/",case_control[1,4]),
                               paste0(case_control[2,3],"/",case_control[2,4])) )
        
    
    names(Table) <- c("Population","CMRB","Low_PRS","High_PRS")
    Table
    Row <- cbind(valores,Table)
    COEFFS <- rbind(COEFFS,Row)
  }
  values_tresh <- cbind(values_tresh,COEFFS)
  names(values_tresh) <- c("Percentil","Tresh","BETA","SE","p.value","Population","CMRB","Low_PRS","High_PRS")
  values_tresh$BETA_low <- values_tresh$BETA - 1.96 * values_tresh$SE  
  values_tresh$BETA_high <- values_tresh$BETA + 1.96 * values_tresh$SE  
  values_tresh$OR <- signif(exp(values_tresh$BETA),3)  
  values_tresh$OR_low <-  signif(exp(values_tresh$BETA_low),3)
  values_tresh$OR_high <-  signif(exp(values_tresh$BETA_high),3)
  
  values_tresh$OR_signif <- signif(values_tresh$OR,2)
  values_tresh$CMRB <- rep(Disease,4)
  values_tresh$Population <- rep(ifelse(subset==FALSE,"GP",subset),4)
  
  values_tresh
}
###Plot PRSs
Onset_curves_prevalences <- function(PRS,Disease,subset="AD",age=0,onset=0){
    DATA <- load_data(PRS)
    subset <- as.character(subset)
    if (subset == "Disease") {
        data_test <- DATA[DATA[[Disease]] == 1 ]
    } else if (subset== "AD") {
        data_test <- DATA[DATA$AD == 1 ]
    }
    
    ##Filter by Age
    data <- data_test[data_test$AGE >= age]
    
    ## cargar AAO AD
    #AAO <-  fread ("../Data/CLINICAL_COVAR_QC_participant.tsv",sep="\t")
    #which(str_detect(colnames(AAO), regex("Date G40|Date G20|Date E11|Date F32|Date G470|Date G43|Date G47")))
    AAO <-  fread ("../Data/CLINICAL_COVAR_QC_participant.tsv",sep="\t")[,c(1,5,62,68,74,76,86,80)]  
    colnames(AAO) <- c("FID","Birth","Date_T2D","Date_MDD","Date_AD","Date_EPI","Date_MH","Date_Insomnia")
    
    #AAO AD
    AAO$Year_AD <- as.numeric(format(AAO$Date_AD,'%Y'))
    AAO$AAO_AD <- AAO$Year_AD - AAO$Birth
    #AAO CMRB
    AAO$Date_CMRB <- ifelse(str_detect(AAO[[paste0("Date_",Disease)]],"Code"),
                            "",
                            ifelse(is.na(AAO[[paste0("Date_",Disease)]]),
                                   "",
                                   AAO[[paste0("Date_",Disease)]]))
    AAO$Date_CMRB <- as.Date(AAO$Date_CMRB)
    AAO$Year_CMRB <- as.numeric(format(AAO$Date_CMRB,'%Y'))
    AAO$AAO_CMRB <- AAO$Year_CMRB - AAO$Birth
    
    ## merge data
    data_AD_AAO_all <- merge(data,AAO, by="FID")
    
    if (subset == "AD") {
        data_AD_AAO_all$AGE_plot <-  data_AD_AAO_all$AAO_AD
    }else if (subset == "Disease"){
        data_AD_AAO_all$AGE_plot <-  data_AD_AAO_all$AAO_CMRB
    }
    data_AD_AAO <- subset(data_AD_AAO_all, data_AD_AAO_all$AGE_plot >= onset )
    
    
    #Deciles
    data_AD_AAO$dec1 <- ifelse(data_AD_AAO$zSCORE <= quantile(data_AD_AAO$zSCORE,0.20), 1, 0)
    data_AD_AAO$dec2_9 <- ifelse(data_AD_AAO$zSCORE > quantile(data_AD_AAO$zSCORE,0.33) & data_AD_AAO$zSCORE <= quantile(data_AD_AAO$zSCORE,0.66), 1, 0)
    data_AD_AAO$dec10 <- ifelse(data_AD_AAO$zSCORE > quantile(data_AD_AAO$zSCORE,0.80), 1, 0)
    
    #agrupate in one varaible
    data_AD_AAO$dec[data_AD_AAO$dec1 == 1] <- "low"
    #data_AD_AAO$dec[data_AD_AAO$dec2_9 == 1] <- "medium"
    data_AD_AAO$dec[data_AD_AAO$dec10 == 1] <- "high"
    
    #phantom_data <- data_AD_AAO[1,]
    #phantom_data$MH <- 0
    #phantom_data$AAO_AD <- 30
    #phantom_data_1 <- phantom_data
    #phantom_data_2 <- phantom_data
    #phantom_data_3 <- phantom_data
    #phantom_data_1$dec <- "high"
    #phantom_data_2$dec <- "medium"
    #phantom_data_3$dec <- "low"
    
    #data_AD_AAO_with0 <- rbind(data_AD_AAO,phantom_data_1,phantom_data_3) 
    data_AD_AAO_with0 <- data_AD_AAO
    data_AD_AAO_with0 <- data_AD_AAO_with0[!is.na(data_AD_AAO_with0$dec)]
    total_low <- sum(data_AD_AAO_with0$dec == "low")
    total_high <-sum(data_AD_AAO_with0$dec == "high")
    
    if (subset== "AD"){
        data_plot <- data_AD_AAO_with0 %>% 
            select(FID,zSCORE,Disease,AGE,AD,AGE_plot,dec) %>% 
            group_by(dec,AGE_plot) %>% 
            summarise(cases=sum(!!as.symbol(  Disease) )  )
    } else if (subset== "Disease" ) {
        data_plot <- data_AD_AAO_with0 %>% 
            select(FID,zSCORE,AD,AGE,AGE_plot,dec) %>% 
            group_by(dec,AGE_plot) %>% 
            summarise(cases=sum(AD )  )
    }
    data_plot$N_total <- ifelse(data_plot$dec == "low", total_low, total_high ) 
    
    data_plot <- data_plot %>%
        mutate (prevalence=(cases/N_total)*100,
                accumulate = cumsum(prevalence))
    
    
    data_plot
    
    
    #plot
    mypal <- brewer.pal(3, "Set1")[1]
    mypal2 <- brewer.pal(3, "Set1")[2]
    
    
    plot1 <- ggplot(data_plot,aes(AGE_plot, accumulate, color=dec)) +
        geom_point() +
        geom_line() +
        theme_bw() + 
        labs(colour=paste("PGS Groups"))  +
        scale_color_manual(values = c("high" = mypal, "low" = mypal2),
                           limits = c("high", "low"),
                           breaks = c("high", "low"),
                           labels = c("High", "Low")) +
        theme_pubclean()
    plot1
    plot2 <- plot1 + theme(legend.title = element_text(color = "black", size = 20),
                           legend.text = element_text(color = "black", size = 15),
                           legend.background = element_rect(fill = "white", colour = "gray30"),
                           axis.text.y = element_text(size=20),
                           axis.text.x = element_text(size=20) ,
                           axis.title = element_text(size=25))
    plot2
    if (subset == "AD") {
        plot3 <- plot2 + 
            xlab("Age At Onset of AD (Years)")  +
            ylab(paste0("% Prevalence of ", Disease))     
    } else if (subset == "Disease") {
        plot3 <- plot2 + 
            xlab(paste0("Age At Onset of ",Disease ,"(Years)") ) +
            ylab(paste0("% Prevalence of AD"))   
    }
    plot3
}

#make plots
p1 <- Onset_curves_prevalences("T2D","T2D",subset="AD",age=65) %>% label_plot(., "A",25)
p2 <- Onset_curves_prevalences("EPI","EPI",subset="AD",age=65)  %>% label_plot(., "B",25)
  
##Arrenge in grid
lay_horizontal <- rbind(c(1,2))
grid_horizontal <- arrangeGrob(grobs = list(p1,p2), layout_matrix = lay_vertical)

#Save Plots
ggsave(file="All_figures/Figure3.pdf",grid_horizontal,width = 12,height = 6.5)
  
