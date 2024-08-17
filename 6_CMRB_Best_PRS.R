#!/usr/bin/env Rscript
#usage 
#./6_CMRB_Best_PRS.R {used_tresh.txt} {final.Phenotipo_QC.fam} {COVARS.tsv } {TRESH VALUE (Run with this empty the firs time)}
## 
arg <- commandArgs(trailingOnly = TRUE)

## Definir Epilepsia
#print(cmrb)
#definir fam
#fam <- args[1]
#print(fam)
#soft <- args[3]
#print(soft)
#suffix <- args[4]
#print(suffix)
PATH_TO_CLUMP <- "PRS_EPI_all_2022/3_Clump_PRS/"
PATH_TO_THRESH <- "PRS_EPI_all_2022/4_Thresholding/"

TRESH_FILE <- arg[1]
FAM_FILE <- arg[2] 
COVAR_FILE <- arg[3]
tresh <- arg[4]

#Librerias a Usar
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(rms)))
suppressWarnings(suppressMessages(library("pROC")))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(stringr)))

if (is.na(tresh) ){
  if(! file.exists(paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_results.txt"))){
    #treshold usados
    cat("Calculating Best PRS \n"  )
    p.threshold <- read.table(TRESH_FILE, colClasses = "character")
    p.threshold <- p.threshold$V1
    # Read in the fam file 
    phenotype <- read.table(FAM_FILE, header=FALSE)[,c(1,2,6)]
    colnames(phenotype) <- c("FID", "IID","PHENO")
    # Read in the PCs
    pcs <- read.table(COVAR_FILE, header=T)[,c(1,1,5:8)] #borrar para usar todas las PC
    # The default output from plink does not include a header
    # To make things simple, we will add the appropriate headers
    colnames(pcs) <- c("FID", "IID", paste0("PC",1:4))  ###cambair a 8 para todas las PC
    # Read in the covariates (here, it is age)
    covariate <- fread(COVAR_FILE, header=T)[,c(1,3,2,4)] #añadir 6 para TD
    colnames(covariate) <- c("FID", "AGE", "SEX","TD") #añadir TD
    covariate$IID <- covariate$FID
    # read covariates PD EPI
    covariate_diagnosis <- fread(COVAR_FILE, header=T)[,c(1,9)]
    colnames(covariate_diagnosis) <- c("FID", "PD")
    covariate_diagnosis$IID <- covariate_diagnosis$FID

           
    # Now merge the files
    pheno_2 <- merge(merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID")),covariate_diagnosis,by=c("FID","IID"))
    pheno_1 <- pheno_2[pheno_2$PHENO != -9,]
    
     
    
    all_rows <- nrow(pheno_1)
    training <- all_rows * 0.8
    set.seed(123456)
    pheno <- sample_n(pheno_1, training )
    pheno_validation <- subset(pheno_1, !(FID %in% pheno$FID))
    pheno_validation$PHENO <- pheno_validation$PHENO - 1 
    
    #recodificar casos y controles
    pheno$PHENO <- pheno$PHENO - 1 
    # We can then calculate the null model (model with PRS) using a logistic regression 
    # (as pheno is cualitative)
    null.model_lrm <- lrm(PHENO ~ 1 + AGE + SEX + PC1+ PC2+ PC3 +PC4 + TD, data=pheno)
    null.model_glm <- glm(PHENO ~ 1 + AGE + SEX + PC1+ PC2+ PC3 +PC4 + TD, data=pheno, family="binomial")
    #summary(null.model_glm)
    #summary model 
    R2.null <- null.model_lrm$stats["R2"]
    fwrite(pheno_1 ,paste0(PATH_TO_THRESH,"EPI_all_2022_covar.txt"),sep='\t')
    
    
    prs.result <- NULL
    for(i in p.threshold){
      # Go through each p-value threshold
      prs <- read.table(paste0(PATH_TO_CLUMP,"SS_PLINK1.9_PRS.",i,".profile"), header=T)[,c(2,2,4,5,6)]
      colnames(prs)<- c("FID","IID","Total","N","SCORE")
      # Merge the prs with the phenotype matrix
      # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
      # relevant columns
      pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
      #setting SCORE to Zscore
      meanControls <- mean(pheno.prs$SCORE[pheno.prs$PHENO == 0])
      sdControls <- sd(pheno.prs$SCORE[pheno.prs$PHENO == 0])
      pheno.prs$zSCORE <- (pheno.prs$SCORE - meanControls)/sdControls  
      pheno_validation.prs <- merge(pheno_validation, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
      #setting SCORE to Zscore
      meanControls <- mean(pheno_validation.prs$SCORE[pheno_validation.prs$PHENO == 0])
      sdControls <- sd(pheno_validation.prs$SCORE[pheno_validation.prs$PHENO == 0])
      pheno_validation.prs$zSCORE <- (pheno_validation.prs$SCORE - meanControls)/sdControls 
      # Now perform a logistic regression on pheno with PRS and the covariates
      # ignoring the FID and IID from our model
      model_glm <- glm(PHENO ~ zSCORE + AGE + SEX + PC1+ PC2+ PC3 +PC4 +TD , data=pheno.prs, family = "binomial")
      model_lrm <- lrm(PHENO ~ zSCORE + AGE + SEX + PC1+ PC2+ PC3 +PC4 +TD , data=pheno.prs)
      #summary model
      # model AUC is obtained as 
      predicted_pheno <- predict(model_glm, newdata= pheno_validation.prs, type= "response")
      AUC <- roc(pheno_validation$PHENO, predicted_pheno)
      AUC.model <- AUC$auc
      R2.model <- model_lrm$stats["R2"]
  
      # AUC of PRS is simply calculated as the model AUC minus the null AUC
      prs.AUC <- AUC.model
      prs.R2 <- R2.model - R2.null
      # We can also obtain more statistics
      prs.deviance <- summary(model_glm)$deviance
      prs.aic <- summary(model_glm)$aic
      prs.pvalue <-anova(model_lrm)["zSCORE",3] #cambiar a 13 TD -  12 para 8PC 
      # We can then store the results
      prs.result <- rbind(prs.result, data.frame(Threshold=i, AUC=prs.AUC, deviance=prs.deviance, AIC=prs.aic,R2_diff=prs.R2,R2=R2.model ,Pvalue=prs.pvalue))
    }
    prs.result$print.p <- round(prs.result$Pvalue, digits = 3)
    prs.result$print.p[!is.na(prs.result$print.p) &
                         prs.result$print.p == 0] <-
      format(prs.result$Pvalue[!is.na(prs.result$print.p) &
                                 prs.result$print.p == 0], digits = 2)
    prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
    fwrite(prs.result, paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_results.txt"),sep="\t")
    }
  else {
    prs.result <- fread(paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_results.txt"))
    
  
  }
  ### Validation
  if(! file.exists(paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_results_validation.txt"))){
    #treshold usados
      p.threshold <- read.table(TRESH_FILE,colClasses = "character")
      p.threshold <- p.threshold$V1    # Read in the fam file 
    phenotype <- read.table(FAM_FILE, header=FALSE)[,c(1,2,6)]
    colnames(phenotype) <- c("FID", "IID","PHENO")
    # Read in the PCs
    pcs <- read.table(COVAR_FILE, header=T)[,c(1,1,5:8)] #borrar para usar todas las PC
    # The default output from plink does not include a header
    # To make things simple, we will add the appropriate headers
    colnames(pcs) <- c("FID", "IID", paste0("PC",1:4))  ###cambair a 8 para todas las PC
    # Read in the covariates (here, it is age)
    covariate <- fread(COVAR_FILE, header=T)[,c(1,3,2,4)] #añadir 6 para TD
    colnames(covariate) <- c("FID", "AGE", "SEX","TD") #añadir TD
    covariate$IID <- covariate$FID
    # read covariates PD EPI
    covariate_diagnosis <- fread(COVAR_FILE, header=T)[,c(1,9)]
    colnames(covariate_diagnosis) <- c("FID", "PD")
    covariate_diagnosis$IID <- covariate_diagnosis$FID
    
    
    # Now merge the files
    pheno_2 <- merge(merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID")),covariate_diagnosis,by=c("FID","IID"))
    pheno_1 <- pheno_2[pheno_2$PHENO != -9,]
    
    
    
    all_rows <- nrow(pheno_1)
    training <- all_rows * 0.8
    set.seed(123456)
    pheno <- sample_n(pheno_1, training )
    pheno_validation <- subset(pheno_1, !(FID %in% pheno$FID))
    pheno_validation$PHENO <- pheno_validation$PHENO - 1 
    
    #recodificar casos y controles
    pheno$PHENO <- pheno$PHENO - 1 
    # We can then calculate the null model (model with PRS) using a logistic regression 
    # (as pheno is cualitative)
    null.model_lrm_validation <- lrm(PHENO ~ 1 + AGE + SEX + PC1+ PC2+ PC3 +PC4 + TD, data=pheno_validation)
    null.model_glm_validation <- glm(PHENO ~ 1 + AGE + SEX + PC1+ PC2+ PC3 +PC4 + TD, data=pheno_validation, family="binomial")
    #summary(null.model_glm)
    #summary model 
    R2.null_validation <- null.model_lrm_validation$stats["R2"]
    fwrite(pheno_1 ,paste0(PATH_TO_THRESH,"EPI_all_2022_covar_validation.txt"),sep='\t')
    
    
    prs.result_validation <- NULL
    for(i in p.threshold){
        # Go through each p-value threshold
      prs <- read.table(paste0(PATH_TO_CLUMP,"SS_PLINK1.9_PRS.",i,".profile"), header=T)[,c(2,2,4,5,6)]
      colnames(prs)<- c("FID","IID","Total","N","SCORE")
      # Merge the prs with the phenotype matrix
      # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
      # relevant columns
      pheno.prs <- merge(pheno_validation, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
      #setting SCORE to Zscore
      meanControls <- mean(pheno.prs$SCORE[pheno.prs$PHENO == 0])
      sdControls <- sd(pheno.prs$SCORE[pheno.prs$PHENO == 0])
      pheno.prs$zSCORE <- (pheno.prs$SCORE - meanControls)/sdControls  
            # Now perform a logistic regression on pheno with PRS and the covariates
      # ignoring the FID and IID from our model
      model_glm <- glm(PHENO ~ zSCORE + AGE + SEX + PC1+ PC2+ PC3 +PC4 +TD , data=pheno.prs, family = "binomial")
      model_lrm <- lrm(PHENO ~ zSCORE + AGE + SEX + PC1+ PC2+ PC3 +PC4 +TD , data=pheno.prs)
      #summary model
      # model AUC is obtained as 
      R2.model_validation <- model_lrm$stats["R2"]
      
      # AUC of PRS is simply calculated as the model AUC minus the null AUC
      prs.R2 <- R2.model_validation - R2.null_validation
      # We can also obtain more statistics
      prs.deviance <- summary(model_glm)$deviance
      prs.aic <- summary(model_glm)$aic
      prs.pvalue <-anova(model_lrm)["zSCORE",3] #cambiar a 13 TD -  12 para 8PC 
      # We can then store the results
      prs.result_validation <- rbind(prs.result_validation, data.frame(Threshold=i, R2_diff=prs.R2,R2=R2.model_validation ,Pvalue=prs.pvalue))
    }
    prs.result_validation$print.p <- round(prs.result_validation$Pvalue, digits = 3)
    prs.result_validation$print.p[!is.na(prs.result_validation$print.p) &
                                    prs.result_validation$print.p == 0] <-
      format(prs.result_validation$Pvalue[!is.na(prs.result_validation$print.p) &
                                            prs.result_validation$print.p == 0], digits = 2)
    prs.result_validation$print.p <- sub("e", "*x*10^", prs.result_validation$print.p)
    fwrite(prs.result_validation, paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_results_validation.txt"),sep="\t")
  }
  else {
    prs.result_validation <- fread(paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_results_validation.txt"))
    
    
  }
  
# Best result is:
max_r2 <- prs.result[which.max(prs.result$R2_diff),]
max_r2_validation <- prs.result_validation[which.max(prs.result_validation$R2_diff),]

cat('Best P value treshold is:\n ')
print(max_r2)
print(max_r2_validation)





# Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
p <- ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
    geom_text(
        aes(label = paste(print.p)),
        vjust = -1.5,
        hjust = 0,
        angle = 45,
        cex = 4,
        parse = T
    )  +
    scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
    xlab(expression(italic(P) - value ~ threshold ~ (italic(Pvalue)[T]))) +
    ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
    geom_bar(aes(fill = -log10(Pvalue)), stat = "identity") +
    scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-4,
        name = bquote(atop(-log[10] ~ model, italic(Pvalue) - value),)) +
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size =
                                        18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust =
                                    1))
ggsave(plot=p , paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_eval_tresholds.pdf"))

#Validation
p1 <- ggplot(data = prs.result_validation, aes(x = factor(Threshold), y = R2)) +
  geom_text(
    aes(label = paste(print.p)),
    vjust = -1.5,
    hjust = 0,
    angle = 45,
    cex = 4,
    parse = T
  )  +
  scale_y_continuous(limits = c(0, max(prs.result_validation$R2) * 1.25)) +
  xlab(expression(italic(P) - value ~ threshold ~ (italic(Pvalue)[T]))) +
  ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
  geom_bar(aes(fill = -log10(Pvalue)), stat = "identity") +
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-4,
    name = bquote(atop(-log[10] ~ model, italic(Pvalue) - value),)) +
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size =
                                  18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust =
                                 1))
ggsave(plot=p1 , paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_eval_tresholds_validation.pdf"))

to_plot <- fread(paste0(PATH_TO_THRESH,"EPI_all_2022_covar.txt"))
q <- ggplot(data = to_plot, aes(x = AGE, group=as.factor(PHENO), fill=as.factor(PHENO))) +
  geom_density(alpha = 0.3, adjust=3) + 
  labs(fill="Phenotype")
 

ggsave(plot=q , paste0(PATH_TO_THRESH,"EPI_all_2022_PRS_AGE_cohorts.pdf"))


}

tresh <- arg[4]
print(tresh)


if ( is.na(tresh)){
  stop("Indicate best Treshold to save data (-t).\n")
}
if (! is.na(tresh) ){

  # Read in the fam file 
  phenotype <- read.table(FAM_FILE, header=FALSE)[,c(1,2,6)]
  colnames(phenotype) <- c("FID", "IID","PHENO")
  # Read in the PCs
  pcs <- read.table(COVAR_FILE, header=T)[,c(1,1,5:8)] #borrar para usar todas las PC
  # The default output from plink does not include a header
  # To make things simple, we will add the appropriate headers
  colnames(pcs) <- c("FID", "IID", paste0("PC",1:4))  ###cambair a 8 para todas las PC
  # Read in the covariates (here, it is age)
  covariate <- fread(COVAR_FILE, header=T)[,c(1,3,2,4)] #añadir 6 para TD
  colnames(covariate) <- c("FID", "AGE", "SEX","TD") #añadir TD
  covariate$IID <- covariate$FID
  # read covariates PD EPI
  covariate_diagnosis <- fread(COVAR_FILE, header=T)[,c(1,9)]
  colnames(covariate_diagnosis) <- c("FID", "PD")
  covariate_diagnosis$IID <- covariate_diagnosis$FID
  
  
  # Now merge the files
  pheno_2 <- merge(merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID")),covariate_diagnosis,by=c("FID","IID"))
  pheno_1 <- pheno_2[pheno_2$PHENO != -9,]
  pheno <- pheno_1
  #recodificar casos y controles
  pheno$PHENO <- pheno$PHENO - 1   
  PRS_TRESH <- arg[4] 
  print(PRS_TRESH)
  prs <- read.table(paste0(PATH_TO_CLUMP,"SS_PLINK1.9_PRS.",PRS_TRESH,".profile"), header=T)[,c(2,2,4,5,6)]
  colnames(prs)<- c("FID","IID","Total","N","SCORE")
  
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns
  pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
  #setting SCORE to Zscore
  meanControls <- mean(pheno.prs$SCORE[pheno.prs$PHENO == 0])
  sdControls <- sd(pheno.prs$SCORE[pheno.prs$PHENO == 0])
  pheno.prs$zSCORE <- (pheno.prs$SCORE - meanControls)/sdControls 
  fwrite(pheno.prs ,paste0(cmrb,"/",fam,"/",cmrb,"_",PRS_TRESH,"_",soft,"_PRS.covar.txt"),sep='\t')
  print("Data_saved")
}

