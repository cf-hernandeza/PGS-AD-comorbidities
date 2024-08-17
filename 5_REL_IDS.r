#!/usr/bin/env Rscript
####USAGE######
# ./5_REL_IDS.r  {ID1} {ID2} {COVARS}
#
#
#
#
args <- commandArgs(trailingOnly = TRUE)

### load libraries
library(data.table)
library(dplyr)


ID1=args[1]
ID2=args[2]
COVAR_FILE=args[3]

REL1 <- fread(ID1, header=FALSE)  
colnames(REL1) <- c("ID")
REL1$nrow <- c(1:length(REL1$ID))


REL2 <- fread(ID2, header=FALSE)  
colnames(REL2) <- c("ID")
REL2$nrow <- c(1:length(REL2$ID))

covar <- fread(COVAR_FILE)
covar$SCORE <- covar$PD * 2 +  as.numeric(covar$EPI )

REL1 <- merge(REL1,covar,by="ID")
REL2 <- merge(REL2,covar,by="ID")

data <- merge(REL1,REL2,by="nrow")


to_remove <- ifelse(data$SCORE.x == data$SCORE.y, 
                    ifelse(data$age.x >= data$age.y, data$ID.y, data$ID.x),
                    ifelse(data$SCORE.x > data$SCORE.y,data$ID.y,data$ID.x)
                    )


fwrite (as.data.frame(to_remove), "PRS_EPI_all_2022/4_Thresholding/Sample_QC/REL_IDS_to_remove.txt")

