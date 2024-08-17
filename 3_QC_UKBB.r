#!/usr/bin/env Rscript
####USAGE######
# ./3_QC_UKB.r [Comorbidity] {final GWAS SS} {.pvar file no header}
#
#
#
#
args <- commandArgs(trailingOnly = TRUE)

library(data.table)
library(magrittr)
library(stringr)


# Function for finding the complementary allele
complement <- function(x) {
    switch (
        x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}
CMRB=args[1]
GWAS_SS=args[2]
PVAR_FILE=args[3]
# Read in summary statistic data (require data.table v1.12.0+)
print("Reading SS file")
SS <- fread(GWAS_SS,colClasses=c(CHR="character")) %>%
    setnames(., colnames(.), c("ID", "Marker", "CHR", "BP", "A1", "A2","Freq1","Beta","P-value")) %>%
    # And immediately change the alleles to upper cases
    .[,c("A1","A2"):=list(toupper(A1), toupper(A2))] 
#str(height)
for (chr in 22) {
  # Read in bim file  A1=minor A2=major
    print(paste("Reading BIM file of chromosome",chr,sep=" "))
    bim <- fread(PVAR_FILE,colClasses=c(V1="character")) %>%
        # Note: . represents the output from previous step
        # The syntax here means, setnames of the data read from
        # the bim file, and replace the original column names by 
        # the new names
       setnames(., colnames(.), c("CHR", "BP", "ID", "A1", "A2")) %>%
        # And immediately change the alleles to upper cases
        .[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
    bim$Position <- str_c(bim$CHR,bim$BP,sep=":")
    colnames(bim) <- c("CHR", "BP", "ID", "B.A1", "B.A2","Position") 
    
    #filter SS to merge
    SS_filtered <- SS[SS$CHR == chr]
    SS_filtered$Position <- str_c(SS_filtered$CHR,SS_filtered$BP,sep=":") 
    # Merge summary statistic with target
    print("Processing")
    info <- merge(bim, SS_filtered, 
              by = c("CHR","BP","Position"),
              suffixes = c("_bim","_SS") )
    
    # Get SNPs that have the same alleles across base and target
    info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
    
    # Identify SNPs that are complementary between base and target
    info$C.A1 <- sapply(info$B.A1, complement)
    info$C.A2 <- sapply(info$B.A2, complement)
    info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
    # Update the complementary alleles in the bim file
    # This allow us to match the allele in subsequent analysis
    complement.snps <- bim$ID %in% info.complement$ID_bim
    bim[complement.snps,]$B.A1 <-
        sapply(bim[complement.snps,]$B.A1, complement)
    bim[complement.snps,]$B.A2 <-
        sapply(bim[complement.snps,]$B.A2, complement)
    
    # identify SNPs that need recoding
    info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
    # Update the recode SNPs
    recode.snps <- bim$ID %in% info.recode$ID_bim
    tmp <- bim[recode.snps,]$B.A1
    bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
    bim[recode.snps,]$B.A2 <- tmp
    
    # identify SNPs that need recoding & complement
    info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
    # Update the recode + strand flip SNPs
    com.snps <- bim$ID %in% info.crecode$ID_bim
    tmp <- bim[com.snps,]$B.A1
    bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
    bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
    
    # Output updated bim file
    print(paste("Writing pvar File updated of chromosome",chr,sep=" "))
    fwrite(bim[,c("ID", "B.A1")], paste("PRS_",CMRB,"/2_QC_target_data/",CMRB,"_chr",chr,"_pvar_updated",sep=""),quote = F,row.names = F , col.names=F, sep="\t")
    
    # Output mismatch file
    print(paste("Writing QCded SNPs of chromosome",chr,sep=" "))
    mismatch <-
        bim$ID[!(bim$ID %in% info.match$ID_bim |
                  bim$ID %in% info.complement$ID_bim | 
                  bim$ID %in% info.recode$ID_bim |
                  bim$ID %in% info.crecode$ID_bim)]
    write.table(mismatch, paste("PRS_",CMRB,"/2_QC_target_data/",CMRB,"_chr",chr,"_snplist_pvar.mismatch",sep=""), quote=F, row.names=F, col.names=F)
    
}
