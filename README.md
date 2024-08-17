# PGS-AD-comorbidities
Pipeline with methodology to replicate manuscript "Polygenic Score Analysis Reveals Distinct Genetic Risk Profiles in the Presentation of Alzheimer's Disease Comorbidities"
This tutorial allow to replicate the analysis of PGS of Epilepsy in Alzheimer disease Patients from UK biobank. 
## Files
0. *EPI_GWAS* (NOT INCLUDED): Summary statistics from the last meatanalysis of epilepsy .
1. *over chain*: chain file to liftover GWAS SS.
2. COVARS.tsv*: Simulated Covariables to use in PRS calcualtion
3. *genetic_sex.tsv*: Simulated genetic sex to use in PRS calcualtion
4. *PCA.tsv*: Simulated PCA to use in PRS calcualtion
5. *releated_ids.tsv*: Simulated related individuals to use in PRS calcualtion
6. *CMRBs_PRS_pheno_covar.tsv.gz*: Results for each comorbidity to replicate Manuscript figures to use in PRS calcualtion
## Dependencies
1. [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
2. [PLINK 2](https://www.cog-genomics.org/plink/2.0/)
3. [Liftover](https://genome.sph.umich.edu/wiki/LiftOver)
4. [R](https://www.r-project.org/)

## Directory structure
	PGS-PD-Comorbidities
	├── PRS_EPI_all_2022
	│   ├── 1_SS_GWAS_EPI_all_2022
	|   ├── 2_QC_target_data     
	│   ├── 3_Clump_PRS
	│   ├── 4_Tresholding
 	│   └── Metadata
  	|       ├── COVARS.tsv
  	|       ├── genetic_sex.tsv
  	|       ├── PCA.tsv
	│       └── releated_ids.tsv
 	├── Figures
	│   ├── T2D_PRS_pheno_covar.tsv.gz
 	│   ├── MDD_PRS_pheno_covar.tsv.gz
	│   ├── MH_PRS_pheno_covar.tsv.gz
	│   ├── EPI_PRS_pheno_covar.tsv.gz
	│   ├── Height_PRS_pheno_covar.tsv.gz
	│   └── Figure*.R 
	├── 1_GWAS_liftover.sh
	├── 2_GWAS_QC.sh
 	├── 3_QC_UKBB.r
	├── 4_PCA_only_ukb.r
	├── 5_REL_IDS.r
	└── 6_CMRB_Best_PRS.R
## Step 0: Obtain Datasets
To generate the RPS we need base data (GWAS SS) and target data (genotypes). The genotypes used correspond to the data imputed from the UK biobank. Due to the availability of the dataset we will use sample data from 1000 genomes.
### 1. Download WGS data for chr 22
```
PATH_TO_QC=PRS_EPI_all_2022/2_QC_target_data
wget -P $PATH_TO_QC -c	http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz 
plink2 --vcf ${PATH_TO_QC}/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --export bgen-1.2 --out ${$PATH_TO_QC}/1000Genomes_chr22 
```
## Step 1: Quality controls for GWAS summary statistics
For a conplete detail of Quality controls (QC) procedures refers to [Choi *et al*, 2020](https://www.nature.com/articles/s41596-020-0353-1). Summary statistics from Epilepsy can be downloaded from [here](https://www.epigad.org/download/final_sumstats.zip), or with the code below. GWAS of epilepsy can be found in [ILAE *et al*, 2022](https://www.nature.com/articles/s41588-023-01485-w#article-info).  

### 1. Downlaoad GWAS EPI 
```
GWAS_SS_PATH=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022
mkdir -p $GWAS_SS_PATH
wget -P $GWAS_SS_PATH -c https://www.epigad.org/download/final_sumstats.zip
unzip -d $GWAS_SS_PATH $GWAS_SS_PATH/final_sumstats.zip 

```
### 2 .Adjust GWAS SS
GWAS SS file need to be formated in this way
A Space delimieted file with a header like this "CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value"
1. CHR
2. BP
3. MarkerName (If not present, fill with NA)
4. Allele1 (Effect allele gwas)
5. Allele2 (Non-effect allele)
6. Freq1 (If not present, fill with NA)
7. Beta
8. P-value

In order to run this tutorial with other comorbidities, please follow the recommendations above, this is only an example of processing for epilepsy GWAS.
```
GWAS_SS_PATH=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022
{
echo "CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value"
cat $GWAS_SS_PATH/ILAE3_TRANS_all_epilepsy_final.tbl | sed '1d' | grep '^22' | sed 's/\t/ /g' | awk -F' ' '{print $1,$2,$3,$4,$5,$6,$12,$10}'  
}  > $GWAS_SS_PATH/GWAS_SS_EPI_all_2022_RAW_hg19.txt
```
### 3. Liftover GWAS SS to GRCh38
Most of the GWAS are published in GRCh37 version. And the imputed data from UK Biobank is in GRCh38, so it is necessary to update the GWAS SS potions to GRCh48 version with Liftover. The following bash script can be used to generate the GWAS SS in GRCh38.
```
#!/usr/bin/bash
set -e
#### ARGUMENTS #####
# -c Comorbidity GWAS used
# -d Directory
# -X Indicate if GWAS SS has X chromosome
# -C Path to over cahin file to liftover


while getopts "c:d:X:C:" flag ;
do
    case "${flag}" in
        c) CMRB=${OPTARG};;
        d) WORK_START=${OPTARG};;
        X) xchrms=${OPTARG};;
        C) CHAIN=${OPTARG};;
    esac
done
if [[ -z $WORK_START  ]] ;then WORK_START=$(pwd) ; fi
test() { 
    cd $WORK_START
}
test || {
echo "Enter a valid Directory" 1>&2
false
}
if [[ -z $CMRB ]]; then echo "Indicate some parameter for Comorbilities (-c)" ; exit 0 ; fi
test() { 
    cd PRS_$CMRB
}
test || {
echo "Enter a valid Comorbidity" 1>&2
false
}
if [[ -z $xchrms ]] ;then echo "-X not present: indicate if GWAS SS have X chromosome (true or false)"; exit 0 ; fi
if [[ $xchrms == true ]] ;then 
    CHRS=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X   ) ;### variable util para loops
elif [[ $xchrms == false ]] ; then
    CHRS=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22  ) ;### variable util para loops
else 
    echo "invalid -X valid options are: true or false"
fi

if [[ -z $CHAIN ]] ;then echo "-C not present: indicate path to hg19ToHg38.over.chain.gz "; exit 0 ; fi


WORKD=$(pwd)
echo Working in $WORKD 

### generate .bed file to liftover
cd $WORKD/1_SS_GWAS_${CMRB}

printf "Formating GWAS SS to match hg19 to liftover: "
cat GWAS_SS_${CMRB}_RAW_hg19.txt | awk -F' ' '{if (NR>1) print "chr"$1 "\t" ($2-1) "\t" $2 "\t" $3","$4","$5","$6","$7","$8 }' | sed 's/chr23/chrX/g'  > GWAS_SS_hg19.bed
printf "DONE\n"

printf "Liftover to hg38: "
liftOver GWAS_SS_hg19.bed $CHAIN GWAS_SS_hg38.bed GWAS_SS_unlifted.bed
printf "DONE\n"

printf "Formating hg38 to GWAS SS pipeline PRS: "
cat GWAS_SS_hg38.bed | sed 's/,/\t/g' | sed 's/\t/ /g' | awk -F' ' '{printf $1; for (i=3; i <= NF; i++) printf FS$i; print NL }' | sed 's/^chr//g' > GWAS_SS_${CMRB}_RAW_hg38.txt 
printf "DONE\n"
```
```
./1_GWAS_liftover.sh -c EPI_all_2022 -X true -C hg19ToHg38.over.chain.gz
```
### 4. Perform QC steps in GWAS SS 
The GWAS SS quality controls consist of eliminating duplicate and ambiguous SNPs. The following bash script allows to generate the final GWAS SS.
```
#!/usr/bin/bash
set -e
#### ARGUMENTS #####
# -c Comorbidity GWAS used
# -d Directory
# -X Indicate if GWAS SS has X chromosome

while getopts "c:d:X:" flag ;
do
    case "${flag}" in
        c) CMRB=${OPTARG};;
        d) WORK_START=${OPTARG};;
        X) xchrms=${OPTARG};;

    esac
done
if [[ -z $WORK_START  ]] ;then WORK_START=$(pwd) ; fi
test() { 
    cd $WORK_START
}
test || {
echo "Enter a valid Directory" 1>&2
false
}
if [[ -z $CMRB ]]; then echo "Indicate some parameter for Comorbilities (-c)" ; exit 0 ; fi
test() { 
    cd PRS_$CMRB
}
test || {
echo "Enter a valid Comorbidity" 1>&2
false
}
if [[ -z $xchrms  ]] ;then echo "-X not present: indicate if GWAS SS have X chromosome (true or false)"; exit 0 ; fi


WORKD=$(pwd)
echo Working in $WORKD 

################ set Variables ##############################
if [[ $xchrms == true ]] ;then 
    CHRS=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X   ) ;### variable util para loops
elif [[ $xchrms == false ]] ; then
    CHRS=$(echo 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22  ) ;### variable util para loops
else 
    echo "invalid -X valid options are: true or false"
fi


########################################## QC GWAS SS
echo -e "\nQC base data, SS_GWAS_$CMRB\n"
cd $WORKD/1_SS_GWAS_$CMRB

############################ QC SS Comorbidities
##### Format
#header CHR BP MarkerName Allele1 Allele2 Freq1 Beta P-value
# 1 CHR
# 2 BP
# 3 MarkerName (If not present, fill with NA)
# 4 Allele1 (effect allele)
# 5 Allele2 (non-effect allele)
# 6 Freq1 (If not present, fill with NA)
# 7 Beta
# 8 P-value

cat GWAS_SS_${CMRB}_RAW_hg38.txt |  sort  -nk1,1 -k2,2 > ${CMRB}_sorted.txt
# filter
cat ${CMRB}_sorted.txt | awk '{if ($6 >= 0.01) {print $0}}' > ${CMRB}_filtered.txt

#Alleles to upercase
cat ${CMRB}_filtered.txt | awk '{print $1,$2,$3,toupper($4) ,toupper($5) ,$6 ,$7 ,$8 }' > ${CMRB}_filtered_touper.txt

#recode SNPS
cat ${CMRB}_filtered_touper.txt | awk -F' ' '{print " " $1":"$2":"$4":"$5 FS $3 FS $1 FS $2 FS $4 FS $5 FS $6 FS $7 FS $8}' > ${CMRB}_recode.txt
## remove duplciated
cat ${CMRB}_recode.txt | awk '{seen[$1$2]++; if(seen[$1$2]==1){ print}}' > ${CMRB}_nodup.txt

# Remove ambiguous SNPs
awk '!( ($5=="A" && $6=="T") || \
    ($5=="T" && $6=="A") || \
    ($5=="G" && $6=="C") || \
    ($5=="C" && $6=="G")) {print}' ${CMRB}_nodup.txt > ${CMRB}_no_ambiguos.txt

{
echo " ID MarkerName CHR BP Allele1 Allele2 Freq1 Beta P-value"
cat ${CMRB}_no_ambiguos.txt | sed 's/ 23:/X:/g' | sed 's/ 23 / X /g' 
} > final_${CMRB}_QC.txt 

```
```
./2_GWAS_QC.sh -c EPI_all_2022 -X true
```

## Step 2: Perform QC in UKBB data 
The reference and alternative alleles of the genotype must match those reported in the GWAS SS, so we must generate a quality control of the UK biobank in order to work with the SNPs in a unified way.
### 1. Filter BGEN files of UKB by allele frequency
Due to the magnitude of the imputed UKBB data we first filtered by allelic frequency greater than 1% and then performed the corresponding quality controls on smaller files.
The next step was performed in the UK Biobank Research Analysis Platform (RAP), due to the availability of the data, we provide example data from 1000 genomes
```
PATH_TO_QC=PRS_EPI_all_2022/2_QC_target_data
BGEN_PREFIX=${PATH_TO_QC}/1000Genomes_chr22

mkdir -p $PATH_TO_QC
cat ${BGEN_PREFIX}.sample | sed 's/0 0 0 0/0 0 0 D/g'> ${PATH_TO_QC}/sample_chr22.sample
plink2  --bgen ${BGEN_PREFIX}.bgen ref-first \
		--sample ${PATH_TO_QC}/sample_chr22.sample \
		--oxford-single-chr 22  \
		--set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 missing \
		--maf 0.01 \
		--make-pfile --out ${PATH_TO_QC}/1000Genomes_common_chr22
grep -v '#' ${PATH_TO_QC}/1000Genomes_common_chr22.pvar > ${PATH_TO_QC}/1000Genomes_common_chr22_noHeader.pvar
```
### 2. QC UK biobank 
Once the BGEN files have been filtered, it is necessary to identify the common SNPs between the UK bibank and the GWAS SS. To do this, using the .pvar file (provided in this tutorial) and the GWAS QC file generated in step 1, we identify the common SNPs and generate a list to update the .pvar file with the alleles reported in the GWAS SS. The following R script is the one in charge of what is described above
```
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
```
```
PATH_TO_QC=PRS_EPI_all_2022/2_QC_target_data
GWAS=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022/final_EPI_all_2022_QC.txt
PVAR=${PATH_TO_QC}/1000Genomes_common_chr22_noHeader.pvar
./3_QC_UKBB.r EPI_all_2022 ${GWAS} ${PVAR}
```

### 3. Extract common positions between UKBB and Epilepsy GWAS 
With the files generated in the previous step, we proceed to filter the UKbiobank imputed data and update the reported alleles. 
```
CMRB=EPI_all_2022
PATH_TO_QC=PRS_EPI_all_2022/2_QC_target_data
PFILE_PREFIX=${PATH_TO_QC}/1000Genomes_common_chr22

plink2  --pfile ${PFILE_PREFIX} --exclude ${PATH_TO_QC}/${CMRB}_chr22_snplist_pvar.mismatch --ref-allele force ${PATH_TO_QC}/${CMRB}_chr22_pvar_updated --make-pfile --out ${PFILE_PREFIX}_swapped
plink2  --pfile ${PFILE_PREFIX}_swapped --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 missing --make-bed --out ${PFILE_PREFIX}_recoded
```
## Step 3: PRS Calcualtions Clumping
With the common SNPs we proceed to clumping, which consists of selecting independent variants in the haplotype blocks, considering the GWAS SS pvalue. We will use the clumped data to calculate the PRS using different pvalue thresholds.  
### 1. Clumping
Consists of taking the independent alleles taking into consideration the GWAS Pvalue.
```
GWAS=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022/final_EPI_all_2022_QC.txt
PATH_TO_CLUMP=PRS_EPI_all_2022/3_Clump_PRS
BFILE_PREFIX=${PATH_TO_QC}/1000Genomes_common_chr22_recoded

mkdir -p $PATH_TO_CLUMP
plink   --bfile ${BFILE_PREFIX} --clump-p1 1 --clump-r2 0.1 --clump-kb 500 --clump ${GWAS} --clump-snp-field ID --clump-field P-value --out ${PATH_TO_CLUMP}/clumped_SS_chr22
awk 'NR!=1{print $3}' ${PATH_TO_CLUMP}/clumped_SS_chr22.clumped | sed '/^$/d' >  ${PATH_TO_CLUMP}/SS_chr22.valid.snp
plink2  --bfile ${BFILE_PREFIX} --extract ${PATH_TO_CLUMP}/SS_chr22.valid.snp --make-bed --out ${PATH_TO_CLUMP}/clumped_SS_chr22
```
#### 2. Merge clumped data by chromosome in a unique file
When we work with all the data we must unify the chromosomes in a single file using the dociog below, for this case we will only use chromosome 22, omit this execution.
```
#ls | grep ".pvar" | grep "clumped" > TEMP1.txt
#awk -F'.' '{print $1 }' TEMP1.txt > merge_list.chrALL.txt
#rm TEMP1.txt
#plink2 --pmerge-list merge_list.chrALL.txt --make-bed --out clumped_SS_chr.all 
```
### 3. Define SNPs and their Pvalue 
It takes the GWAS SS and extracts the SNPs ID along with their pvalue, to be used to define the snps that meet the pvalue treshold. additionally it generates the list of the pvalues included in each treshold AND filters the GWAS SS with the valid SNPs obtained in step 2.
```
GWAS=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022/final_EPI_all_2022_QC.txt
PATH_TO_CLUMP=PRS_EPI_all_2022/3_Clump_PRS

#Generates SNP and Pvalue in one file
cat $GWAS | awk -F' ' '{print $1,$9}' | grep -wFf ${PATH_TO_CLUMP}/SS_chr22.valid.snp  > ${PATH_TO_CLUMP}/SNP.pvalue

#define threshold to use in Thresholding step
echo "0.00001 0 0.00001" > ${PATH_TO_CLUMP}/range_list
echo "0.0001 0 0.0001" >> ${PATH_TO_CLUMP}/range_list
echo "0.001 0 0.001" >> ${PATH_TO_CLUMP}/range_list
echo "0.01 0 0.01" >> ${PATH_TO_CLUMP}/range_list
echo "0.05 0 0.05" >> ${PATH_TO_CLUMP}/range_list
echo "0.1 0 0.1" >> ${PATH_TO_CLUMP}/range_list
echo "0.2 0 0.2" >> ${PATH_TO_CLUMP}/range_list
echo "0.3 0 0.3" >> ${PATH_TO_CLUMP}/range_list
echo "0.4 0 0.4" >> ${PATH_TO_CLUMP}/range_list
echo "0.5 0 0.5" >> ${PATH_TO_CLUMP}/range_list
echo "0.8 0 0.8" >> ${PATH_TO_CLUMP}/range_list
echo "1 0 1" >> ${PATH_TO_CLUMP}/range_list 

#Filter GWAS SS with valid SNPs
grep -wFf ${PATH_TO_CLUMP}/SS_chr22.valid.snp $GWAS > ${PATH_TO_CLUMP}/clumped_final_SS_common_SNPS_all.txt

```
### 5. Calculate PRS for all Pvalue trehshold
Takes the files generated above to calculate the PRS for each of the pvalue treshold defined above
```
GWAS=PRS_EPI_all_2022/1_SS_GWAS_EPI_all_2022/final_EPI_all_2022_QC.txt
PATH_TO_CLUMP=PRS_EPI_all_2022/3_Clump_PRS
CLUMPED_PREFIX=${PATH_TO_CLUMP}/clumped_SS_chr22

plink --bfile ${CLUMPED_PREFIX} --score ${PATH_TO_CLUMP}/clumped_final_SS_common_SNPS_all.txt 1 5 8 header --q-score-range ${PATH_TO_CLUMP}/range_list ${PATH_TO_CLUMP}/SNP.pvalue --out ${PATH_TO_CLUMP}/SS_PLINK1.9_PRS
```
## Step 4: PRS threshholding
Thresholding consists of taking the PRS calculated for the different tresholds and selecting the one that best explains the variability of the data.
### 1. Sample QC
The quality controls of the individuals were performed by extracting the available information from the cohort browser of dnanexus, for this example we will make use of synthesized clinical data to perform the tresholding and sample QC steps.
```
PATH_TO_CLUMP=PRS_EPI_all_2022/3_Clump_PRS/
PATH_TO_THRESH=PRS_EPI_all_2022/4_Thresholding
METADATA=PRS_EPI_all_2022/Metadata/
mkdir -p $PATH_TO_THRESH
mkdir -p $PATH_TO_THRESH/Sample_QC
PATH_TO_SAMPLEQC=$PATH_TO_THRESH/Sample_QC
```
#### a. Ancestry QC
```
cat $PATH_TO_CLUMP/clumped_SS_chr22.psam | grep -v '#'| awk -F' ' '{print $1}' | sort  > ${PATH_TO_SAMPLEQC}/Phenotipo.Ids     
mkdir -p ${PATH_TO_SAMPLEQC}/PCA_results     
./4_PCA_only_ukb.r EPI_all_2022 ${METADATA}/PCA.tsv
grep -wf ${PATH_TO_SAMPLEQC}/keep_ancestry_eur.txt ${PATH_TO_SAMPLEQC}/Phenotipo.Ids > ${PATH_TO_SAMPLEQC}/Phenotipo_QC_EUR.Ids 
```
#### b. Discordant Genetic sex QC
```
cut -f2 ${METADATA}/COVARS.tsv | paste -d'\t' - ${METADATA}/genetic_sex.tsv > ${PATH_TO_SAMPLEQC}/TEMP  
cat ${PATH_TO_SAMPLEQC}/TEMP | awk -F'\t' '{if ($1 != $3) print  }' | cut -f2 | grep -v "ID" > ${PATH_TO_SAMPLEQC}/remove_sex_discordant.ids 
```
#### c. Releted individual QC
```
cat  ${METADATA}/releated_ids.tsv | awk -F'\t' '{if ($5 > 0.0442 ) print $1 FS $2}' | grep -v "ID"   >   ${PATH_TO_SAMPLEQC}/rel_IDs_all_filter.txt
cut -f1 ${PATH_TO_SAMPLEQC}/rel_IDs_all_filter.txt > ${PATH_TO_SAMPLEQC}/rel_ID1.txt 
cut -f2 ${PATH_TO_SAMPLEQC}/rel_IDs_all_filter.txt > ${PATH_TO_SAMPLEQC}/rel_ID2.txt 
./5_REL_IDS.r ${PATH_TO_SAMPLEQC}/rel_ID1.txt ${PATH_TO_SAMPLEQC}/rel_ID2.txt ${METADATA}/COVARS.tsv 
cat ${PATH_TO_SAMPLEQC}/REL_IDS_to_remove.txt | sort | uniq | grep -v "remove" > ${PATH_TO_SAMPLEQC}/TEMP 
sed -e "s/\r//" ${PATH_TO_SAMPLEQC}/TEMP  > ${PATH_TO_SAMPLEQC}/remove_REL.ids 
```
#### d. Generate final QCded Ids
```
cat ${PATH_TO_SAMPLEQC}/Phenotipo_QC_EUR.Ids | grep -vwFf ${PATH_TO_SAMPLEQC}/remove_REL.ids  | grep -vwFf ${PATH_TO_SAMPLEQC}/remove_sex_discordant.ids   > ${PATH_TO_THRESH}/Phenotipo_QC_EUR_final.Ids 
```
### 2. Correct .PSAM file and generate .FAM to work in PLINK1.9
```
PATH_TO_CLUMP=PRS_EPI_all_2022/3_Clump_PRS/
PATH_TO_THRESH=PRS_EPI_all_2022/4_Thresholding/
METADATA=PRS_EPI_all_2022/Metadata/

awk -F'\t' '{print $1}' ${PATH_TO_THRESH}/Phenotipo_QC_EUR_final.Ids > ${PATH_TO_THRESH}/keep_Ids 
cat ${METADATA}/COVARS.tsv | cut -f1,10 |  awk -F'\t' 'BEGIN {print "#IID\tPHENO1"} {if(NR > 1 ) print $1 FS $2 + 1}' > ${PATH_TO_THRESH}/case_controls.txt 
cat ${METADATA}/COVARS.tsv | cut -f1,2 | awk -F'\t' 'BEGIN {print "#IID\tSEX"} {if(NR > 1 ) print $1 FS $2 +1 }' > ${PATH_TO_THRESH}/update_sex.txt 

plink2 --psam ${PATH_TO_CLUMP}/clumped_SS_chr22.psam \
	--pheno ${PATH_TO_THRESH}/case_controls.txt  \
	--keep ${PATH_TO_THRESH}/keep_Ids \
	--update-sex ${PATH_TO_THRESH}/update_sex.txt \
	--make-just-psam --out ${PATH_TO_THRESH}/final.Phenotipo_QC

cat ${PATH_TO_THRESH}/final.Phenotipo_QC.psam | awk -F'\t' '{if (NR != 1) print $1 "\t" $1 "\t0\t0\t" $2 "\t" $3  }'  > ${PATH_TO_THRESH}/final.Phenotipo_QC.fam    
```
### 3. Generate a list of pvalue threshold obtaind (Sometimes not all Pvalues are presetn in GWAS SS)
```
PATH_TO_CLUMP=PRS_EPI_all_2022/3_Clump_PRS/
PATH_TO_THRESH=PRS_EPI_all_2022/4_Thresholding/
SUFFIX_PRS=SS_PLINK1.9_PRS.

ls ${PATH_TO_CLUMP}${SUFFIX_PRS}* | grep 'profile' | sed -e 's/PRS_EPI_all_2022\/3_Clump_PRS\/SS_PLINK1.9_PRS.//g' | sed 's/.profile//g' > ${PATH_TO_THRESH}/used_tresh.txt 
```
### 4. Identify which Pvalue Threshold best describe the data
```
PATH_TO_CLUMP=PRS_EPI_all_2022/3_Clump_PRS/
PATH_TO_THRESH=PRS_EPI_all_2022/4_Thresholding/
METADATA=PRS_EPI_all_2022/Metadata/

./6_CMRB_Best_PRS.R ${PATH_TO_THRESH}/used_tresh.txt ${PATH_TO_THRESH}/final.Phenotipo_QC.fam ${METADATA}/COVARS.tsv
```
### 4. Save PRS at Pvalue Threshold obtained (Indicate the pvalue obtained in the last step, in this case 0.01)
```
PATH_TO_CLUMP=PRS_EPI_all_2022/3_Clump_PRS/
PATH_TO_THRESH=PRS_EPI_all_2022/4_Thresholding/
METADATA=PRS_EPI_all_2022/Metadata/

./6_CMRB_Best_PRS.R ${PATH_TO_THRESH}/used_tresh.txt ${PATH_TO_THRESH}/final.Phenotipo_QC.fam ${METADATA}/COVARS.tsv 0.01
```
## Step 5. Replicate Figures
The following R scripts allow you to produce the results of the manuscript. I make available the final result obtained by changing the IDS coding.
```
cd Figures
gunzip *_PRS_pheno_covar.tsv.gz

./Figure_2.R
./Figure_3.R
./Figure_4.R
./Figure_5.R
./Figure_Sup_1.R
./Figure_Sup_2.R
./Figure_Sup_3.R
./Figure_Sup_4.R
./Figure_Sup_5.R
./Figure_Sup_6.R
./Figure_Sup_7.R
./Figure_Sup_8.R
```




