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
