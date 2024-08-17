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


