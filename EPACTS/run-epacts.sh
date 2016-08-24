#!/bin/bash
module load vcftools
module load EPACTS
module load R

echo ------------------------------------ STEP1 CREATE PLINK FILES ------------------------------------
gunzip all_filtered.multianno.recode.vcf.gz
vcftools --vcf all_filtered.multianno.recode.vcf --plink --out all_filtered.multianno.recode

echo ------------------------------------ STEP2 MODIFY THE PED FILE ------------------------------------
python modify_ped.py all_filtered.multianno.recode.ped

echo ------------------------------------ STEP3 PREPARE INPUT FILE FOR EPACTS ------------------------------------
bgzip all_filtered.multianno.recode.vcf
tabix -pvcf -f all_filtered.multianno.recode.vcf.gz

echo ------------------------------------ STEP4 MAKE GROUP FOR VARIANTS ------------------------------------
./make_group.py -i all_filtered.multianno.recode.vcf.gz -o all_filtered.multianno.recode.grp

echo ------------------------------------ STEP5 MAKE KINSHIP MATRIX ------------------------------------
#epacts make-kin --vcf all_filtered.multianno.recode.vcf.gz --ped all_filtered_with_head.ped --min-maf 0.01 --min-CallRate 0.95 --out all_filtered.multianno.recode.kinf
epacts make-kin --vcf all_filtered.multianno.recode.vcf.gz --ped all_filtered_with_head.ped --min-CallRate 0.95 --out all_filtered.multianno.recode.kinf

make -f all_filtered.multianno.recode.kinf.Makefile -j 8

#STEP6 PERFORM TEST 
#epacts group --groupf all_filtered.multianno.recode.grp --vcf all_filtered.multianno.recode.vcf.gz --max-maf 0.05 --ped all_filtered_with_head.ped \
#	--kin all_filtered.multianno.recode.kinf --test emmaxVT --out EUR-emmaxV


#epacts group --groupf all_filtered.multianno.recode.grp --vcf all_filtered.multianno.recode.vcf.gz --ped all_filtered_with_head.ped \
#	--kin all_filtered.multianno.recode.kinf --test emmaxVT --out emmaxVT

#epacts group --groupf all_filtered.multianno.recode.grp --vcf all_filtered.multianno.recode.vcf.gz --ped all_filtered_with_head.ped \
#	--kin all_filtered.multianno.recode.kinf --test skat --out skato --skat-o --skat-adjust

#python select_by_gene.py -i all_filtered.multianno.recode.vcf.gz -o selected.txt -p all_filtered_with_head.ped -g OXR1 RPE65 GPR1 PDE6B
