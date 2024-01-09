#!/bin/bash
Pheno=$1
for i in {1..22};
do
plink2 \
--bfile /public/home2/nodecw_group/UKB_gene_v3_imp_qc/UKB_gene_v3_imp_qc_chr${i} \
--keep /public/home/geyijun/resources/UKB_genetic/british_id.txt \
--exclude /public/home2/nodecw_group/UKB_gene_v3_imp_qc_excludesnp/snp_chr${i}.txt \
--glm hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,beta,se,tz,p \
--covar-variance-standardize \
--geno 0.05 \
--mind 0.05 \
--maf 0.001 \
--hwe 1e-6 \
--1 \
--pheno ${Pheno}.txt \
--covar cov.txt \
--out ${Pheno}_chr${i}_snp_based.csv \
--threads 30 \
;done