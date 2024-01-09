#!/bin/bash
module load miniconda3/RSAIGE
Pheno=$1
for i in {1..22}; do
Rscript step2_SPAtests.R  \
--bedFile=/mnt/storage/home1/Huashan1/UKB_WES_data/qcstep5/unrelated_0_0884/ukb_wes_chr${i}_sample_qc_final_unrelated.bed  \
--bimFile=/mnt/storage/home1/Huashan1/UKB_WES_data/qcstep5/unrelated_0_0884/ukb_wes_chr${i}_sample_qc_final_unrelated.bim  \
--famFile=/mnt/storage/home1/Huashan1/UKB_WES_data/qcstep5/unrelated_0_0884/ukb_wes_chr${i}_sample_qc_final_unrelated.fam  \
--sparseGRMFile=/mnt/storage/home1/Huashan1/UKB_WES_data/SAIGE/Huashan1/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx  \
--sparseGRMSampleIDFile=/mnt/storage/home1/Huashan1/UKB_WES_data/SAIGE/Huashan1/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  \
--AlleleOrder=alt-first \
--minMAF=0 \
--minMAC=20 \
--LOCO=FALSE \
--is_fastTest=TRUE \
--GMMATmodelFile=${Pheno}_s1_chr${i}_10PC_both.rda  \
--varianceRatioFile=${Pheno}_s1_chr${i}_10PC_both.varianceRatio.txt \
--pCutoffforFirth=0.05 \
--is_Firth_beta=TRUE \
--is_output_moreDetails=TRUE \
--SAIGEOutputFile=${Pheno}_chr${i}_single.txt 
done