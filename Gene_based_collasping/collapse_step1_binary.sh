#!/bin/bash
module load miniconda3/RSAIGE
Pheno=$1
for i in {1..22};do
Rscript step1_fitNULLGLMM.R   \
--sparseGRMFile=/mnt/storage/home1/Huashan1/UKB_WES_data/SAIGE/Huashan1/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
--sparseGRMSampleIDFile=/mnt/storage/home1/Huashan1/UKB_WES_data/SAIGE/Huashan1/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
--plinkFile=/mnt/storage/home1/Huashan1/UKB_WES_data/qcstep5/unrelated_0_0884/ukb_wes_chr${i}_sample_qc_final_unrelated \
--useSparseGRMtoFitNULL=FALSE   \
--useSparseGRMforVarRatio=TRUE \
--phenoFile=${Pheno}.csv \
--phenoCol=${Pheno} \
--covarColList=Sex,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--sampleIDColinphenoFile=eid \
--isCovariateOffset=FALSE \
--traitType=binary       \
--nThreads=30   \
--isCateVarianceRatio=FALSE \
--outputPrefix=${Pheno}_s1_chr${i}_10PC_both	\
--IsOverwriteVarianceRatioFile=TRUE \
--invNormalize=TRUE 
done
