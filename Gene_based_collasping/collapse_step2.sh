#!/bin/bash
module load miniconda3/RSAIGE
Pheno=$1
for i in {1..22}; do 
Rscript step2_SPAtests.R  \
--bedFile=/public/home2/nodecw_group/UKB_WES_data/qcstep5/unrelated_0_0884/ukb_wes_chr${i}_sample_qc_final_unrelated.bed  \
--bimFile=/public/home2/nodecw_group/UKB_WES_data/qcstep5/unrelated_0_0884/ukb_wes_chr${i}_sample_qc_final_unrelated.bim  \
--famFile=/public/home2/nodecw_group/UKB_WES_data/qcstep5/unrelated_0_0884/ukb_wes_chr${i}_sample_qc_final_unrelated.fam  \
--sparseGRMFile=/public/home2/nodecw_group/UKB_WES_data/SAIGE/Huashan1/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx  \
--sparseGRMSampleIDFile=/public/home2/nodecw_group/UKB_WES_data/SAIGE/Huashan1/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  \
--groupFile=/public/home2/nodecw_group/UKB_WES_data/snpEff/output/SAIGE_group_file/lof_missense_five/SnpEff_gene_group_chr${i}.txt \
--annotation_in_groupTest="lof,missense:lof,missense"  \
--AlleleOrder=alt-first  \
--minMAF=0  \
--minMAC=0.5  \
--maxMAF_in_groupTest=0.00001,0.0001,0.001,0.01  \
--is_output_markerList_in_groupTest=TRUE  \
--LOCO=FALSE  \
--is_fastTest=TRUE  \
--SAIGEOutputFile=${Pheno}_chr${i}_10PC_both.txt  \
--GMMATmodelFile=${Pheno}_s1_chr${i}_10PC_both.rda  \
--varianceRatioFile=${Pheno}_s1_chr${i}_10PC_both.varianceRatio.txt 
done
