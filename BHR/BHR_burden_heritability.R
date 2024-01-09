####WES BHR_burden_heretability_aggregate####
library(bhr)
library(tidyr)
library(readr)
library(data.table)

###
single_result_path="/BHR/BHR_input/"
###
setwd(single_result_path)
filelist <- list.files(single_result_path)
filelist<-grep("_plof.txt",filelist,value = T)
pheno_list<-c(0)
for (i in 1:length(filelist)) {
  pheno_path<-filelist[i]
  pheno_name<-gsub("_plof.txt","",pheno_path)
  
  ###lof
  BHR_input_path<-paste("/BHR/BHR_input/",pheno_name,"_plof.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  baseline <- read.table("/BHR/ms_baseline_oe5.txt")
  
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-5)
  pheno_name_bin1<-paste(pheno_name,"bin1_plof",sep="_")
  assign(pheno_name_bin1,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name_bin1)
  
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-4 & BHR_input$AF>=1e-5)
  pheno_name_bin2<-paste(pheno_name,"bin2_plof",sep="_")
  assign(pheno_name_bin2,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name_bin2)
  
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-3 & BHR_input$AF>=1e-4)
  pheno_name_bin3<-paste(pheno_name,"bin3_plof",sep="_")
  assign(pheno_name_bin3,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name_bin3)
  
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-2 & BHR_input$AF>=1e-3)
  pheno_name_bin4<-paste(pheno_name,"bin4_plof",sep="_")
  assign(pheno_name_bin4,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name_bin4)
  
  ###missense
  BHR_input_path<-paste("/BHR/BHR_input/",pheno_name,"_missense.txt",sep = "")
  BHR_input<-as.data.frame(fread(BHR_input_path))
  
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-5)
  pheno_name_bin1<-paste(pheno_name,"bin1_missense",sep="_")
  assign(pheno_name_bin1,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name_bin1)
  
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-4 & BHR_input$AF>=1e-5)
  pheno_name_bin2<-paste(pheno_name,"bin2_missense",sep="_")
  assign(pheno_name_bin2,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name_bin2)
  
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-3 & BHR_input$AF>=1e-4)
  pheno_name_bin3<-paste(pheno_name,"bin3_missense",sep="_")
  assign(pheno_name_bin3,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name_bin3)
  
  BHR_input_rare<-subset(BHR_input,BHR_input$AF<1e-2 & BHR_input$AF>=1e-3)
  pheno_name_bin4<-paste(pheno_name,"bin4_missense",sep="_")
  assign(pheno_name_bin4,BHR_input_rare)
  pheno_list<-c(pheno_list,pheno_name_bin4)
}
pheno_list<-pheno_list[-1]

result_aggregate<-data.frame(pheno_name=as.character(0),aggregate_h2=as.numeric(0),aggregate_h2_se=as.numeric(0),maf_func_bins=as.character(0))
for (i in 1:length(filelist)) {
  pheno_path<-filelist[i]
  pheno_name<-gsub("_plof.txt","",pheno_path)
  
  ##all
  BHR_aggregate<-BHR(mode = "aggregate", 
                     ss_list = list(get(paste(pheno_name,"bin1_plof",sep = "_")),get(paste(pheno_name,"bin2_plof",sep = "_")),get(paste(pheno_name,"bin3_plof",sep = "_")),get(paste(pheno_name,"bin4_plof",sep = "_")),
                                    get(paste(pheno_name,"bin1_missense",sep = "_")),get(paste(pheno_name,"bin2_missense",sep = "_")),get(paste(pheno_name,"bin3_missense",sep = "_")),get(paste(pheno_name,"bin4_missense",sep = "_"))),
                     trait_list = list(pheno_name),
                     annotations = list(baseline),
                     num_blocks = 100, 
                     genomewide_correction=TRUE,
                     overdispersion=TRUE)
  aggregate_h2<-BHR_aggregate$aggregated_mixed_model_h2
  aggregate_h2_se<-BHR_aggregate$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,aggregate_h2,aggregate_h2_se,"plof+missense")
  colnames(result)<-c("pheno_name","aggregate_h2","aggregate_h2_se","maf_func_bins")
  result_aggregate<-rbind(result_aggregate,result)
  
  ##lof
  BHR_aggregate<-BHR(mode = "aggregate", 
                     ss_list = list(get(paste(pheno_name,"bin1_plof",sep = "_")),get(paste(pheno_name,"bin2_plof",sep = "_")),get(paste(pheno_name,"bin3_plof",sep = "_")),get(paste(pheno_name,"bin4_plof",sep = "_"))),
                     trait_list = list(pheno_name),
                     annotations = list(baseline),
                     num_blocks = 100, 
                     genomewide_correction=TRUE,
                     overdispersion=TRUE)
  aggregate_h2<-BHR_aggregate$aggregated_mixed_model_h2
  aggregate_h2_se<-BHR_aggregate$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,aggregate_h2,aggregate_h2_se,"plof")
  colnames(result)<-c("pheno_name","aggregate_h2","aggregate_h2_se","maf_func_bins")
  result_aggregate<-rbind(result_aggregate,result)
  
  ##missense
  BHR_aggregate<-BHR(mode = "aggregate", 
                     ss_list = list(get(paste(pheno_name,"bin1_missense",sep = "_")),get(paste(pheno_name,"bin2_missense",sep = "_")),get(paste(pheno_name,"bin3_missense",sep = "_")),get(paste(pheno_name,"bin4_missense",sep = "_"))),
                     trait_list = list(pheno_name),
                     annotations = list(baseline),
                     num_blocks = 100, 
                     genomewide_correction=TRUE,
                     overdispersion=TRUE)
  aggregate_h2<-BHR_aggregate$aggregated_mixed_model_h2
  aggregate_h2_se<-BHR_aggregate$aggregated_mixed_model_h2se
  result<-cbind(pheno_name,aggregate_h2,aggregate_h2_se,"missense")
  colnames(result)<-c("pheno_name","aggregate_h2","aggregate_h2_se","maf_func_bins")
  result_aggregate<-rbind(result_aggregate,result)
}
result_aggregate<-result_aggregate[-1,]
write_path<-paste("/BHR/BHR_output/",pheno_name,"_aggregate.txt",sep = "")
write.table(result_aggregate,write_path,sep = "\t",row.names = F,quote = F)