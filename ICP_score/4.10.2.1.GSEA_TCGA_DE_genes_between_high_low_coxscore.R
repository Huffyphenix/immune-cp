############## GSEA enrichment of DE genes between high(higher than 50th percentile) and low(lower than 50th percentile) cox score
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(GSEABase)
library(enrichplot)

# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
DE_path <- file.path(immune_res_path,"TCGA_GSVAScore/cancer_specific/3.survival_with_GSVA_score.new/DE_analysis")
res_path <- file.path(DE_path,"GSEA")


# load data ---------------------------------------------------------------

gene_DE_between_high_low.all_features.cox_score.OS <-
  readr::read_rds(file.path(DE_path,"gene_DE_between_high_low.all_features.cox_score.OS.rds.gz"))

gene_DE_between_high_low.all_features.cox_score.OS$DE_res[[1]]


# function ----------------------------------------------------------------
fn_data_deal <- function(data){
  
}

fn_GSEA <- function(){
  
}
