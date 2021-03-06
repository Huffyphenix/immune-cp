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
  data %>%
    .$log2FC -> .fc
  names(.fc) <- data$entrez_id
  sort(.fc,decreasing = TRUE) -> .fc
  .fc
}

fn_GSEA <- function(data){
  data %>%
    fn_data_deal() -> .data_ready
  
  .gsea_res <- gseKEGG(geneList     = .data_ready,
                       organism     = 'hsa',
                       nPerm        = 1000,
                       minGSSize    = 120,
                       pvalueCutoff = 1,
                       verbose      = FALSE)
  
  gseaplot2(kk2, geneSetID = "?") #use self-define gene list
  
}

gene_DE_between_high_low.all_features.cox_score.OS %>%
  dplyr::mutate(GSEA_res = purrr::map(DE_res,fn_GSEA))