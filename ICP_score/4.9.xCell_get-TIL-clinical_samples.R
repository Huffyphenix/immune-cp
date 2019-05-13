##################################
# TIL of clinical samples
# by using xCell
#################################
library(magrittr)
library(tidyverse)
library(caret)
library(glmnet)
library(broom)
library(pROC)

theme_set(theme_classic())

# data path ---------------------------------------------------------------
# server 1 basic path
basic_path <- file.path("/home/huff/project")

# E zhou basi path
basic_path <- file.path("F:/我的坚果云")

# Hust basi path
basic_path <- file.path("S:/坚果云/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical/class_metastic_type")

# load data ---------------------------------------------------------------
exp_data <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data/mRNA_exp","all_FPKM_expression_2.txt"))
# gene id transfer --------------------------------------------------------
ensembl_gene <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","Homo_sapiens.gene_info.geneid.symbol.ensembl"))
exp_data %>%
  dplyr::rename("ens_id"="gene_id") %>%
  dplyr::inner_join(ensembl_gene,by="ens_id") %>%
  dplyr::select(-ens_id,-GeneID) %>%
  dplyr::select(Symbol,1:(ncol(exp_data)-1)) -> exp_data


# get expression matrix for xCell -----------------------------------------
exp_data.m <- as.matrix(exp_data[,-1])
rownames(exp_data.m) <- exp_data$Symbol
exp_data.m[is.na(exp_data.m)] <- 0


# run xCell tool to calculate the xCell score for each sample --------------------------------
library(xCell)

fn_xCell <- function(.x,cell){
  exprMatrix = .x
  xCellAnalysis(exprMatrix,cell.types.use = cell)
}

cell_list <- c("B-cells","CD4+ T-cells","CD8+ T-cells","aDC","iDC","DC","Neutrophils","Macrophages")

clinical_xcell_results <- fn_xCell(exp_data.m,cell_list)

clinical_xcell_results %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(Run = colnames(clinical_xcell_results)) %>%
  readr::write_tsv(file.path("F:/我的坚果云/immune_checkpoint/result_20171025/ICP_score","xCell_TIL_of_clinical_data.tsv"))
