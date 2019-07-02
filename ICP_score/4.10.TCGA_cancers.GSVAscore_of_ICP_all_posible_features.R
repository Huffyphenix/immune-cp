########################## get GSVA score of each cancers
########################## features are: expression site of ICPs, gene family
########################## and get the correlation between GSVA score and TIL, mutation burden
library(magrittr)
library(tidyverse)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(GSVA)
# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"TCGA_GSVAScore")

# load data ---------------------------------------------------------------
## all TCGA expression data 
exp_data <- readr::read_rds(file.path(TCGA_path,"pancan33_expr.rds.gz"))
gene_list <- readr::read_tsv(file.path(gene_list_path, "ICPs_all_info_class.tsv")) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Immune","Mainly_exp_on_Immune"),"Mainly_exp_on_Immune",Exp_site)) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Tumor","Mainly_exp_on_Tumor" ),"Mainly_exp_on_Tumor",Exp_site.1)) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("N"),"Not_sure",Exp_site.1)) 

# 1.get GSVA score of all possible features of ICP ----------------------------------------------
# 1.1. get gene feature from gene list ----
genelist <- list()
#### gene list feature by gene exp site #####
for(expsite in c("Both_exp_on_Tumor_Immune","Mainly_exp_on_Immune","Mainly_exp_on_Tumor")){
  genelist[[expsite]] <- gene_list %>%
    dplyr::filter(Exp_site.1 == expsite) %>%
    .$symbol
}

#### gene list feature by gene function in immune system #####
for(fun_class in unique(gene_list$functionWithImmune)){
  genelist[[fun_class]] <- gene_list %>%
    dplyr::filter(functionWithImmune == fun_class) %>%
    .$symbol
}
#### gene list feature by ligand-receptor pairs #####
for(p in unique(gene_list$Recepter_pairs)){
  if(!is.na(p)){
    genelist[[p]] <- gene_list %>%
      dplyr::filter(Recepter_pairs == p) %>%
      .$symbol
  }
  
}
#### gene list feature by gene family #####
for(f in unique(gene_list$family)){
  genelist[[f]] <- gene_list %>%
    dplyr::filter(family == f) %>%
    .$symbol
}
genelist[["All_gene"]] <- gene_list$symbol

# 1.2.calculation function -----------------------------------------------
fn_GSVA <- function(tissue, exp){
  print(paste("GSVA",tissue))
  #### do GSVA ####
  if (length(genelist)>0) {
    exp <- as.matrix(exp)
    rownames.exp <- exp[,1]
    exp <- exp[,-1]
    exp <- apply(exp, 2, as.numeric)
    rownames(exp) <- rownames.exp
    res.gsva <- gsva(exp,genelist, mx.diff = FALSE, 
                     method = c("gsva"),
                     kcdf = c("Gaussian"), 
                     verbose = FALSE, parallel.sz = 1)
    res.gsva %>%
      as.data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::mutate(feature = rownames(res.gsva)) %>%
      tidyr::gather(-feature, key = "barcode",value = "GSVA_score") %>%
      tidyr::spread(key = "feature", value = "GSVA_score") -> gsva.score
  } else{
    gsva.score <- tibble::tibble()
  }
  
  # fn_compare_TN(gsva.score, cancer_types, result_path = file.path(res_path,"TN_compare"))
  
  # fn_survival.calculate(gsva.score, cancer_types)
  
  # fn_heatmap(res.gsva)
  gsva.score
}

fn_get_TCGA_sample_class <- function(.x,class){
  sample <- colnames(.x)[substr(colnames(.x),14,14) %in% class]
  
  data <- .x[,c("symbol",sample)]
  
  
  data
}
# 1.3.calculation ---------------------------------------------------------

# 1.3.1 tumor and normal samples togather ----------
# exp_data %>%
#   # head(1) %>%
#   dplyr::mutate(exp_filter = purrr::map(expr,.f=function(.x){
#     .x %>%
#       # dplyr::filter(symbol %in% gene_list$symbol) %>%
#       dplyr::select(-entrez_id)
#   })) %>%
#   dplyr::mutate(exp_data = purrr::map(exp_filter,fn_get_TCGA_sample_class,class=c("1","0"))) %>%
#   dplyr::select(-expr,-exp_filter) %>%
#   dplyr::mutate(GSVA = purrr::map2(cancer_types, exp_data, fn_GSVA)) %>%
#   dplyr::select(-exp_data) -> GSVA.score
# GSVA.score %>%
#   readr::write_rds(file.path(res_path,"TCGA_cancer_specific.allsamples(T-N)_GSVA.score_ICPs_features.rds.gz"),compress = "gz")

# 1.3.1 only tumor samples ----------
exp_data %>%
  # head(1) %>%
  dplyr::mutate(exp_filter = purrr::map(expr,.f=function(.x){
    .x %>%
      # dplyr::filter(symbol %in% gene_list$symbol) %>%
      dplyr::select(-entrez_id)
  })) %>%
  dplyr::mutate(exp_data = purrr::map(exp_filter,fn_get_TCGA_sample_class,class=c("0"))) %>%
  dplyr::select(-expr,-exp_filter) %>%
  dplyr::mutate(GSVA = purrr::map2(cancer_types, exp_data, fn_GSVA)) %>%
  dplyr::select(-exp_data) -> GSVA.score.onlytumor
GSVA.score.onlytumor %>%
  readr::write_rds(file.path(res_path,"TCGA_cancer_specific.onlyTumor_GSVA.score_ICPs_features.rds.gz"),compress = "gz")



