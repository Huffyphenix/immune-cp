############## DE analysis between high and low cox score 
############## score estimate by sum(coef[coxph]*GSVA score)
library(magrittr)
library(tidyverse)
# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
score_path <- file.path(immune_res_path,"TCGA_GSVAScore/cancer_specific/3.survival_with_GSVA_score.new")
res_path <- file.path(score_path,"DE_analysis")


# load data ---------------------------------------------------------------

# expression data
expr <- readr::read_rds(file.path(TCGA_path, "pancan33_expr.rds.gz"))
PFS_cox_score.all_features <- readr::read_tsv(file.path(score_path, "3.PFS_cox_score(GSVA*coef_of_multi_cox)-all_features.tsv"))
OS_cox_score.all_features <- readr::read_tsv(file.path(score_path, "3.OS_cox_score(GSVA*coef_of_multi_cox)-all_features.tsv"))

# function to do DE analysis ---------------------------------------------
fn_DE <- function(.expr,score){
  .expr %>%
    # head(50) %>%
    dplyr::select(-symbol) %>%
    tidyr::gather(-entrez_id,key="barcode",value="exp") %>%
    dplyr::inner_join(score,by="barcode") %>%
    dplyr::select(-cancer_types) %>%
    dplyr::group_by(entrez_id) %>%
    dplyr::mutate(group = ifelse(cox_score > quantile(cox_score,0.5),"high","low")) %>%
    dplyr::ungroup() %>%
    tidyr::nest(-entrez_id) %>%
    dplyr::mutate(FC = purrr::map2(data,entrez_id,.f=function(.x,.y){
      print(.y)
      .x %>%
        dplyr::filter(group == "high") %>%
        .$exp %>%
        mean() -> high_mean
      .x %>%
        dplyr::filter(group == "low") %>%
        .$exp %>%
        mean() -> low_mean
      if(max(.x$exp)!=0){
        broom::tidy(
          wilcox.test(exp ~ group, data = .x, alternative = "two.sided")
        ) %>%
          dplyr::mutate(high_mean=high_mean,low_mean=low_mean,log2FC = log2((high_mean+1)/(low_mean+1)))
      }else{
        tibble::tibble()
      }
    })) %>%
    dplyr::select(-data) %>%
    tidyr::unnest()
}

# calculation ----
# PFS 
library(multidplyr)
library(parallel)
# start <- proc.time() # Start clock

# cl <- 16
# cluster <- create_cluster(cores = cl)
# expr %>%
#   partition(cluster = cluster) %>%
#   multidplyr::cluster_library("magrittr") %>%
#   multidplyr::cluster_library("tidyverse") %>%
#   multidplyr::cluster_library("stringr") %>%
#   multidplyr::cluster_assign_value("fn_DE",fn_DE) %>%
#   multidplyr::cluster_assign_value("PFS_cox_score.all_features",PFS_cox_score.all_features) %>%
#   dplyr::filter(cancer_types != "LAML") %>%
#   dplyr::mutate(DE_res = purrr::map(expr,fn_DE,score=PFS_cox_score.all_features)) %>%
#   collect() %>%
#   dplyr::as_tibble() %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-PARTITION_ID, -expr) -> gene_DE_between_high_low.all_features.cox_score.PFS
# on.exit(parallel::stopCluster(cluster))
# 
# gene_DE_between_high_low.all_features.cox_score.PFS %>%
#   readr::write_rds(file.path(res_path,"gene_DE_between_high_low.all_features.cox_score.PFS.rds.gz"),compress = "gz")

# OS
# start <- proc.time() # Start clock
cl <- 16
cluster <- create_cluster(cores = cl)
expr %>%
  partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("tidyverse") %>%
  multidplyr::cluster_library("stringr") %>%
  multidplyr::cluster_assign_value("fn_DE",fn_DE) %>%
  multidplyr::cluster_assign_value("OS_cox_score.all_features",OS_cox_score.all_features) %>%
  dplyr::filter(cancer_types != "LAML") %>%
  dplyr::mutate(DE_res = purrr::map(expr,fn_DE,score=OS_cox_score.all_features)) %>%
  collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -expr) -> gene_DE_between_high_low.all_features.cox_score.OS
on.exit(parallel::stopCluster(cluster))

gene_DE_between_high_low.all_features.cox_score.OS %>%
  readr::write_rds(file.path(res_path,"gene_DE_between_high_low.all_features.cox_score.OS.rds.gz"),compress = "gz")

