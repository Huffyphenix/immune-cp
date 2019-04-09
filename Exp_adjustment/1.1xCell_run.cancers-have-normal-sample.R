# Get immune fraction from immune + stromal data ---------------------------------

library(magrittr)

# data path ---------------------------------------------------------------
basic_path <- "/home/huff/project"

xcell_path <- file.path(basic_path,"data/TCGA/immune_infiltration/xCell")

purity_path <- file.path(basic_path,"data/TCGA_tumor_purity_from_ncomms9971")
immune_path <- file.path(basic_path,"immune_checkpoint")
tcga_path <- file.path(immune_path,"data/TCGA_data")


# load data ---------------------------------------------------------------

TCGA_purity <- readr::read_tsv(file.path(purity_path,"ncomms9971-s2.txt")) %>%
  dplyr::select(`Cancer type`,`Sample ID`,`CPE`) %>%
  dplyr::mutate(`Sample ID`=substr(`Sample ID`,1,15))

all_TCGA_exp <- readr::read_rds(file.path(tcga_path,"pancan33_expr.rds.gz")) 

### Importantly, xCell performs best with heterogenous dataset. Thus it is recommended to use all data combined in one run, and not break down to pieces (especially not cases and control in different runs).
### So there is no need to run samples not in the xCell data set.

cancers_with_paired_normal <- c("KIRC","KIRP","HNSC","ESCA","STAD","BRCA","KICH","THCA","LIHC","COAD","LUAD","PRAD","BLCA","LUSC")
all_TCGA_exp %>%
  dplyr::filter(cancer_types %in% cancers_with_paired_normal) -> all_TCGA_exp

fn_get_subset_samples <- function(.cancer,.x){
  .x %>%
    dplyr::filter(!symbol =="?") %>%
    dplyr::select(-entrez_id) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate_all(mean)%>%  # rows with same gene symbols will be averaged 
    unique() %>%
    t()  %>%
    as.data.frame()-> .x.t
  .x.t %>%
    # dplyr::select(barcode,everything()) %>%
    t() %>%
    as.data.frame()-> .x.fiter
  rownames(.x.fiter) <- t(.x.t["symbol",])
  .x.fiter <- .x.fiter[,-1]
  
  print(.cancer)
  .x.fiter
}


all_TCGA_exp %>%
  # head(1) %>%
  dplyr::mutate(expr_filter=purrr::map2(cancer_types,expr,.f=fn_get_subset_samples)) %>%
  dplyr::select(-expr) -> smaple_exp.matrix

# run xCell tool to calculate the xCell score for each sample --------------------------------
library(xCell)

fn_xCell <- function(.x){
  exprMatrix = .x
  xCellAnalysis(exprMatrix)
}
library(multidplyr)
start <- proc.time() # Start clock
library(parallel)
cl <- 23
cluster <- create_cluster(cores = cl)
smaple_exp.matrix %>%
  partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("tidyverse") %>%
  multidplyr::cluster_library("stringr") %>%
  multidplyr::cluster_library("xCell") %>%
  multidplyr::cluster_assign_value("fn_xCell",fn_xCell) %>%
  dplyr::mutate(xCell = purrr::map(expr_filter,fn_xCell)) %>%
  collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -expr_filter) -> Pan14_xCell_results
on.exit(parallel::stopCluster(cluster))

Pan14_xCell_results %>%
  readr::write_rds(file.path(xcell_path,"Pan14_xCell_results.rds.gz"),compress = "gz")
