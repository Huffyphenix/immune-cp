# Get immune fraction from immune + stromal data ---------------------------------

library(magrittr)

# data path ---------------------------------------------------------------
xcell_path <- "/project/huff/huff/data/TCGA/immune_infiltration/xCell"
purity_path <- "/project/huff/huff/data/TCGA_tumor_purity_from_ncomms9971"
immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")

smaple_exp.matrix <- 
  readr::read_rds(file.path(xcell_path,"Pancan21_expr.for_xCell.rds.gz"))


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
  dplyr::select(-PARTITION_ID, -expr_filter) -> Pan21_xCell_results
on.exit(parallel::stopCluster(cluster))

Pan21_xCell_results %>%
  readr::write_rds(file.path(xcell_path,"Pan21_xCell_results.rds.gz"),compress = "gz")