library(magrittr)
path <- "/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp"
files.names <- list.files("/extraspace/TCGA/TCGA_protein/", pattern="*re_20160627")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"

# load rns
rppa_name_symbol <- readr::read_rds(file.path(tcga_path, "rppa_name_symbol.rds.gz"))

cancer_types <- 
  files.names %>% 
  stringr::str_split(pattern = "\\_", simplify = T) %>% .[,1]
cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)

fn_rppa_expr <- function(.x, .rns){
  path <- "/extraspace/TCGA/TCGA_protein/"
  readr::read_tsv(file = file.path(path, .x), progress = F) %>% 
    dplyr::select(- dplyr::starts_with("X")) %>% 
    dplyr::left_join(.rns, by = "protein") 
}

cancers_names %>%
  head(1) %>% 
  dplyr::mutate(
    protein_expr = purrr::map(names, fn_rppa_expr, .rns = rppa_name_symbol)
  ) 

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_rppa_expr", fn_rppa_expr)  %>%
  multidplyr::cluster_assign_value("rppa_name_symbol", rppa_name_symbol) %>% 
  dplyr::mutate(protein_expr = purrr::map(names, fn_rppa_expr, .rns = rppa_name_symbol)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan_rppa_expr
on.exit(parallel::stopCluster(cluster))

pancan_rppa_expr %>% readr::write_rds(file.path(tcga_path, "pancan_rppa_expr.rds.gz"), compress = "gz")
  