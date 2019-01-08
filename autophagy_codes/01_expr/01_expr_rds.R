path <- "/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp"
files.names <- list.files(path=path ,pattern="_20160513")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
library(magrittr)

process_raw_data <- function(.x){
  # print(file.path(path, .x))
  path <- "/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp"
  d <- readr::read_tsv(file = file.path(path, .x), progress = F)
  d %>% dplyr::select(- dplyr::starts_with("X")) %>%
    # split gene and id
    tidyr::separate(col = gene, into = c("symbol", "entrez_id"), sep = "\\|")
}

cancer_types <- files.names %>% stringr::str_split(pattern = "\\_", simplify = T) %>% .[,1]
# files.names  %>% purrr::map(process_raw_data, path = path) -> te
# names(te) <- cancer_types

cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)
cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("process_raw_data", process_raw_data) -> cancers_names_shards
cancers_names_shards %>%
  dplyr::mutate(expr = purrr::map(.x = names, process_raw_data)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(names, PARTITION_ID)) -> pancan_expr
parallel::stopCluster(cluster)


# pancan_expr %>% readr::write_rds(path = file.path(tcga_path, 'pancan_expr_20160513.rds.gz'), compress = "gz")
