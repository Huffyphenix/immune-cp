library(magrittr)
files.names <- list.files(path = "/extraspace/TCGA/Mutation", pattern = ".txt")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
cancer_types <- files.names %>% stringr::str_split(pattern = "\\.|\\_|-", simplify = T) %>% .[,2]
cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)

process_raw_data <- function(.x){
  print(.x)
  path <- "/extraspace/TCGA/Mutation"
  d <- readr::read_tsv(file = file.path(path, .x), progress = F)
  names(d)[1] <- "symbol"
  d %>% 
    dplyr::select(- dplyr::starts_with("X"))
}

cancers_names %>% 
  head(1) %>%
  dplyr::mutate(snv = purrr::map(names, process_raw_data))


cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("process_raw_data", process_raw_data) %>% 
  dplyr::mutate(snv = purrr::map(names, process_raw_data)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan_snv
on.exit(parallel::stopCluster(cluster))
pancan_snv %>% readr::write_rds(path = file.path(tcga_path, 'pancan33_snv.rds.gz'), compress = "gz")
