`%>%` <- magrittr::`%>%`

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
files.names <- list.files("/extraspace/yye1/share_data/TCGA_CNV_GISTIC", pattern = "GDAC_*")

cancer_types <- files.names %>% stringr::str_split(pattern = "\\.|\\_", simplify = T) %>% .[,2]
cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)

process_raw_data <- function(.x){
  print(.x)
  path <- "/extraspace/yye1/share_data/TCGA_CNV_GISTIC"

  d <- 
    readr::read_tsv(file = file.path(path, .x, "all_data_by_genes.txt"), progress = F) %>% 
    dplyr::select(- dplyr::starts_with("X")) %>% 
    dplyr::rename(symbol = `Gene Symbol`) %>% 
    dplyr::select(-`Locus ID`, -Cytoband)
}

cancers_names %>%
  head(1) %>%
  dplyr::mutate(cnv = purrr::map(names, process_raw_data)) -> te


cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("process_raw_data", process_raw_data) %>% 
  dplyr::mutate(cnv = purrr::map(names, process_raw_data)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan_cnv
on.exit(parallel::stopCluster(cluster))
pancan_cnv %>% readr::write_rds(path = file.path(tcga_path, 'pancan34_cnv.rds.gz'), compress = "gz")
