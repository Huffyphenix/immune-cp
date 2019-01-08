`%>%` <- magrittr::`%>%`
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
files.names <- list.files("/extraspace/yye1/share_data/normalized.ccle.tcga.exp", pattern = "*.csv")
cancer_types <- files.names %>% stringr::str_split(pattern = "\\.|\\_", simplify = T) %>% .[,1]
cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)


process_raw_data <- function(.x){
  print(.x)
  path <- "/extraspace/yye1/share_data/normalized.ccle.tcga.exp"
  
  d <- 
    readr::read_csv(file = file.path(path, .x), progress = F) %>% 
    dplyr::select(- dplyr::starts_with("X"))  
  names(d)[1] <- "refseq_nm"
  d
}

cancers_names %>%
  dplyr::slice(7) %>% 
  dplyr::mutate(cnv = purrr::map(names, process_raw_data)) -> te

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  # dplyr::slice(7) %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_assign_value("process_raw_data", process_raw_data) %>% 
  multidplyr::cluster_assign_value("%>%", magrittr::`%>%`) %>% 
  dplyr::mutate(expr = purrr::map(names, process_raw_data)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan_ccle_expr
on.exit(parallel::stopCluster(cluster))
pancan_ccle_expr %>% readr::write_rds(path = file.path(tcga_path, 'pancan_ccle26_normalized_expr.rds.gz'), compress = "gz")








