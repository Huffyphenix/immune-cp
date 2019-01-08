library(magrittr)
path <- "/extraspace/yye1/tumor_purity_summary"
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
rds_files <- list.files(path = path, pattern = "txt")

cancers_names <- 
  tibble::tibble(names = rds_files) %>% 
  dplyr::mutate(cancer_types = stringr::str_split(rds_files, "[^[:alnum:]]+", simplify = T) %>% .[,2]) %>% 
  dplyr::mutate(cancer_types = stringr::str_to_upper(cancer_types))

fn_process_agp <- function(.x, path = path){
  readr::read_tsv(file.path(path, .x), progress = F) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(sampleid = stringr::str_replace_all(sampleid, "\\.", "-")) %>% 
    dplyr::rename(barcode = sampleid)
}
cancers_names %>% 
  head(1) %>% 
  dplyr::mutate(methy = purrr::map(.x = names, .f = fn_process_agp, path = path)) %>% 
  dplyr::select(-names)


cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_process_agp", fn_process_agp)  %>%
  multidplyr::cluster_assign_value("path", path) %>% 
  dplyr::mutate(methy = purrr::map(.x = names, .f = fn_process_agp, path = path)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan24_tumor_purity
on.exit(parallel::stopCluster(cluster))

pancan24_tumor_purity %>% readr::write_rds(file.path(tcga_path, "pancan24_tumor_purity.rds.gz"), compress = "gz")
