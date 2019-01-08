library(magrittr)
path <- "/extraspace/yye1/share_data/TCGA_methy_mostNegatively"
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"

rds_files <- list.files(path = path, pattern = "rds")

cancers_names <- 
  tibble::tibble(names = rds_files) %>% 
  dplyr::mutate(cancer_types = stringr::str_split(rds_files, "_", simplify = T) %>% .[,1])

fn_process_rds <- function(.x, path = path){
  readr::read_rds(file.path(path, .x)) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(symbol = stringr::str_split(gene, "_", simplify = T) %>% .[,2]) -> .d
  names(.d) %>% stringr::str_replace_all("\\.", "-") -> names(.d)
  .d
}
cancers_names %>% 
  head(1) %>% 
  dplyr::mutate(methy = purrr::map(.x = names, .f = fn_process_rds, path = path)) %>% 
  dplyr::select(-names)

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_process_rds", fn_process_rds)  %>%
  multidplyr::cluster_assign_value("path", path) %>% 
  dplyr::mutate(methy = purrr::map(.x = names, .f = fn_process_rds, path = path)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan_meth
on.exit(parallel::stopCluster(cluster))

pancan_meth %>% readr::write_rds(file.path(tcga_path, "pancan33_meth.rds.gz"), compress = "gz")