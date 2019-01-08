library(magrittr)
path <- "/extraspace/TCGA/TCGA_exp_DataPortal/miRNA_exp"
files.names <- list.files("/extraspace/TCGA/TCGA_exp_DataPortal/miRNA_exp", pattern="*_miRNA_each_exp_20160513")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"

mir_acc_name <- 
  readr::read_rds(file.path(tcga_path, "mirna_acc_name.rds.gz")) %>% 
  dplyr::select(gene = ID, name = Name)
cancers_names <- 
  tibble::tibble(names = files.names) %>% 
  dplyr::mutate(cancer_types = stringr::str_split(names, "_", simplify = T) %>% .[,1])

fn_process_raw <- function(.x, path = path, man = mir_acc_name){
  readr::read_tsv(file.path(path, .x), progress = F) %>% 
    dplyr::select(- dplyr::starts_with("X")) %>% 
    dplyr::mutate(gene = stringr::str_replace_all(gene, "mature,", "")) %>% 
    dplyr::filter(stringr::str_detect(gene, "MIMAT")) -> .d
  
  .d %>% 
    dplyr::left_join(man, by = "gene")
}

cancers_names %>% 
  head(1) %>% 
  dplyr::mutate(mirna = purrr::map(.x = names, .f = fn_process_raw, path = path, man = mir_acc_name)) %>% 
  dplyr::select(-names)

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_process_raw", fn_process_raw)  %>%
  multidplyr::cluster_assign_value("mir_acc_name", mir_acc_name)  %>%
  multidplyr::cluster_assign_value("path", path)  %>%
  dplyr::mutate(
    mirna = purrr::map(
    .x = names, 
    .f = fn_process_raw, 
    path = path, 
    man = mir_acc_name)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(names, PARTITION_ID)) -> pancan33_mirna_expr
parallel::stopCluster(cluster)

pancan33_mirna_expr %>% readr::write_rds(file.path(tcga_path, "pancan33_mirna_expr.rds.gz"), compress = "gz")



