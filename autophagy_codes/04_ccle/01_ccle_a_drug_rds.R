library(magrittr)
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
drug_path <- "/extraspace/yye1/share_data/DrugData"
ctrp_drug <- readr::read_rds(file.path(drug_path, "CTRP_Drug_Exp.spearman.rds")) %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(symbol = GeneSymbol)
gdsc_drug <- readr::read_rds(file.path(drug_path, "GDSC_Drug_Exp.spearman.rds")) %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(symbol = GDSC.exp.GENE_SYMBOLS)


fn_gdsc_reform <- function(.x){
  .x %>% 
    dplyr::mutate(drug_name = stringr::str_replace(string = drug, pattern = "_FDR", replacement = "")) %>% 
    dplyr::mutate(drug = ifelse(stringr::str_detect(string = drug, pattern = "_FDR"), "fdr", "sprm")) %>% 
    tidyr::spread(key = drug, value = sprm) %>% 
    dplyr::select(drug_name, cor_sprm = sprm, fdr)
}
gdsc_drug %>% 
  head(100) %>% 
  dplyr::filter(!is.na(symbol)) %>% 
  # dplyr::group_by(symbol) %>% dplyr::filter(n() > 1)
  tidyr::gather(key = drug, value = sprm, -symbol) %>% 
  tidyr::nest(-symbol, .key = drug) %>% # .$data %>% .[[1]] -> .x
  dplyr::mutate(drug = purrr::map(.x = drug, .f = fn_gdsc_reform))

fn_ctrp_reform <- function(.x){
  .x %>% 
    dplyr::mutate(drug_name = stringr::str_replace(string = drug, pattern = "_p", replacement = "")) %>% 
    dplyr::mutate(drug = ifelse(stringr::str_detect(string = drug, pattern = "_p"), "p_val", "sprm")) %>% 
    tidyr::spread(key = drug, value = sprm) %>% 
    dplyr::select(drug_name, cor_sprm = sprm, p_val)
}
ctrp_drug %>% 
  head(1200) %>%
  dplyr::filter(symbol != "") %>% 
  dplyr::distinct(symbol, .keep_all = T) %>% 
  # dplyr::group_by(symbol) %>% dplyr::filter(n() > 1) -> te
  tidyr::gather(key = drug, value = sprm, -symbol) %>% 
  tidyr::nest(-symbol, .key = drug) %>%  #.$drug %>% .[[1]] -> .x
  dplyr::mutate(drug = purrr::map(.x = drug, .f = fn_ctrp_reform)) 
  
  

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gdsc_drug %>% #head(2) %>% 
  dplyr::filter(!is.na(symbol)) %>% 
  tidyr::gather(key = drug, value = sprm, -symbol) %>%
  tidyr::nest(-symbol, .key = drug) %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_gdsc_reform", fn_gdsc_reform) %>%
  # multidplyr::cluster_assign_value("fn_gdsc_multidplyr", fn_gdsc_multidplyr)%>%
  # dplyr::rowwise() %>% 
  # tidyr::gather(key = drug, value = sprm, -symbol) %>% 
  # tidyr::nest(-symbol, .key = drug) %>%
  dplyr::mutate(drug = purrr::map(.x = drug, .f = fn_gdsc_reform)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> gdsc_drug_reform
gdsc_drug_reform %>% readr::write_rds(path = file.path(tcga_path, "drug_gdsc_exp_spearman.rds.gz"), compress = "gz")


ctrp_drug %>% #head(2) %>%
  dplyr::filter(symbol != "") %>% 
  dplyr::distinct(symbol, .keep_all = T) %>% 
  tidyr::gather(key = drug, value = sprm, -symbol) %>%
  tidyr::nest(-symbol, .key = drug) %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_ctrp_reform", fn_ctrp_reform) %>%
  # multidplyr::cluster_assign_value("fn_gdsc_multidplyr", fn_gdsc_multidplyr)%>%
  # dplyr::rowwise() %>% 
  # tidyr::gather(key = drug, value = sprm, -symbol) %>% 
  # tidyr::nest(-symbol, .key = drug) %>%
  dplyr::mutate(drug = purrr::map(.x = drug, .f = fn_ctrp_reform)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> ctrp_drug_reform
ctrp_drug_reform %>% readr::write_rds(path = file.path(tcga_path, "drug_ctrp_exp_spearman.rds.gz"), compress = "gz")
on.exit(parallel::stopCluster(cluster))
