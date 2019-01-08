library(magrittr)
library(methods)
path <- "/extraspace/TCGA/TCGA_protein_TCPA_L4"
files.names <- list.files("/extraspace/TCGA/TCGA_protein_TCPA_L4", pattern = "*-L4.csv")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"

# load rns
rppa_name_symbol <- readr::read_rds(file.path(tcga_path, "rppa_name_symbol.rds.gz"))

cancer_types <- 
  files.names %>% 
  stringr::str_split(pattern = "-", simplify = T) %>% .[,2]
cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)

fn_rppa_expr <- function(.x){
  path <- "/extraspace/TCGA/TCGA_protein_TCPA_L4"
  
  readr::read_csv(file = file.path(path, .x), progress = F) %>% 
    # dplyr::select(- dplyr::starts_with("X")) %>%
    dplyr::select(-c(Cancer_Type, Sample_Type, SetID)) %>% 
    dplyr::rename(barcode = Sample_ID) %>% 
    tidyr::gather(key = protein, value = rppa, -barcode) %>% 
    tidyr::spread(key = barcode, value = rppa) %>% 
    dplyr::mutate(symbol = stringr::str_split(protein, "_", simplify = T)[,1])
}

cancers_names %>%
  head(1) %>% 
  dplyr::mutate(
    protein_expr = purrr::map(names, fn_rppa_expr)
  ) 

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_rppa_expr", fn_rppa_expr)  %>%
  multidplyr::cluster_assign_value("rppa_name_symbol", rppa_name_symbol) %>% 
  dplyr::mutate(protein_expr = purrr::map(names, fn_rppa_expr)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan33_rppa_expr
parallel::stopCluster(cluster)

pancan33_rppa_expr %>% readr::write_rds(file.path(tcga_path, "pancan33_rppa_expr_l4.rds.gz"), compress = "gz")
