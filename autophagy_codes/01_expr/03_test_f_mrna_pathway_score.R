`%>%` <- magrittr::`%>%`
library(ggplot2)
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
test_f_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/03_test_f_mrna_pathway"

expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
gene_list <- readr::read_rds(file.path(test_f_path, "rds_03_test_f_mrna_pathway_score.rds.gz"))

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr


fn_z_score <- function(.d){
  .d %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    tidyr::replace_na(replace = list(expr = 0)) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(z_score = (expr - median(expr)) / sd(expr)) %>% 
    dplyr::ungroup() -> d_z
  
  d_z %>% 
    dplyr::group_by(barcode) %>% 
    dplyr::summarise(expr = sum(expr)) %>% 
    dplyr::ungroup() -> d_s
    
  d_s %>% 
    ggplot(aes(x = expr)) +
    geom_histogram()
}

gene_list_expr %>% 
  dplyr::mutate(normlized = purrr::map(.x = filter_expr, .f = fn_z_score))
