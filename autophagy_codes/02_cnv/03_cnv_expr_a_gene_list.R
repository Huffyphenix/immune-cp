`%>%` <- magrittr::`%>%`
library(ggplot2)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
cnv_path <- "/home/cliu18/liucj/projects/6.autophagy/03_cnv"

# load data
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))
gene_list_expr <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"))
gene_list_cnv <- readr::read_rds(path = file.path(cnv_path, ".rds_02_cnv_a_gene_list.rds.gz"))

cnv_expr <- 
  gene_list_cnv %>% 
  dplyr::inner_join(gene_list_expr, by = "cancer_types")

fn_merge_cnv_expr <- function(.x, .y){
  #.x <- te$filter_cnv[[1]]
  .x %>% 
    tidyr::gather(key = barcode, value = cnv, -symbol) %>% 
    dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 16)) %>% 
    dplyr::select(-barcode) %>% 
    dplyr::mutate(type = stringr::str_sub(sample, start = 14, end = 16)) %>% 
    dplyr::mutate(sample = stringr::str_sub(sample, start = 1, end = 12)) %>% 
    dplyr::filter(type %in% c("01A")) %>% 
    dplyr::select(-type) ->
    .x_cnv
  #.y <- te$filter_expr[[1]]
  .y %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 16)) %>% 
    dplyr::select(-barcode) %>% 
    dplyr::mutate(type = stringr::str_sub(sample, start = 14, end = 16)) %>% 
    dplyr::mutate(sample = stringr::str_sub(sample, start = 1, end = 12)) %>% 
    dplyr::filter(type %in% c("01A", "11A")) %>% 
    dplyr::mutate(type = dplyr::recode(type, "01A" = "Tumor", "11A" = "Normal")) %>% 
    tidyr::spread(key = type, value = expr) %>% 
    tidyr::drop_na() ->
    .y_expr
  
  .x_cnv %>% 
    dplyr::inner_join(.y_expr, by = c("symbol", "sample")) ->
    .merge
  
  .merge %>% 
    dplyr::group_by(symbol, cnv) %>% 
    dplyr::count() %>% 
    dplyr::group_by(symbol) %>%  
    dplyr::summarise(s = sum(cnv * n))  %>% 
    dplyr::arrange(s) %>% 
    dplyr::left_join(gene_list, by = "symbol") %>%
    dplyr::select(1,2,3) %>% 
    dplyr::filter(pathway == "autophagesome formation-core") %>% 
    print(n = Inf)
  
  .merge %>% 
    tidyr::drop_na() %>% 
    dplyr::filter(symbol == "MAP1LC3A") %>% #dplyr::pull(cnv) %>% table
    tidyr::gather(key = type, value = expr, Normal, Tumor) %>% 
    dplyr::mutate(cnv = dplyr::case_when(cnv < 0 ~ -1, cnv > 0 ~ 1, TRUE ~ 0)) %>% 
    ggplot(aes(x = as.factor(cnv), y = log2(expr + 1))) +
    geom_boxplot(aes(color =type))
  
  .merge %>% 
    tidyr::drop_na() %>% 
    dplyr::filter(symbol == "MAP1LC3C") %>% 
    tidyr::gather(key = type, value = expr, Normal, Tumor) %>% 
    ggplot(aes(x = as.factor(type), y = log2(expr + 1))) +
    geom_boxplot(aes(color = type))
  
}

cnv_expr %>% 
  dplyr::filter(cancer_types == "BRCA") -> te
  dplyr::mutate(merge_cnv_expr = dplyr::map2(filter_cnv, filter_expr, .f = fn_merge_cnv_expr))

