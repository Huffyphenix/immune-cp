library(ggplot2)
library(magrittr)


# Path
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
methy_path <- "/home/cliu18/liucj/projects/6.autophagy/06_methylation"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
methy_box <- file.path(methy_path, "boxplot")
cnv_path <- "/home/cliu18/liucj/projects/6.autophagy/03_cnv"

methy <- readr::read_rds(file.path(methy_path, ".rds_02_methy_a_gene_list_methy.rds.gz"))
expr <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"))
cnv <- readr::read_rds(path = file.path(cnv_path, ".rds_02_cnv_a_gene_list_raw.rds.gz"))

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
}
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}

fn_transform <- function(.d){
  # .d <- te$filter_methy[[1]]
  .d %>% 
    # dplyr::select(-2) %>%
    tidyr::gather(key = barcode, value = value, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::select(-type) %>% 
    dplyr::mutate(sample = fun_barcode(barcode)) %>%
    dplyr::distinct(symbol, sample, .keep_all = T) %>% 
    dplyr::select(-barcode)
}

methy %>% 
  dplyr::mutate(filter_methy = purrr::map(.x = filter_methy, .f = fn_transform)) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(b_val = value) -> methy_df

expr %>% 
  dplyr::mutate(filter_expr = purrr::map(.x = filter_expr, .f = fn_transform)) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(expr = value) %>% 
  dplyr::mutate(expr = log2(expr + 1)) -> expr_df

cnv %>% 
  dplyr::mutate(filter_cnv = purrr::map(.x = filter_cnv, .f = fn_transform)) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(! is.na(value)) %>% 
  dplyr::rename(cnv = value) -> cnv_df

df_expr_methy <-
  methy_df %>% 
  dplyr::inner_join(expr_df, by = c("cancer_types", "symbol", "sample"))

df_expr_methy %>% 
  dplyr::group_by(cancer_types, symbol) %>% 
  dplyr::do(broom::glance(lm(expr ~ b_val, data = .))) %>% 
  dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cancer_types, symbol, ars = adj.r.squared, p.value, fdr) -> df_expr_methy_cor

df_expr_methy_cor %>%   
  ggplot(aes(x = ars, y = -log10(fdr))) +
  geom_point() +
  facet_wrap(~ cancer_types) 

df_expr_cnv <- 
  cnv_df %>% 
  dplyr::inner_join(expr_df, by = c("cancer_types", "symbol", "sample"))

df_expr_cnv %>%
  dplyr::group_by(cancer_types, symbol) %>% 
  dplyr::do(broom::glance(lm(expr ~ cnv, data = .))) %>% 
  dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cancer_types, symbol, ars = adj.r.squared, p.value, fdr) -> df_expr_cnv_cor

df_expr_cnv %>% 
  dplyr::filter(cancer_types == "BRCA", symbol == "ATG5") %>% 
  dplyr::group_by(cancer_types, symbol) %>% 
  dplyr::do(broom::tidy(lm(expr ~ cnv, data = .))) %>% 
  dplyr::filter(term == "cnv") %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cancer_types, symbol, estimate) -> df_expr_cnv_coef

df_expr_cnv_cor %>%   
  ggplot(aes(x = ars, y = -log10(fdr))) +
  geom_point() +
  facet_wrap(~ cancer_types) 

df_expr_cnv_methy <- 
  df_expr_methy %>% 
  dplyr::inner_join(cnv_df, by = c("cancer_types", "symbol", "sample"))

df_expr_cnv_methy %>% 
  # dplyr::filter(cancer_types == "BRCA", symbol == "ATG5") %>%
  dplyr::group_by(cancer_types, symbol) %>% 
  dplyr::do(broom::glance(lm(expr ~ cnv + b_val, data = .))) %>% 
  dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cancer_types, symbol, ars = adj.r.squared, p.value, fdr) -> df_expr_cnv_methy_cor

df_expr_cnv_methy_cor %>% 
  ggplot(aes(x = ars, y = -log10(fdr))) +
  geom_point()




save.image(file = file.path(methy_path, ".rda_meth_a_expr.rda"))














