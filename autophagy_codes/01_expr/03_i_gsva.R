library(magrittr)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
#output path
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

#
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
gene_list_expr <- readr::read_rds(file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"))
gene_list_fc_pvalue_simplified <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_a_atg_lys_fc_pvalue_simplified.rds.gz"))

gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

gene_sets <- list(
  atg = dplyr::filter(gene_list, type == "atg"),
  lys = dplyr::filter(gene_list, type == "lys"),
  pathway = dplyr::filter(gene_list, type == "pathway"),
  atg_core = dplyr::filter(gene_list, pathway == "autophagesome formation-core"),
  atg_form = dplyr::filter(gene_list, pathway == "autophagesome formation"),
  lys_com = dplyr::filter(gene_list, pathway == "lysosome composition"),
  lys_deg = dplyr::filter(gene_list, pathway == "lysosome degradation")
) %>% purrr::map("symbol")

# use GSVA to calculate score for every part
library(GSVA)

fn_gsva <- function(.x, .y, gene_sets = gene_sets){
  # .x <- .te$cancer_types
  # .y <- .te$filter_expr[[1]]
  print(.x)
  .y %>% 
    tidyr::drop_na() %>% 
    dplyr::select( -entrez_id) -> .d
  
  .d_mat <- as.matrix(.d[,-1])
  rownames(.d_mat) <- .d$symbol
  
  .es_dif <- gsva(.d_mat, gene_sets, method = "gsva", mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
  
  .es_dif$es.obs %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(set = rownames(.es_dif$es.obs), .before = 1) -> .d_es
}

expr %>% 
  dplyr::mutate(gsva = purrr::map2(.x = cancer_types, .y = expr, .f = fn_gsva, gene_sets = gene_sets)) %>% 
  dplyr::select(-expr) -> gene_list_gsva

readr::write_rds(x = gene_list_gsva, path = file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"), compress = "gz")

calculate_fc_pvalue <- function(.x, .y) {
  .y %>% tibble::add_column(cancer_types = .x, .before = 1) -> df
  
  # get cancer types and get # of smaple >= 10
  samples <-
    tibble::tibble(barcode = colnames(df)[-c(1:3)]) %>%
    dplyr::mutate(
      sample = stringr::str_sub(
        string = barcode,
        start = 1,
        end = 12
      ),
      type = stringr::str_split(barcode, pattern = "-", simplify = T)[, 4] %>% stringr::str_sub(1, 2)
    ) %>%
    dplyr::filter(type %in% c("01", "11")) %>%
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("01" = "Tumor", "11" = "Normal"),
      warn_missing = F
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
    dplyr::ungroup()
  sample_type_summary <- table(samples$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) ||
      any(sample_type_summary < c(10, 10))) {
    return(NULL)
  }
  
  # filter out cancer normal pairs
  df_f <-
    df %>%
    dplyr::select(c(1, 2, 3), samples$barcode) %>%
    tidyr::gather(key = barcode, value = expr, -c(1, 2, 3)) %>%
    dplyr::left_join(samples, by = "barcode")
  
  # pvalue & fdr
  df_f %>%
    dplyr::group_by(cancer_types, set) %>%
    tidyr::drop_na(expr) %>%
    dplyr::do(broom::tidy(t.test(expr ~ type, data = .))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    dplyr::select(cancer_types, set,  p.value, fdr) -> df_pvalue
  
  # log2 fold change mean
  df_f %>%
    dplyr::group_by(cancer_types, set,  type) %>%
    tidyr::drop_na(expr) %>%
    dplyr::summarise(mean = mean(expr)) %>%
    tidyr::spread(key = type, mean) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fc = (Tumor + 0.1) / (Normal + 0.1)) -> df_fc
  
  df_fc %>%
    dplyr::inner_join(df_pvalue, by = c("cancer_types", "set")) %>%
    dplyr::mutate(n_normal = sample_type_summary[1], n_tumor = sample_type_summary[2]) -> res
  return(res)
}

purrr::map2(.x = gene_list_gsva$cancer_types,
            .y = gene_list_gsva$gsva,
            .f = calculate_fc_pvalue) -> gene_list_gsva_fc

library(ggplot2)
gene_list_gsva_fc %>% 
  dplyr::bind_rows() %>% 
  ggplot(aes(x = set, y = cancer_types, size = fc)) +
  geom_point()

fn_diff <- function(.gsva){
  # .gsva <- .te$gsva[[1]]
  
  .samples <-
    tibble::tibble(barcode = colnames(.gsva)[-1]) %>% 
    dplyr::mutate(
      sample = stringr::str_sub(
        string = barcode,
        start = 1,
        end = 12
      ),
      type = stringr::str_split(barcode, pattern = "-", simplify = T)[, 4] %>% stringr::str_sub(1,2)
    ) %>% 
    dplyr::filter(type %in% c("01", "11")) %>% 
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("01" = "Tumor", "11" = "Normal"),
      warn_missing = F
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
    dplyr::ungroup()
  
  .sample_type_summary <- table(.samples$type) %>% as.numeric()
  
  if (gtools::invalid(.sample_type_summary) ||
      any(.sample_type_summary < c(10, 10))) {
    return(NULL)
  }
  
  .gsva %>% 
    tidyr::gather(key = barcode, value = es, - set) %>% 
    dplyr::inner_join(.samples, by = "barcode") -> .gsva_cr
  
  .gsva_cr %>% 
    dplyr::group_by(set) %>% 
    tidyr::drop_na(es) %>% 
    dplyr::do(broom::tidy(t.test(es ~ type, data = .))) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::select(set, estimate1, estimate2, statistic, p.value, fdr)
  
}
gene_list_gsva %>%
  dplyr::mutate(gsva_diff = purrr::map(.x = gsva, .f = fn_diff)) -> gene_list_gsva_diff

# no expression difference between tumor and normal
gene_list_gsva_diff %>% 
  dplyr::filter(! purrr::map_lgl(.x = gsva_diff, .f= is.null)) %>% 
  tidyr::unnest(gsva_diff) %>% 
  dplyr::filter(fdr < 0.05) %>% 
  ggplot(aes(x = cancer_types, y = set)) +
  geom_point(aes(size = -log10(fdr)))

# autophagy gene and lysosome gene are negative correlated
gene_list_gsva %>% 
  dplyr::filter(cancer_types == "THCA") %>% 
  dplyr::pull(gsva) %>% 
  .[[1]] %>% 
  tidyr::gather(key = barcode, value = gsva, - set) %>% 
  tidyr::spread(key = set, value = gsva) %>% 
  dplyr::filter(stringr::str_sub(barcode, start = 14, end = 15) == "01") %>% 
  dplyr::select(- barcode) %>% 
  cor() -> .te_cor
ggcorrplot::ggcorrplot(corr = .te_cor)

gene_list_gsva %>% 
  dplyr::filter(cancer_types == "PRAD") %>% 
  dplyr::pull(gsva) %>% 
  .[[1]] %>% 
  tidyr::gather(key = barcode, value = gsva, - set) %>% 
  tidyr::spread(key = set, value = gsva) %>% 
  dplyr::filter(stringr::str_sub(barcode, start = 14, end = 15) == "01") %>% 
  ggplot(aes(x = reorder(barcode, atg_core), y = atg_core)) +
  geom_point()

save.image(file = file.path(expr_path_a, ".rda_03_i_gsva.rda"))
load(file = file.path(expr_path_a, ".rda_03_i_gsva.rda"))

