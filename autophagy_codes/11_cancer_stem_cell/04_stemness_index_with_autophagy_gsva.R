library(magrittr)
library(ggplot2)
tcga_dir <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
pcbc_dir <- "/home/cliu18/liucj/projects/6.autophagy/synapse/PCBC"
csc_dir <- "/home/cliu18/liucj/projects/6.autophagy/08_cancer_stem_cell"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

expr_sci_score <- readr::read_rds(path = file.path(csc_dir, ".rds_03_sample_score_expr_sci_score.rds.gz"))
gsva_score <- readr::read_rds(file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"))

gsva_score %>% dplyr::left_join(expr_sci_score, by = "cancer_types") -> scores_gsva_sci

fn_tn <- function(.s){
  .s %>% stringr::str_sub(start = 14, end = 15)
}

fn_test_corr <- function(.p, .d){
  # .p <- "atg_core"
  .v <- rlang::sym(.p)
  
  .dd <- .d %>% dplyr::filter(high_low %in% c("high", "low"))
  
  tryCatch(
    t.test(rlang::eval_bare(.v) ~ high_low, data = .dd),
    error = function(e){1}, warning = function(e){1}
    ) %>% 
    broom::tidy() -> .t_test
  
  tryCatch(
    cor.test(
    x = .d %>% dplyr::pull(csi) , 
    y = .d %>% dplyr::pull(rlang::eval_bare(.v)), 
    method = "spearman"
    ),
    error = function(e){1}, warning = function(e){1}
    ) %>% 
    broom::tidy() -> .cor_test
  
  tibble::tibble(
    pathway = .p, 
    t_diff = ifelse(tibble::has_name(.t_test, "estimate"), .t_test$estimate, 0), # diff is high - low
    t_pval = ifelse(tibble::has_name(.t_test, "p.value"), .t_test$p.value, 1), 
    cor_coef = .cor_test$estimate, 
    cor_pval = .cor_test$p.value
    )
}

fn_cor_gsva_sci <- function(.gsva, .sci){
  # .gsva <- .te$gsva[[1]]
  # .sci <- .te$sci[[1]]
  
  # use 1.5 sigma as cutoff
  .sigma <- c(mean(.sci$score) - 1 * sd(.sci$score), mean(.sci$score) + 1.5 * sd(.sci$score))
  
  .gsva %>% 
    tidyr::gather(key = barcode, value = gsva, -set) %>% 
    tidyr::spread(key = set, value = gsva) %>% 
    dplyr::inner_join(.sci, by = "barcode") %>% 
    dplyr::rename(csi = score) %>% 
    dplyr::mutate(type = fn_tn(barcode)) %>% 
    dplyr::filter(type != "11") %>% 
    # dplyr::mutate(high_low = dplyr::case_when(
    #   csi >= .sigma[2] ~ "high",
    #   csi <= .sigma[1] ~ "low",
  #   TRUE ~ "mid")) -> .d
    dplyr::mutate(high_low = dplyr::case_when(
      csi >= 0 ~ "high",
      TRUE ~ "low")) -> .d
  
  # .dt %>% dplyr::count(high_low) %>% print()
  .d %>% dplyr::count(type) %>% print()

  .gsva$set %>% 
    purrr::map_df(.f = fn_test_corr, .d = .d) %>% 
    dplyr::arrange(cor_coef)
}

scores_gsva_sci %>%
  # dplyr::filter(cancer_types == "BRCA") %>% 
  dplyr::mutate(pathway_diff = purrr::map2(.x = gsva, .y = sci, .f = fn_cor_gsva_sci)) ->
  scores_gsva_sci_cor

readr::write_rds(x = scores_gsva_sci_cor, path = file.path(csc_dir, ".rds_04_stemness_index_scores_gsva_sci_cor.rds.gz"), compress = "gz")


save.image(file = file.path(csc_dir, '.rds_04_stemness_index.rda'))
load(file = file.path(csc_dir, '.rds_04_stemness_index.rda'))

#-----------------------------------------------------------
scores_gsva_sci_cor <- readr::read_rds(path = file.path(csc_dir, ".rds_04_stemness_index_scores_gsva_sci_cor.rds.gz"))

scores_gsva_sci_cor %>% 
  tidyr::unnest(pathway_diff) %>% 
  dplyr::filter(cor_pval < 0.05) %>%
  ggplot(aes(cancer_types, y = pathway)) +
  geom_point(aes(size = -log10(cor_pval), color = cor_coef)) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue") +
  theme(
    axis.text.x = element_text(angle = 90)
  )

scores_gsva_sci_cor %>% 
  tidyr::unnest(pathway_diff) %>% 
  dplyr::filter(t_pval < 0.05) %>%
  ggplot(aes(cancer_types, y = pathway)) +
  geom_point(aes(size = -log10(t_pval), color = t_diff)) +
  scale_color_gradient2(high = "red", mid = "white", low = "blue") +
  theme(
    axis.text.x = element_text(angle = 90)
  )
