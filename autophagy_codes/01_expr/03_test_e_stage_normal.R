library(magrittr)
library(ggplot2)

load(file = file.path(stage_path, ".rda_03_b_stage_gene_expr.rda"))
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
stage_path <- file.path(expr_path, "03_b_stage")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- file.path(expr_path, "03_a_gene_expr")
gene_list <- readr::read_rds(file.path(expr_path, "rds_03_a_atg_lys_gene_list.rds.gz"))
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))


fn_paired_sample <- function(.x, .y, .z){
  .te$cancer_types -> .x
  .te$filter_expr[[1]] -> .y
  .te$stage[[1]] -> .z
  .y %>% tibble::add_column(cancer_types = .x, .before = 1) -> df
  
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
  
  df %>% 
    tidyr::gather(key = barcode, value = expr, -c(cancer_types, symbol, entrez_id)) %>% 
    dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 12)) %>% 
    dplyr::mutate(type = stringr::str_sub(barcode, start = 14, end = 15)) %>% 
    dplyr::filter(type %in% c("01", "11")) %>% 
    dplyr::mutate(type = dplyr::recode(type, "01" = "Tumor", "11" = "Normal")) %>% 
    dplyr::inner_join(dplyr::rename(.z, sample = barcode), by = "sample") %>% 
    dplyr::filter(symbol == "BECN1", type == "Tumor") -> .te
  
  .te %>% ggplot(aes(x = stage, y = log2(expr + 1))) + geom_boxplot()
  
  df_f <-
    df %>%
    dplyr::select(c(1, 2, 3), samples$barcode) %>%
    tidyr::gather(key = barcode, value = expr, -c(1, 2, 3)) %>%
    dplyr::left_join(samples, by = "barcode") 
  
  .z %>% 
    dplyr::rename(sample = barcode) %>% 
    dplyr::inner_join(df_f, by = "sample") %>% 
    dplyr::select(sample, stage, cancer_types, symbol, expr, type) -> stage_ready
  
  stage_ready %>% 
    dplyr::group_by(stage, symbol, type) %>% 
    dplyr::count() %>% 
    dplyr::filter(symbol == "AGA")
  
  
  stage_ready %>%
    dplyr::filter(stage != "Stage IV") %>% 
    dplyr::group_by(stage, symbol,  type) %>%
    tidyr::drop_na(expr) %>%
    dplyr::summarise(mean = mean(expr)) %>%
    tidyr::spread(key = type, mean) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fc = (Tumor + 0.1) / (Normal + 0.1)) -> stage_ready_fc
  
  stage_ready %>% 
    dplyr::filter(stage != "Stage IV") %>% 
    dplyr::group_by(stage, symbol) %>% 
    tidyr::drop_na(expr) %>% 
    dplyr::do(broom::tidy(t.test(log2(expr + 1) ~ type, data = .))) -> stage_ready_pval
  
  stage_ready_pval %>% 
    dplyr::select(stage, symbol, p.value) %>% 
    dplyr::group_by(stage) %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::ungroup() %>% 
    dplyr::inner_join(stage_ready_fc, by = c("stage", "symbol")) %>% 
    dplyr::filter(fdr < 0.05, abs(log2(fc)) > log2(1.5)) %>% 
    dplyr::mutate(fc = ifelse(fc < 1/4, 1/4, ifelse(fc > 4, 4, fc))) -> draw_ready

  stage_ready %>%
    dplyr::filter(stage != "Stage IV") %>% 
    tidyr::drop_na(expr) %>%
    dplyr::mutate(expr = log2(expr + 1)) %>% 
    dplyr::filter(symbol== "TFE3") %>% 
    ggplot(aes(x = type, y = expr)) +
    geom_boxplot() +
    facet_grid(symbol~stage) -> p
  ggsave(filename = "test_stage_boxplot.pdf", plot = p, device = "pdf", path = stage_path, width = 10, height = 45)  

  gene_list %>% 
    dplyr::filter(symbol %in% unique(draw_ready$symbol)) -> gene_rank
    
  draw_ready %>% 
    ggplot(aes(x = stage, y = symbol)) +
    geom_point(aes(color = log(fc), size = -log10(fdr))) +
    scale_size_continuous(name = "FDR") +
    scale_color_gradient2(name = "FC", low="blue", mid = "white", high = "red") +
    scale_y_discrete(limit = gene_rank$symbol) +
    ggthemes::theme_gdocs()+
    theme(axis.text.y = element_text(color = gene_rank$color)) -> p
  ggsave(filename = "test_brca_stage_point.pdf", plot = p, device = "pdf", path = stage_path, width = 5, height = 20)  
  
}

expr_stage %>% 
  dplyr::filter(cancer_types == "BRCA") -> .te
  dplyr::mutate(merged_clean = purrr::pmap(cancer_types, filter_expr, stage, fn_paired_sample)) 
