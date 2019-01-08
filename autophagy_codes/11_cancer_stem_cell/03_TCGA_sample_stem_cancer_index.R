library(magrittr)
tcga_dir <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
pcbc_dir <- "/home/cliu18/liucj/projects/6.autophagy/synapse/PCBC"
csc_dir <- "/home/cliu18/liucj/projects/6.autophagy/08_cancer_stem_cell"

expr <- readr::read_rds(path = file.path(tcga_dir, "pancan33_expr.rds.gz"))
res_model <- readr::read_rds(file.path(pcbc_dir, "02_OCLR_res_model.rds.gz")) 

gene_order <- res_model$gene_order %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(symbol = value)

filter_gene_list <- function(.x, gene_list) {.x
  
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol") %>% 
    dplyr::filter(entrez_id != "728661")
}

expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_order)) %>%
  dplyr::select(-expr) -> expr_sci_rd

fn_sci <- function(.x, .model = res_model){
  .x %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    tidyr::spread(key = symbol, value = expr) -> .d
  
  setNames(as.list(rep(0, length(.model$gene_order))), .model$gene_order) -> .replace
  
  .d %>% 
    dplyr::select(.model$gene_order) %>% 
    tidyr::replace_na(replace = .replace) %>% 
    as.matrix() -> .d_mat

  log2(.d_mat + 0.01) %>% scale(scale = F) -> .mat
  .mat[is.nan(.mat)] <- 0
  
  
  .score <- .mat  %*% .model$model$w
  
  hist(.score)
  
  # spearman correlation
  # .score <- apply(.mat, 1, function(.m){cor(.m, .model$model$w, method = "spearman")})
  # .score %>% hist()
  
  tibble::tibble(barcode = .d$barcode, score = as.vector(.score))
}

expr_sci_rd %>% 
  dplyr::mutate(sci = purrr::map(.x = filter_expr, .f = fn_sci, .model = res_model)) %>% 
  dplyr::select(-filter_expr) -> 
  expr_sci_score

readr::write_rds(expr_sci_score, path = file.path(csc_dir, ".rds_03_sample_score_expr_sci_score.rds.gz"), compress = "gz")

expr_sci_score %>% 
  dplyr::mutate(percent = purrr::map_dbl(sci, .f = function(.s){sum(.s$score > 0)/ length(.s$score)}) ) %>% 
  print(n = Inf)
#
save.image(file = file.path(csc_dir, ".rda_03_TCGA_sample_stem_cancer_index.rda"))

