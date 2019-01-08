library(ggplot2)
library(magrittr)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
#output path
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
hypo_out_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/03_f_hypoxia"
hypo_info <- readr::read_rds(path = "/extraspace/yye1/analysis/Hypoxia/PSM/summary/All.signatures.summary.rds") %>% dplyr::as_tibble()
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

hypo_mrna <- 
  hypo_info %>% 
  dplyr::filter(signatures == "mRNA") %>% 
  dplyr::rename(symbol = feature.sig, 
                fdr = fdr.sig, 
                normaxia = mean0.sig, 
                hypoxia =  mean1.sig) 

filter_pattern <- function(fc, fdr) {
  if ((fc > 2) && (fdr < 0.05)) {
    return(1)
  } else if ((fc < 1/2) && (fdr < 0.05)) {
    return(-1)
  } else {
    return(0)
  }
}


get_pattern <- function(.x){
  .x %>% 
    # dplyr::filter(Normal > 10 & Tumor > 10) %>%
    dplyr::mutate(expr_pattern = purrr::map2_dbl(fc, fdr, filter_pattern)) %>%
    dplyr::select(cancer_types, symbol, expr_pattern) %>%
    tidyr::spread(key = cancer_types, value = expr_pattern) %>% 
    # dplyr::mutate_all(.funs = dplyr::funs(ifelse(is.na(.), 0, .))) 
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .)))
  
}

# gene rank by up and down
get_gene_rank <- function(pattern){
  pattern %>% 
    dplyr::rowwise() %>%
    dplyr::do(
      symbol = .$symbol,
      rank =  unlist(.[-1], use.names = F) %>% sum(),
      up = (unlist(.[-1], use.names = F) == 1) %>% sum(),
      down = (unlist(.[-1], use.names = F) == -1) %>% sum()
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest() %>%
    dplyr::mutate(up_p = up / 15, down_p = down / 15, none = 1 - up_p - down_p) %>% 
    dplyr::arrange(rank)
}

# cancer types rank by gene
get_cancer_types_rank <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(.))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
}



hypo_mrna %>% 
  dplyr::inner_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(fc = hypoxia / normaxia) %>% 
  dplyr::select(symbol, fdr, normaxia, hypoxia, fc, pathway, color,cancer_types =  class)  -> hypo_ready

hypo_ready %>% get_pattern() -> pattern
pattern %>% get_cancer_types_rank() -> cancer_rank
pattern %>% get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(status, replace = c('a' = "#e41a1c", "l" = "#377eb8", "i" = "#4daf4a", "p" = "#984ea3"))) %>% 
  # dplyr::filter(status != "l") %>% 
  dplyr::arrange( rank) -> gene_rank

hypo_ready %>% 
  dplyr::mutate(fc = ifelse(fc > 4, 4, fc)) %>% 
  dplyr::mutate(fc = ifelse(fc < 1/ 4, 1/4, fc)) %>% 
  dplyr::mutate(fdr = -log10(fdr)) %>% 
  dplyr::mutate(fdr = ifelse(fdr > 15, 15, fdr)) %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = fdr, color = log2(fc))) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(name = "FDR") +
  scale_color_gradient2(name = "FC", low = "blue", mid = "white", high = "red") +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.y = element_text(color = gene_rank$color)
  ) -> p
ggsave(filename = "hypoxia_autophagy2.pdf", plot = p, device = "pdf", path = hypo_out_path, height = 10, width = 8)
