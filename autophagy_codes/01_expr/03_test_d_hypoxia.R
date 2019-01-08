library(magrittr)
library(ggplot2)

expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- file.path(expr_path, "03_a_gene_expr")
hypo_path <- "/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/"
hypo_out_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/03_f_hypoxia"

gene_list <- readr::read_rds(file.path(expr_path, "rds_03_at_ly_comb_gene_list.rds.gz"))
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

# test KICH hypoxia
fn_test_pval <- function(.d){
  if(.d$type %>% unique() %>% length() !=2 ){return(1)}
  broom::tidy(t.test(log2(expr + 1) ~ type, data = .d)) %>% .$p.value
}
fn_fc <- function(.d){
  if(.d$type %>% unique() %>% length() !=2 ){return(1)}
  .d %>% 
    dplyr::group_by(type) %>% 
    dplyr::summarise(me = mean(expr)) %>% 
    tidyr::spread(key = type, value = me) %>% 
    dplyr::mutate(fc = hypoxic / normoxic) %>% 
    .$fc
}
fn_fc_pval <- function(symbol, .d){
  print(symbol)
  pval <- fn_test_pval(.d)
  fc <- fn_fc(.d)
  tibble::tibble(pval = pval, fc = fc)
}
fn_boxplot <- function(symbol, .d){
  .d %>% 
    ggplot(aes(x = type, y = expr)) + geom_boxplot()
}
fn_hypo_classification <- function(cancer_types, filter_expr){
  print(cancer_types)
  hypo_path <- "/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/"
  hypo_file <- file.path(hypo_path,paste(cancer_types, "Hypoxia.stratification.txt", sep = "."))
  if(!file.exists(hypo_file)){
    return(tibble::tibble())
  }
  
  hypo <- readr::read_tsv(
    hypo_file, 
  col_names = c("type", "barcode"),
  skip = 1) %>%
    dplyr::mutate(barcode = stringr::str_replace_all(barcode, pattern = "\\.", "-"))
  
  filter_expr %>% 
    tidyr::gather(key = barcode, value = expr, -c(symbol, entrez_id)) %>% 
    dplyr::inner_join(hypo, by = "barcode") %>% 
    dplyr::filter(type %in% c("hypoxic", "normoxic")) %>% 
    tidyr::drop_na() %>% 
    tidyr::nest(-symbol, .key = data) %>% 
    dplyr::mutate(fc_pval = purrr::map2(symbol, data, .f = fn_fc_pval)) %>% 
    tidyr::unnest(fc_pval, .drop = F) %>% 
    tibble::add_column(cancer_types = cancer_types, .before = "symbol") %>% 
    dplyr::mutate(fdr = p.adjust(pval, method = "fdr"))
}

gene_list_expr %>% 
  dplyr::transmute(res = purrr::map2(cancer_types, filter_expr, fn_hypo_classification)) %>%
  tidyr::unnest() -> hypo_pval_fc
readr::write_rds(hypo_pval_fc, path = file.path(hypo_out_path, ".rds_03_f_hypoxia_hypo_pval_fc.rds.gz"), compress = "gz")

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
    dplyr::mutate(expr_pattern = purrr::map2_dbl(fc, fdr, filter_pattern)) %>%
    dplyr::select(cancer_types, symbol, expr_pattern) %>%
    tidyr::spread(key = cancer_types, value = expr_pattern) %>% 
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .)))
  
}
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
    dplyr::mutate(up_p = up / 14, down_p = down / 14, none = 1 - up_p - down_p) %>% 
    dplyr::arrange(rank)
}
get_cancer_types_rank <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(abs(.)))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
}


hypo_pval_fc %>%
  dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.nan(.), 0, .))) %>%  get_pattern() -> pattern

pattern %>% get_cancer_types_rank() -> cancer_types_rank
pattern %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  # dplyr::filter(up+down > 2) %>%
  dplyr::mutate(color = plyr::revalue(type, replace = c("Lysosome" = "black", "Autophagy" = "red"))) %>% 
  dplyr::arrange(color, rank) -> gene_rank


hypo_pval_fc %>% 
  dplyr::filter(fdr < 0.05, abs(log2(fc)) > log2(2)) %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = -log10(fdr), col = log2(fc))) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-3, 3, length.out = 5),
    labels = c("<= -3", "-1.5", "0", "1.5", ">= 3"),
    name = "log2 FC"
  ) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_x_discrete(limit = cancer_types_rank$cancer_types) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  theme(axis.text.y = element_text(color = gene_rank$color)) -> p

ggsave(
  filename = "fig_03_f_hypoxia_all.pdf",
  plot = p,
  device = "pdf",
  width = 15,
  height = 30,
  path = hypo_out_path
)


save.image(file = file.path(hypo_out_path, ".rda_03_f_hypoxia.rda"))
load(file = file.path(hypo_out_path, ".rda_03_f_hypoxia.rda"))









