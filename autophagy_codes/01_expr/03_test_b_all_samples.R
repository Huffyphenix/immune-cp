library(magrittr)
library(ggplot2)
# processed path
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr <- readr::read_rds(file.path(tcga_path, "pancan_expr_20160513.rds.gz"))


#output path
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"

# Read gene list
# Gene list was compress as rds

gene_list <- readr::read_rds(file.path(expr_path, "rds_03_atg_gene_list.rds.gz"))

expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr


filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

calculate_fc_pvalue_all_samples <- function(.x, .y) {
  print(.x)
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
    # dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
    dplyr::ungroup()
  
  sample_type_summary <- table(samples$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) || length(sample_type_summary) != 2 ||
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
    dplyr::group_by(cancer_types, symbol, entrez_id) %>%
    tidyr::drop_na(expr) %>%
    dplyr::do(broom::tidy(t.test(expr ~ type, data = .))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    dplyr::select(cancer_types, symbol, entrez_id, p.value, fdr) -> df_pvalue
  
  # log2 fold change mean
  df_f %>%
    dplyr::group_by(cancer_types, symbol, entrez_id, type) %>%
    tidyr::drop_na(expr) %>%
    dplyr::summarise(mean = mean(expr)) %>%
    tidyr::spread(key = type, mean) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fc = (Tumor + 0.1) / (Normal + 0.1)) -> df_fc
  
  df_fc %>%
    dplyr::inner_join(df_pvalue, by = c("cancer_types", "symbol", "entrez_id")) %>%
    dplyr::mutate(n_normal = sample_type_summary[1], n_tumor = sample_type_summary[2]) -> res
  return(res)
}

purrr::map2(.x = gene_list_expr$cancer_types,
            .y = gene_list_expr$filter_expr,
            .f = calculate_fc_pvalue_all_samples) -> gene_list_fc_pvalue
names(gene_list_fc_pvalue) <- gene_list_expr$cancer_types

gene_list_fc_pvalue %>% dplyr::bind_rows() -> gene_list_fc_pvalue_simplified

expr_path <- file.path(expr_path, "03_test_b_all_samples")

readr::write_rds(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(expr_path, "rds_03_a_at_gene_list_fc_pvalue_simplified.rds.gz"),
  compress = "gz"
)
readr::write_tsv(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(expr_path, "rds_03_a_at_gene_list_fc_pvalue_simplified.tsv")
)



gene_list_fc_pvalue_simplified %>% 
  dplyr::left_join(gene_list, by = "symbol") -> gene_fc_pvalue

gene_fc_pvalue_autophagy <- 
  gene_fc_pvalue %>% 
  dplyr::filter(type == "Autophagy")

filter_fc_pval <- function(.x){
  .x %>% 
    dplyr::filter(abs(log2(fc)) >= log2(1.5), fdr <= 0.05) %>%
    dplyr::mutate(p.value = -log10(p.value)) %>% 
    dplyr::mutate(p.value = ifelse(p.value > 15, 15, p.value)) %>% 
    dplyr::mutate(fc = ifelse(fc < 1/8, 1/8, ifelse(fc > 8, 8, fc)))
}

# filter_pattern based on
filter_pattern <- function(fc, p.value) {
  if ((fc > 1.5) && (p.value < 0.05)) {
    return(1)
  } else if ((fc < 2 / 3) && (p.value < 0.05)) {
    return(-1)
  } else {
    return(0)
  }
}

# get up down pattern
get_pattern <- function(.x){
  .x %>% 
    dplyr::mutate(expr_pattern = purrr::map2_dbl(fc, p.value, filter_pattern)) %>%
    dplyr::select(cancer_types, symbol, expr_pattern) %>%
    tidyr::spread(key = cancer_types, value = expr_pattern)
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
    dplyr::arrange(rank) 
}

# cancer types rank by gene
get_cancer_types_rank <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(abs(.)))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(-rank)
}

# rect plot
plot_rect_pattern <- function(.x_filter, gene_rank, cancer_types_rank){
  ggplot(.x_filter, aes(x = cancer_types, y = symbol, fill = log2(fc))) +
    geom_tile(color = "black") +
    scale_fill_gradient2(
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
    scale_x_discrete(limit = cancer_types_rank$cancer_types, expand = c(0, 0)) +
    theme(
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.key = element_rect(fill = "white", colour = "black")
    ) -> p
  print(p)
  return(p)
}

# point plot
plot_fc_pval_pattern <- function(.x_filter, gene_rank, cancer_types_rank){
  ggplot(.x_filter, aes(x = cancer_types, y = symbol)) +
    geom_point(aes(size = p.value, col = log2(fc))) +
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
    scale_size_continuous(
      limit = c(-log10(0.05), 15),
      range = c(1, 6),
      breaks = c(-log10(0.05), 5, 10, 15),
      labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$"))
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
    ) -> p
  print(p)
  return(p)
}

# cancer count
plot_cancer_count <- function(.x_filter, gene_rank, cancer_types_rank){
  ggplot(
    dplyr::mutate(
      .x_filter,
      alt = ifelse(log2(fc) > 0,  "up", "down")
    ),
    aes(x = symbol, fill = factor(alt))
  ) +
    geom_bar(color = NA, width = 0.5) +
    scale_fill_manual(
      limit = c("down", "up"),
      values = c("blue", "red"),
      guide = FALSE
    ) +
    scale_y_continuous(
      limit = c(-0.1, 12.5),
      expand = c(0, 0),
      breaks = seq(0, 12, length.out = 5)
    ) +
    scale_x_discrete(limit = gene_rank$symbol, expand = c(0.01, 0.01)) +
    theme(
      panel.background = element_rect(
        colour = "black",
        fill = "white",
        size = 1
      ),
      panel.grid.major = element_line(linetype = "dashed", color = "lightgray"),
      axis.title = element_blank(),
      axis.ticks.x = element_blank(),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.key = element_rect(fill = "white", colour = "black")
    ) +
    coord_flip() -> p
  print(p)
  return(p)
}

at_filter <- 
  gene_fc_pvalue_autophagy %>% 
  filter_fc_pval()

at_cancer_rank <- 
  gene_fc_pvalue_autophagy %>% 
  get_pattern() %>% 
  get_cancer_types_rank() 
# dplyr::filter(!cancer_types %in% c("KICH", "LUSC"))

at_gene_rank <- 
  gene_fc_pvalue_autophagy %>% 
  get_pattern() %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::filter(up + down >= 2) %>% 
  dplyr::arrange(rank)

p <- plot_rect_pattern(at_filter, at_gene_rank, at_cancer_rank) + 
  theme(axis.text.y = element_text(color = at_gene_rank$color))
ggsave(
  filename = "fig_04_at_core_expr_rect.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 20,
  path = expr_path
)
readr::write_rds(
  p,
  path = file.path(expr_path, "fig_04_at_core_expr_rect.pdf.rds.gz"),
  compress = "gz"
)


p <- plot_fc_pval_pattern(at_filter, at_gene_rank, at_cancer_rank) + 
  theme(axis.text.y = element_text(color = at_gene_rank$color))
ggsave(
  filename = "fig_05_at_expr_fc_pval.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 20,
  path = expr_path
)
readr::write_rds(
  p,
  path = file.path(expr_path, "fig_05_at_expr_fc_pval.pdf.rds.gz"),
  compress = "gz"
)

