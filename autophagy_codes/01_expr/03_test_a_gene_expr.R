library(magrittr)
library(ggplot2)

# split autophagy gene and lysosome genes
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr"
load(file = file.path(expr_path, "rda_00_gene_expr.rda"))
rm(expr)

# gene_list
# gene_list_fc_pvalue_simplified
# gene_list_fc_pvalue_simplified_filter
# 
# gene_expr_pattern
# 
# gene_rank
# cancer_types_rank

gene_list_fc_pvalue_simplified %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::select(-c(desc, ensembl_gene_id)) -> gene_fc_pvalue

gene_fc_pvalue_autophagy <- 
  gene_fc_pvalue %>% 
  dplyr::filter(type == "Autophagy")

gene_fc_pvalue_lysosome <-
  gene_fc_pvalue %>% 
  dplyr::filter(type == "Lysosome")

# filter fc and pval
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

plot_rect_point_count <- function(.x){
  # .x is gene_list_fc_pvalue_simplified
  # pattern
  type = unique(.x$type)
  
  # filter out not significant
  .x_filter <- .x %>% filter_fc_pval()
    
  # get pattern
  .x %>% get_pattern() -> pattern
  
  # gene rank
  pattern %>% get_gene_rank() -> gene_rank
  
  # cancer rank
  pattern %>% get_cancer_types_rank() -> cancer_types_rank
  
  #plot_rect_pattern
  p_rect <- plot_rect_pattern(.x_filter, gene_rank = gene_rank, cancer_types_rank = cancer_types_rank )
  #plot_fc_pval_pattern
  p_point <- plot_fc_pval_pattern(.x_filter, gene_rank = gene_rank, cancer_types_rank = cancer_types_rank)
  #plot_cancer_count
  p_count <- plot_cancer_count(.x_filter, gene_rank = gene_rank, cancer_types_rank = cancer_types_rank)
  
  tibble::tibble(type = type, p_rect = list(p_rect), p_point = list(p_point), p_count = list(p_count))
}
plot_rect_point_count(gene_fc_pvalue_autophagy) -> autophagy_plot
plot_rect_point_count(gene_fc_pvalue_lysosome) -> lysosome_plot
#################################################################

# autophagy
at_filter <- 
  gene_fc_pvalue_autophagy %>% 
  filter_fc_pval()

at_gene_rank <- 
  gene_fc_pvalue_autophagy %>% 
  get_pattern() %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(desc, replace = c("Initation" = "red", "Nucleation" = "green", "Elongation" = "blue", "Fusion" = "yellow"))) %>% 
  dplyr::mutate(color = ifelse(desc %in% c("Initation", "Nucleation", "Fusion", "Elongation"), color, "black")) %>% 
  dplyr::arrange(desc(complex), rank)

at_cancer_rank <- 
  gene_fc_pvalue_autophagy %>% 
  get_pattern() %>% 
  get_cancer_types_rank()

p <- plot_rect_pattern(at_filter, at_gene_rank, at_cancer_rank) + 
  theme(axis.text.y = element_text(color = at_gene_rank$color))
ggsave(
  filename = "fig_04_at_expr_rect.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 20,
  path = expr_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_04_at_expr_rect.pdf.rds.gz"),
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
  path = file.path(out_path, "fig_05_at_expr_fc_pval.pdf.rds.gz"),
  compress = "gz"
)


# autophagy
ly_filter <- 
  gene_fc_pvalue_lysosome %>% 
  filter_fc_pval()

ly_gene_rank <- 
  gene_fc_pvalue_lysosome %>% 
  get_pattern() %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(desc, replace = c("Initation" = "red", "Nucleation" = "green", "Elongation" = "blue", "Fusion" = "yellow"))) %>% 
  dplyr::mutate(color = ifelse(desc %in% c("Initation", "Nucleation", "Fusion", "Elongation"), color, "black")) %>% 
  dplyr::arrange(desc(complex), rank)

ly_cancer_rank <- 
  gene_fc_pvalue_lysosome %>% 
  get_pattern() %>% 
  get_cancer_types_rank()

p <- plot_rect_pattern(ly_filter, ly_gene_rank, ly_cancer_rank) + 
  theme(axis.text.y = element_text(color = ly_gene_rank$color))
ggsave(
  filename = "fig_04_ly_expr_rect.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 20,
  path = expr_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_04_ly_expr_rect.pdf.rds.gz"),
  compress = "gz"
)

p <- plot_fc_pval_pattern(ly_filter, ly_gene_rank, ly_cancer_rank) + 
  theme(axis.text.y = element_text(color = ly_gene_rank$color))
ggsave(
  filename = "fig_05_ly_expr_fc_pval.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 20,
  path = expr_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_05_ly_expr_fc_pval.pdf.rds.gz"),
  compress = "gz"
)



save.image(file = file.path(expr_path, "rda_01_gene_expr_detail.rda.gz"), compress = T)
