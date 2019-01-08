library(magrittr)
library(ggplot2)
# processed path

tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
#output path
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
# Read gene list
# Gene list was compress as rds
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz")) 


#######################
# filter out genes ---------------
#######################
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr

readr::write_rds(x = gene_list_expr, path = file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"), compress = "gz")


#################################
# Caculate p-value and fold-change.
##################################
calculate_fc_pvalue <- function(.x, .y) {
  .y %>%
    tibble::add_column(cancer_types = .x, .before = 1) -> df
  
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
            .f = calculate_fc_pvalue) -> gene_list_fc_pvalue
names(gene_list_fc_pvalue) <- gene_list_expr$cancer_types

gene_list_fc_pvalue %>% dplyr::bind_rows() -> gene_list_fc_pvalue_simplified

readr::write_rds(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(expr_path_a, ".rds_03_a_atg_lys_fc_pvalue_simplified.rds.gz"),
  compress = "gz"
)
readr::write_tsv(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(expr_path_a, "tsv_03_a_atg_lys_fc_pvalue_simplified.tsv")
)

# check autophagy gene expression distribution-----------------------------
gene_list_fc_pvalue_simplified %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  tidyr::drop_na(Tumor, Normal) %>% 
  dplyr::filter(status %in% c( "p")) %>%
  tidyr::gather(key = nt, value = expr, c(Tumor, Normal)) %>% 
  ggplot(aes(x = log2(expr + 1), color = nt)) +
  geom_density(aes(y = ..density..)) +
  ggthemes::scale_color_gdocs(name = "Type") +
  ggthemes::theme_gdocs() +
  facet_wrap( ~ cancer_types) +
  labs(
    x = "Expression (log2)",
    y = "Frenquncy"
  ) -> p
ggsave(filename = "03_a_expression_distribution_promotion.pdf", plot = p, device = "pdf", path = expr_path_a, width = 6, height = 6)

gene_list_fc_pvalue_simplified %>% 
  # dplyr::filter(cancer_types == "KICH") %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  tidyr::drop_na(Tumor, Normal) %>% 
  tidyr::gather(key = type, value = expr, Normal, Tumor) %>% 
  dplyr::group_by(cancer_types, pathway) %>% 
  dplyr::do(broom::tidy(t.test(expr ~ type, data = .))) %>% 
  dplyr::filter(p.value < 0.5)
  



# draw picture------------------------------------------------------------------------
gene_list_fc_pvalue_simplified %>% 
  dplyr::left_join(gene_list, by = "symbol") -> gene_fc_pvalue


filter_fc_pval <- function(.x){
  .x %>% 
    # dplyr::filter(Normal > 10 & Tumor > 10) %>%
    dplyr::filter(abs(log2(fc)) >= log2(3 / 2), fdr <= 0.05) %>%
    dplyr::mutate(p.value = -log10(p.value)) %>% 
    dplyr::mutate(p.value = ifelse(p.value > 15, 15, p.value)) %>% 
    dplyr::mutate(fdr = -log10(fdr)) %>% 
    dplyr::mutate(fdr = ifelse(fdr > 15, 15, fdr)) %>% 
    dplyr::mutate(fc = ifelse(fc < 1/8, 1/8, ifelse(fc > 8, 8, fc)))
}

# filter_pattern based on ----
filter_pattern <- function(fc, fdr) {
  if ((fc > 3/ 2) && (fdr < 0.05)) {
    return(1)
  } else if ((fc < 2/ 3) && (fdr < 0.05)) {
    return(-1)
  } else {
    return(0)
  }
}

# get up down pattern ----
get_pattern <- function(.x){
  .x %>% 
    # dplyr::filter(Normal > 10 & Tumor > 10) %>%
    dplyr::mutate(expr_pattern = purrr::map2_dbl(fc, fdr, filter_pattern)) %>%
    dplyr::select(cancer_types, symbol, expr_pattern) %>%
    tidyr::spread(key = cancer_types, value = expr_pattern) %>% 
    # dplyr::mutate_all(.funs = dplyr::funs(ifelse(is.na(.), 0, .))) 
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .)))
    
}

# gene rank by up and down ----
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

# cancer types rank by gene ----
get_cancer_types_rank <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(abs(.)))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
}

# point plot ----
plot_fc_pval_pattern <- function(.x_filter, gene_rank, cancer_types_rank){
  ggplot(.x_filter, aes(x = cancer_types, y = symbol)) +
    geom_point(aes(size = fdr, col = log2(fc))) +
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
      labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
      name = "FDR"
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


# rect plot ----
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
# cancer count ----
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

#------------------------------------

# Autophagy and Lysosome combined plot----

gene_fc_pvalue %>% filter_fc_pval() -> gene_fc_pvalue_filter
gene_fc_pvalue %>% get_pattern() -> gene_fc_pvalue_pattern

gene_fc_pvalue_pattern %>% get_cancer_types_rank() -> cancer_rank

gene_list$process %>% unique() -> process
factor(process, levels = c("ULK-complex", "PI3K-III-complex", "ATG5-ATG12-ATG16", "ATG9-cycling-system", "LC3", "TRAPP", "Membrane-delivery", "Mitophagy", "Positive", "Negative", "Lys_comp", "Lys_deg"), ordered = T) -> process
setNames(ggthemes::gdocs_pal()(12), process) -> rcb_color

gene_fc_pvalue_pattern %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(process, rcb_color)) %>% 
  dplyr::mutate(process = factor(process, levels = c("ULK-complex", "PI3K-III-complex", "ATG5-ATG12-ATG16", "ATG9-cycling-system", "LC3", "TRAPP", "Membrane-delivery", "Mitophagy", "Positive", "Negative", "Lys_comp", "Lys_deg"), ordered = T)) %>% 
  dplyr::arrange(dplyr::desc(process), rank) -> gene_rank

CPCOLS <- c("#000080", "#F8F8FF", "#CD0000")

gene_fc_pvalue_filter %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = fdr, col = log2(fc))) +
  scale_color_gradient2(
    low = CPCOLS[1],
    mid = CPCOLS[2],
    high = CPCOLS[3],
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
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
    name = "FDR"
  ) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
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
    axis.text.y = element_text(color = gene_rank$color),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  labs(subtitle = "All Autophagy & Lysosome genes differential expression in Tumor vs. Normal ") -> p

ggsave(
  filename = "fig_process_01_atg_lys.pdf",
  plot = p,
  device = "pdf",
  width = 6,
  height = 35,
  path = expr_path_a
)

# autophagy related -------------------------------------------------------

gene_fc_pvalue %>% 
  get_pattern() %>% 
  dplyr::filter(symbol %in% dplyr::pull(dplyr::filter(gene_list, status != "l"), "symbol")) -> gene_fc_pvalue_pattern

# Cancer and gene rank
gene_fc_pvalue_pattern %>% get_cancer_types_rank() -> cancer_rank

gene_fc_pvalue_pattern %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(process, rcb_color)) %>% 
  dplyr::mutate(process = factor(process, levels = c("ULK-complex", "PI3K-III-complex", "ATG5-ATG12-ATG16", "ATG9-cycling-system", "LC3", "TRAPP", "Membrane-delivery", "Mitophagy", "Positive", "Negative", "Lys_comp", "Lys_deg"), ordered = T)) %>% 
  dplyr::arrange(dplyr::desc(process), rank) -> gene_rank

CPCOLS <- c("#000080", "#F8F8FF", "#CD0000")

gene_fc_pvalue_filter %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = fdr, col = log2(fc))) +
  scale_color_gradient2(
    low = CPCOLS[1],
    mid = CPCOLS[2],
    high = CPCOLS[3],
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
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
    name = "FDR"
  ) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
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
    axis.text.y = element_text(color = gene_rank$color),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  labs(subtitle = "All Autophagy related genes differential expression in Tumor vs. Normal ") -> p

ggsave(
  filename = "fig_process_02_atg.pdf",
  plot = p,
  device = "pdf",
  width = 6,
  height = 20,
  path = expr_path_a
)

# lysosome ----------------------------------------------------------------
gene_fc_pvalue %>% 
  get_pattern() %>% 
  dplyr::filter(symbol %in% dplyr::pull(dplyr::filter(gene_list, status == "l"), "symbol")) -> gene_fc_pvalue_pattern


gene_fc_pvalue_pattern %>% get_cancer_types_rank() -> cancer_rank

gene_fc_pvalue_pattern %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(process, rcb_color)) %>% 
  dplyr::mutate(process = factor(process, levels = c("ULK-complex", "PI3K-III-complex", "ATG5-ATG12-ATG16", "ATG9-cycling-system", "LC3", "TRAPP", "Membrane-delivery", "Mitophagy", "Positive", "Negative", "Lys_comp", "Lys_deg"), ordered = T)) %>% 
  dplyr::arrange(dplyr::desc(process), rank) -> gene_rank

CPCOLS <- c("#000080", "#F8F8FF", "#CD0000")

gene_fc_pvalue_filter %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = fdr, col = log2(fc))) +
  scale_color_gradient2(
    low = CPCOLS[1],
    mid = CPCOLS[2],
    high = CPCOLS[3],
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
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
    name = "FDR"
  ) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
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
    axis.text.y = element_text(color = gene_rank$color),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  labs(subtitle = "All Lysosome related genes differential expression in Tumor vs. Normal ") -> p

ggsave(
  filename = "fig_process_03_lys.pdf",
  plot = p,
  device = "pdf",
  width = 6,
  height = 18,
  path = expr_path_a
)


# save ----
save.image(file = file.path(expr_path_a, ".rda_03_a_gene_expr.rda"))
load(file = file.path(expr_path_a, ".rda_03_a_gene_expr.rda"))




# not use any more--------------
gene_fc_pvalue_pattern %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::filter(up+down > 2) %>%
  dplyr::arrange(color, rank) -> gene_rank

p <- plot_fc_pval_pattern(gene_fc_pvalue_filter, gene_rank = gene_rank, cancer_types_rank = cancer_rank) + 
  theme(axis.text.y = element_text(color = gene_rank$color))
ggsave(
  filename = "fig_03_a_atg_lys_fc_pval_1_over_3_cancer_color.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height = 20,
  path = expr_path_a
)



# p <- plot_rect_pattern(gene_fc_pvalue_filter, gene_rank = gene_rank, cancer_types_rank = cancer_rank) +
#   theme(axis.text.y = element_text(color = gene_rank$color))
# ggsave(
#   filename = "fig_03_a_at_ly_comb_expr_rect.pdf",
#   plot = p,
#   device = "pdf",
#   width = 8,
#   height = 15,
#   path = expr_path
# )
# p <- plot_cancer_count(.x_filter = gene_fc_pvalue_filter, gene_rank = gene_rank, cancer_types_rank = cancer_rank) + 
#   theme(axis.text.y = element_text(color = gene_rank$color))
# ggsave(
#   filename = "fig_03_a_at_ly_comb_cancer_count.pdf",
#   plot = p,
#   device = "pdf",
#   width = 8,
#   height = 22,
#   path = expr_path
# )



#----------------------------------------------
# Autophagy----
#-----------------------------------------------
gene_fc_pvalue_autophagy <- 
  gene_fc_pvalue %>% 
  dplyr::filter(color %in% c("#DC3912", "#DD4477")) 
    

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
  # dplyr::filter(! symbol %in% c("TP53")) %>%
  # dplyr::filter(core == "Y") %>%
  # dplyr::filter(up + down > 2) %>%
  # dplyr::mutate(color = ifelse(is.na(core), "black", "red")) %>%
  dplyr::arrange(color, rank)

p <- plot_fc_pval_pattern(at_filter, at_gene_rank, at_cancer_rank) + 
  theme(axis.text.y = element_text(color = at_gene_rank$color))
ggsave(
  filename = "fig_03_a_atg_fc_pval_all.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height = 16,
  path = expr_path_a
)
# 
# p <- plot_rect_pattern(at_filter, at_gene_rank, at_cancer_rank) + 
#   theme(axis.text.y = element_text(color = at_gene_rank$color))
# ggsave(
#   filename = "fig_03_a_at_expr_rect.pdf",
#   plot = p,
#   device = "pdf",
#   width = 8,
#   height = 11,
#   path = expr_path
# )
# 
# 
# p <- plot_cancer_count(at_filter, at_gene_rank, at_cancer_rank) + 
#   theme(axis.text.y = element_text(color = at_gene_rank$color))
# ggsave(
#   filename = "fig_03_a_at_cancer_count.pdf",
#   plot = p,
#   device = "pdf",
#   width = 8,
#   height = 22,
#   path = expr_path
# )

#--------------------------------------------
#Lysosome----
#---------------------------------------------
gene_fc_pvalue_lysosome <- 
  gene_fc_pvalue %>% 
  dplyr::filter(status == "l")

ly_filter <- 
  gene_fc_pvalue_lysosome %>% 
  filter_fc_pval()

ly_cancer_rank <- 
  gene_fc_pvalue_lysosome %>% 
  get_pattern() %>% 
  get_cancer_types_rank() 
# dplyr::filter(!cancer_types %in% c("KICH", "LUSC"))

ly_gene_rank <- 
  gene_fc_pvalue_lysosome %>% 
  get_pattern() %>% 
  get_gene_rank() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::filter(up + down > 2) %>%
  # dplyr::mutate(color = ifelse(is.na(core), "black", "red")) %>%
  dplyr::arrange(color, rank)


p <- plot_fc_pval_pattern(ly_filter, ly_gene_rank, ly_cancer_rank) + 
  theme(axis.text.y = element_text(color = ly_gene_rank$color))
ggsave(
  filename = "new_fig_03_a_lys_fc_pval_all_color_2.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height =14,
  path = expr_path_a
)
# p <- plot_rect_pattern(ly_filter, ly_gene_rank, ly_cancer_rank) + 
#   theme(axis.text.y = element_text(color = ly_gene_rank$color))
# ggsave(
#   filename = "fig_03_a_ly_expr_rect.pdf",
#   plot = p,
#   device = "pdf",
#   width = 8,
#   height = 11,
#   path = expr_path
# )
# p <- plot_cancer_count(ly_filter, ly_gene_rank, ly_cancer_rank) + 
#   theme(axis.text.y = element_text(color = ly_gene_rank$color))
# ggsave(
#   filename = "fig_03_a_ly_cancer_count.pdf",
#   plot = p,
#   device = "pdf",
#   width = 8,
#   height = 22,
#   path = expr_path
# )



gene_rank %>% 
  dplyr::filter(up + down >=2) %>% 
  tidyr::gather(key = up_down, value = portion, c(up_p, down_p, none)) %>% 
  dplyr::mutate(color = plyr::revalue(type, replace = c("Autophagy" = "red", "Lysosome" = "black"))) %>% 
  dplyr::select(symbol, rank, up_down, portion, color) %>%
  dplyr::mutate(up_down = forcats::fct_relevel(up_down,"up_p", "down_p", "none")) %>% 
  dplyr::mutate(symbol = forcats::as_factor(symbol)) -> gene_rank_portion

gene_rank_portion %>% head(2) %>% 
  ggplot(mapping = aes(x = factor(1), y = portion, fill = up_down)) +
  geom_bar(stat = "identity",
           width = 1,
           position = "stack") +
  facet_wrap(~ symbol) +
  coord_polar("y") +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(
      size = 11,
      vjust = 0
    ),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    panel.spacing  = unit(0.02, "lines")
  ) -> p

ggsave(
  filename = "fig_03_a_at_ly_gene_portion.pdf",
  plot = p,
  device = "pdf",
  width = 20,
  height = 22,
  path = expr_path
)
