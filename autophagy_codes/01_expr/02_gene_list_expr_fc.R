library(magrittr)
library(dplyr)
library(RDS)
library(gtools)
# processed path
tcga_path = c("/project/huff/huff/immune_checkpoint/data/TCGA_data")
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
result_path<-c("/project/huff/huff/immune_checkpoint/result_20171025")
expr_path<-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")
# Read gene list
# Gene list was compress as rds
gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

#output path
out_path<-c(file.path(result_path,"e_2_DE"))

#######################
# filter out genes
#######################
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}
expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr
gene_list_expr %>%
  readr::write_rds(file.path(expr_path,".rds_03_a_gene_list_expr.rds.gz"), compress = "gz")

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
  path = file.path(out_path, "rds_01_gene_list_fc_pvalue_simplified.rds.gz"),
  compress = "gz"
)
readr::write_tsv(
  x = gene_list_fc_pvalue_simplified,
  path = file.path(out_path, "tsv_01_gene_list_fc_pvalue_simplified.tsv")
)

# write sample pairs number
gene_list_fc_pvalue_simplified %>%
  dplyr::select(cancer_types, n_normal, n_tumor) %>%
  dplyr::distinct() -> pancan_samples_pairs
readr::write_rds(
  x = pancan_samples_pairs,
  path = file.path(out_path, "rds_02_pancan_samples_pairs.rds.gz"),
  compress = "gz"
)
readr::write_tsv(x = pancan_samples_pairs,
                 path = file.path(out_path, "tsv_02_pancan_samples_pairs.tsv"))

###############
#Draw pictures
################
gene_list_fc_pvalue_simplified %>%
  # filter |log2(fc)| >= log2(1.5), fdr <= 0.05
  dplyr::filter(abs(log2(fc)) >= log2(1.5), fdr <= 0.05) %>%
  dplyr::mutate(p.value = -log10(p.value)) -> gene_list_fc_pvalue_simplified_filter

# expression pattern
# significant high expression is 1
# significant low expression is -1
# not significant is 0
get_pattern <- function(fc, p.value) {
  if ((fc > 1.5) && (p.value < 0.05)) {
    return(1)
  } else if ((fc < 2 / 3) && (p.value < 0.05)) {
    return(-1)
  } else {
    return(0)
  }
}

gene_list_fc_pvalue_simplified %>%
  dplyr::mutate(expr_pattern = purrr::map2_dbl(fc, p.value, get_pattern)) %>%
  dplyr::select(cancer_types, symbol, expr_pattern) %>%
  tidyr::spread(key = cancer_types, value = expr_pattern) ->
  gene_expr_pattern

gene_expr_pattern %>%
  dplyr::rowwise() %>%
  dplyr::do(
    symbol = .$symbol,
    rank =  unlist(.[-1], use.names = F) %>% sum(),
    up = (unlist(.[-1], use.names = F) == 1) %>% sum(),
    down = (unlist(.[-1], use.names = F) == -1) %>% sum()
  ) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() %>%
  dplyr::left_join(gene_list, by = "symbol") %>%
  dplyr::arrange(col,rank) -> gene_rank
gene_rank$col %>% as.character() ->gene_rank$col
gene_rank$size %>% as.character() ->gene_rank$size

gene_expr_pattern %>%
  dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(.))) %>%
  tidyr::gather(key = cancer_types, value = rank) %>%
  dplyr::arrange(-rank) -> cancer_types_rank

library(ggplot2)
ggplot(gene_list_fc_pvalue_simplified_filter,
       aes(x = cancer_types, y = symbol, fill = log2(fc))) +
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
    axis.text.y = element_text(color = gene_rank$col),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) -> p;p
ggsave(
  filename = "fig_01_expr_pattern.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 12,
  path = out_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_01_expr_pattern.pdf.rds.gz"),
  compress = "gz"
)


ggplot(gene_list_fc_pvalue_simplified_filter,
       aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = p.value, col = log2(fc))) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-3, 3, length.out = 5),
    labels = c("<= -3", "-1.5", "0", "1.5", ">= 3"),
    name = "Fold change"
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
    axis.text.y = element_text(color = gene_rank$col),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) -> p;p
ggsave(
  filename = "fig_02_expr_pattern_fc_pval.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 12,
  path = out_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_02_expr_pattern_fc_pval.pdf.rds.gz"),
  compress = "gz"
)



ggplot(
  dplyr::mutate(
    gene_list_fc_pvalue_simplified_filter,
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
    axis.text.y = element_text(color = gene_rank$col),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  coord_flip() -> p;p
ggsave(
  filename = "fig_03_expr_pattern_cancer_counts.pdf",
  plot = p,
  device = "pdf",
  width = 4,
  height = 12,
  path = out_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_03_expr_pattern_cancer_counts.pdf.gz"),
  compress = "gz"
)

save.image(file = file.path(out_path, "rda_00_gene_expr.rda"))
load(file = file.path(out_path, "rda_00_gene_expr.rda"))
