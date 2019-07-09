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
gene_list <- readr::read_tsv(file.path(gene_list_path, "ICPs_all_info_class.tsv"))

#output path
out_path<-c(file.path(result_path,"e_2_DE"))
load(file = file.path(out_path, "rda_00_gene_expr.rda"))

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
  dplyr::arrange(site_col,rank) -> gene_rank
gene_rank$site_col %>% as.character() ->gene_rank$site_col
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
    breaks = seq(-6, 6, length.out = 5),
    labels = c("-6", "-3", "0", "3", "6"),
    name = "log2 (FC)"
  ) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_x_discrete(limit = cancer_types_rank$cancer_types, expand = c(0, 0)) +
  my_theme+
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(color = gene_rank$site_col),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) -> p;p
ggsave(
  filename = "fig_01_expr_pattern-Exp_site.pdf",
  plot = p,
  device = "pdf",
  width = 10,
  height = 10,
  path = out_path
)
ggsave(
  filename = "fig_01_expr_pattern-Exp_site.png",
  plot = p,
  device = "png",
  width = 10,
  height = 10,
  path = out_path
)
readr::write_rds(
  p,
  path = file.path(out_path, "fig_01_expr_pattern-Exp_site.pdf.rds.gz"),
  compress = "gz"
)

###################
#point heat map
###################
gene_list_fc_pvalue_simplified_filter %>%
  # dplyr::mutate(fc=ifelse(fc>10,10,fc))%>%
  # dplyr::mutate(fc=ifelse(fc<0.1,0.1,fc))%>%
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(size = p.value, col = log2(fc))) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-6, 6, length.out = 5),
    #limits=c(-4,4),
    labels = c("-6", "-3", "0", "3", "6"),
    name = "Log2 (FC)"
  ) +
  scale_size_continuous(
    name="P value",
    #limit = c(-log10(0.05), 15), #set limit will delete the points out of size: -log10(0.05):15
    range = c(1, 6),
    breaks = c(-log10(0.05), 5, 10, 15),
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$"))
  ) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_x_discrete(limit = cancer_types_rank$cancer_types) +
  my_theme +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.text.y = element_text(color = gene_rank$site_col),
    axis.text.x = element_text(angle = 45,hjust = 1),
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) -> p1;p1
ggsave(
  filename = "fig_02_expr_pattern_fc_pval-Exp_site.pdf",
  plot = p1,
  device = "pdf",
  width = 8,
  height = 10,
  path = out_path
)
ggsave(
  filename = "fig_02_expr_pattern_fc_pval-Exp_site.png",
  plot = p1,
  device = "png",
  width = 8,
  height = 10,
  path = out_path
)
readr::write_rds(
  p1,
  path = file.path(out_path, "fig_02_expr_pattern_fc_pval-Exp_site.pdf.rds.gz"),
  compress = "gz"
)


##########################
#counts pic
###########################
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
  my_theme +
  theme(
    panel.background = element_rect(
      colour = "black",
      fill = "white",
      size = 1
    ),
    panel.grid.major = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8),
    axis.text.y = element_text(color = gene_rank$site_col),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  labs(y="DE gene counts") +
  coord_flip() -> p2;p2
ggsave(
  filename = "fig_03_expr_pattern_gene_counts-Exp_site.pdf",
  plot = p2,
  device = "pdf",
  width = 4,
  height = 10,
  path = out_path
)
readr::write_rds(
  p2,
  path = file.path(out_path, "fig_03_expr_pattern_cancer_counts.pdf.gz"),
  compress = "gz"
)


###cancer counts
ggplot(
  dplyr::mutate(
    gene_list_fc_pvalue_simplified_filter,
    alt = ifelse(log2(fc) > 0,  "up", "down")
  ),
  aes(x = cancer_types, fill = factor(alt))
) +
  geom_bar(color = NA, width = 0.5) +
  scale_fill_manual(
    limit = c("down", "up"),
    values = c("blue", "red"),
    guide = FALSE
  ) +
  scale_y_continuous(
    limit = c(0, 70),
    expand = c(0, 0)
    #   breaks = seq(0, 70, length.out =9)
  ) +
  scale_x_discrete(limit = cancer_types_rank$cancer_types, expand = c(0.01, 0.01)) +
  labs(y="DE gene counts") +
  my_theme +
  theme(
    panel.background = element_rect(
      colour = "black",
      fill = "white",
      size = 1
    ),
    panel.grid.major = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.text = element_text(size = 8),
    legend.key = element_rect(fill = "white", colour = "black")
  ) -> p3;p3
ggsave(
  filename = "fig_03_expr_pattern_cancer_counts-Exp_site.pdf",
  plot = p3,
  device = "pdf",
  width = 8,
  height = 3,
  path = out_path
)
readr::write_rds(
  p3,
  path = file.path(out_path, "fig_03_expr_pattern_cancer_counts-Exp_site.pdf.rds.gz"),
  compress = "gz"
)

##################
# combine point +count plot
##################
gene_list_fc_pvalue_simplified_filter %>%
  dplyr::inner_join(gene_rank,by="symbol") %>%
  dplyr::select(symbol,functionWithImmune,family) %>%
  dplyr::rename("Immunity"="functionWithImmune","Gene family"="family") %>%
  tidyr::gather(-symbol,key="fun",value="Gene class") -> gene_class_ready_draw
  
gene_class_ready_draw %>%
  dplyr::arrange(fun,`Gene class`) %>%
  .$`Gene class` -> gene_class.class.rank
gene_class_ready_draw <- within(gene_class_ready_draw,`Gene class` <- factor(`Gene class`,levels = unique(gene_class.class.rank)))
with(gene_class_ready_draw,levels(`Gene class`))
gene_class_ready_draw %>%
  ggplot(aes(y=symbol,x=fun)) +
  geom_tile(aes(fill = `Gene class`),color="grey",size=1) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_fill_manual(
    name = "Gene class",
    values = c( "#838B8B", "#000000","#0000FF", "#8B2323","#CDAA7D","#8EE5EE","#1C86EE", "#EE3B3B","#EE7600")
  ) +
  my_theme +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    axis.text.y = element_blank(),
    # axis.text.x = element_text(angle = 45,hjust = 1),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    # legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.margin=unit(c(0,0,0,-0), "cm")
  ) -> p2.1;p2.1
p1 + theme(axis.ticks.y = element_blank(),
           axis.text = element_text(color = "black"),
           plot.margin=unit(c(-0,-0,0,-0), "cm")) -> p3.1
p2 + theme(axis.ticks.y = element_blank(),
           axis.text = element_text(color = "black"),
           axis.title.x = element_text(size = 10),
           plot.margin=unit(c(0,0,0,-0), "cm")) +
  labs(y="DE gene\n counts") -> p4.1
p3 + theme(axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.text = element_text(color = "black"),
           axis.title.y = element_text(size = 10),
           plot.margin=unit(c(0,0,-0,0), "cm")) +
  labs(y="DE gene\n counts")-> p1.1
ggarrange(NULL,p1.1,NULL,p2.1,p3.1,p4.1,
          ncol = 3, nrow = 2,  align = "hv", 
          widths = c(3, 11, 5), heights = c(1, 5),
          legend = "left",
          common.legend = TRUE) -> p;p

ggsave(
  filename = "fig_04_combine_expr_pattern-Exp_site.pdf",
  device = "pdf",
  plot = p,
  width = 10,
  height = 10,
  path = out_path
)
ggsave(
  filename = "fig_04_combine_expr_pattern-Exp_site.png",
  device = "png",
  plot = p,
  width = 10,
  height = 10,
  path = out_path
)

################
#receptor and ligand heatmap
################
gene_list_fc_pvalue_simplified %>%
  dplyr::inner_join(gene_list,by="symbol") %>%
  dplyr::filter(!is.na(Recepter_pairs)) %>%
  dplyr::filter(!Recepter_pairs %in% c("Pair 10"))-> receptor_ligand.ready

receptor_ligand.ready %>%
  dplyr::arrange(Recepter_pairs,type) %>%
  dplyr::select(symbol,Recepter_pairs,size,site_col) -> receptor_ligand.rank
receptor_ligand.ready <- within(receptor_ligand.ready,symbol<-factor(symbol,levels=unique(receptor_ligand.rank$symbol)))
with(receptor_ligand.ready,levels(symbol))

receptor_ligand.rank %>%
  unique() -> receptor_ligand.color
# dplyr::inner_join(tibble::tibble(Recepter_pairs=unique(receptor_ligand.rank$Recepter_pairs),
# pair_color = c("#FFB6C1", "#FFBBFF", "#0000FF", "#8B2323", "#CDAA7D", "#8EE5EE", "#76EE00",
# "#D2691E", "#8B008B", "#6E8B3D", "#FF3030", "#006400", "#FFD700", "#EE00EE")),by="Recepter_pairs") 

receptor_ligand.ready %>%
  ggplot(aes(x = cancer_types, y = symbol, fill = log2(fc))) +
  geom_tile(color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "white",
    breaks = seq(-6, 6, length.out = 5),
    labels = c("-6", "-3", "0", "3", "6"),
    name = "log2 (FC)"
  ) +
  my_theme +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(color = "black" #c("green","green","green","green","red","blue","blue"),
                               #face =receptor_ligand.color$size#c("bold.italic","plain","bold.italic","plain","bold.italic","plain","plain")
    ),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_grid(Recepter_pairs~.,scales = "free_y",space = "free") +
  # scale_y_discrete(limit = c("TNFRSF4","TNFSF4","TNFRSF9","TNFSF9","CTLA4","CD80","CD86")) +
  scale_x_discrete(limit = cancer_types_rank$cancer_types, expand = c(0, 0)) -> p;p
ggsave(
  filename = "fig_05_ligand-receptor_heatmap.pdf",
  plot = p,
  device = "pdf",
  width = 6,
  height =9,
  path = out_path
)
ggsave(
  filename = "fig_05_ligand-receptor_heatmap.png",
  plot = p,
  device = "png",
  width = 6,
  height = 9,
  path = out_path
)

save.image(file = file.path(out_path, "rda_00_gene_expr.rda"))
