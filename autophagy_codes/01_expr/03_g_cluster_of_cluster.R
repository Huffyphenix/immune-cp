library(magrittr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
# processed path
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
pancan_color <- readr::read_tsv(file.path(tcga_path, 'PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv')) 
pcc <- pancan_color %>% dplyr::pull(color)
names(pcc) <- pancan_color %>% dplyr::pull(cancer_types)

#output path
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
# Read gene list
# Gene list was compress as rds
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

# readr::write_rds(x = gene_list_expr, path = file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"), compress = "gz")

gene_list_expr <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"))

marker_file <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_marker.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz")) %>% 
  dplyr::left_join(marker_file, by = "symbol") %>% 
  dplyr::mutate(symbol = dplyr::recode(symbol, "ATG101" = "C12orf44"))


atg_gene <- gene_list %>% dplyr::filter(type == "atg") %>%  dplyr::pull(symbol)
lys_gene <- gene_list %>% dplyr::filter(type == "lys") %>% dplyr::pull(symbol)

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}

fn_filter_tumor <- function(filter_expr, genes){
  # filter_expr <- te$filter_expr[[1]]
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::drop_na() %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type != 11) %>% 
    dplyr::select(-type) %>% 
    tidyr::spread(key = barcode, value = expr) %>% 
    dplyr::filter(symbol %in% genes)
}

gene_list_expr %>% 
  # head(1) %>% 
  dplyr::mutate(filter_expr = purrr::map(.x = filter_expr, .f = fn_filter_tumor, genes = atg_gene)) %>% 
  tibble::deframe() -> gene_list_expr_tumor_atg


gene_list_expr_tumor_atg %>% 
  dplyr::bind_cols() %>%
  dplyr::select(-dplyr::matches('symbol\\d+')) %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames(var = "symbol") %>% 
  data.matrix() -> atg_mat


# atg clustering and complexheatmap --------------------
# try nmf methods
samples_mat <- atg_mat %>% t() 

res <- factoextra::hcut(samples_mat, k = 6, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T)

samples_subgroup <- cutree(res, k = 6)

atg_mat_scale <- apply(atg_mat, 1, scale) %>% t()

purrr::map_dbl(.x = gene_list_expr_tumor_atg, .f = function(x){ncol(x) - 1}) -> atg_mat_num

cancer <- rep(names(atg_mat_num), times = atg_mat_num)
cluster_col <- c("#00008B", "#00FF00", "#E31A1C", "#8B1C62", "#00F5FF", "#912CEE") #RColorBrewer::brewer.pal(name = "Set1", n = 6)
names(cluster_col) <- c(1,2,3,4,5,6)

# complexheatmap ----------------------------------------------
ha = HeatmapAnnotation(
  df = data.frame(cluster = samples_subgroup[colnames(atg_mat)], cancer = cancer), 
  gap = unit(c(4,2), "mm"),
  col = list(
    cancer = pcc, 
    cluster = cluster_col)
  )
# dev.off()
# draw(ha, 1: 9744)

Heatmap(
  atg_mat_scale, 
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"), space = "RGB"),
  name = "expression",
  
  # annotation
  top_annotation = ha, 
  
  # show row and columns
  show_row_names = T, 
  show_column_names = FALSE, 
  show_row_dend = F,
  show_column_dend = F, 
  
  row_names_gp = gpar(fontsize = 6),
  # clustering
  clustering_distance_columns = "pearson",
  clustering_method_columns = "ward.D",
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D"
  ) -> ht
ht
pdf(file.path(expr_path_a, "cluster_atg_mrna_01.pdf"), width=10, height = 10)
draw(ht, show_heatmap_legend = F, show_annotation_legend = F)
decorate_annotation(
  "cancer", 
  {grid.text("cancer", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 12))}
  )
decorate_annotation(
  "cluster", 
  {grid.text("cluster", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 12))}
  )
dev.off()

# because of the old version of ComplexHeatmap
# save pcc to local and draw annotation legend by local computer

pcc %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_g_cluster_pcc.rds.gz"), compress = 'gz')

# subgroup-------------------------------------
fn_subtype <- function(value, samples_subgroup){
  # value <- te$value[[1]]
  names(value)[-1] -> samples
  samples_subgroup[samples] %>% 
    tibble::enframe() %>% 
    dplyr::group_by(value) %>% 
    dplyr::count() %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(cluster = value, count = n)
}

gene_list_expr_tumor_atg %>% 
  tibble::enframe() %>% 
  dplyr::mutate(cluster_count = purrr::map(.x = value, .f = fn_subtype, samples_subgroup = samples_subgroup)) %>% 
  dplyr::select(-value) %>% 
  tidyr::unnest() -> cancer_cluster_freq

cancer_cluster_freq %>% 
  dplyr::group_by(name) %>% 
  dplyr::mutate(freq = count / sum(count)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = as.factor(cluster), y = name, fill = freq)) +
  geom_tile() +
  geom_text(aes(label = count)) +
  scale_fill_gradient(
    name = "Percent (%)",
    limit = c(0, 1),
    high = "red",
    low = "white",
    na.value = "white"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      color = 'black'
    ),
    axis.ticks = element_blank(),
    
    panel.background = element_blank(),
    panel.border = element_rect(
      color = "black", 
      size = 0.2, 
      fill = NA
      )
  ) +
  guides(
    fill = guide_legend(
      title = "Percent (%)",
      title.position = 'left',
      title.theme = element_text(angle = 90, vjust = 2, size = 10),
      reverse = T,
      keywidth = 0.6,
      keyheight = 0.7
    )
  ) -> tile_plot

ggsave(filename = "cluster_atg_subtype.pdf", device = "pdf", plot = tile_plot, path = expr_path_a, width = 5, height = 6)


cluster_col %>% 
  tibble::enframe() %>% 
  ggplot(aes(x = as.factor(name), y = 1, fill = name)) +
  geom_tile() +
  scale_x_discrete(
    labels = c("C1", "C2", "C3", "C4", "C5", "C6"),
    position = "top") +
  scale_fill_manual(
    values = unname(cluster_col),
    guide = F) +
  theme(
    axis.title = element_blank(),
    # axis.text = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      size = 16,
      face = "bold"
    ),
    axis.ticks = element_blank(),
    
    panel.background = element_blank(),
    panel.border = element_blank(),
    
    plot.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "mm")
    
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 0.3) -> anno_plot
ggsave(filename = "cluster_atg_subtype_anno.pdf", device = "pdf", plot = anno_plot, path = expr_path_a, width = 4, height = 1)

#----------------------------------------------
# lysosome------------------------------
gene_list_expr %>% 
  # head(1) %>% 
  dplyr::mutate(filter_expr = purrr::map(.x = filter_expr, .f = fn_filter_tumor, genes = lys_gene)) %>% 
  tibble::deframe() -> gene_list_expr_tumor_lys

gene_list_expr_tumor_lys %>% 
  dplyr::bind_cols() %>%
  dplyr::select(-dplyr::matches('symbol\\d+')) %>% 
  as.data.frame() %>% 
  tibble::column_to_rownames(var = "symbol") %>% 
  data.matrix() -> lys_mat

samples_mat <- lys_mat %>% t() 

res <- factoextra::hcut(samples_mat, k = 6, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T)

samples_subgroup <- cutree(res, k = 6)

lys_mat_scale <- apply(lys_mat, 1, scale) %>% t()

purrr::map_dbl(.x = gene_list_expr_tumor_lys, .f = function(x){ncol(x) - 1}) -> lys_mat_num

cancer <- rep(names(lys_mat_num), times = lys_mat_num)
cluster_col <- c("#00008B", "#00FF00", "#E31A1C", "#8B1C62", "#00F5FF", "#912CEE") #RColorBrewer::brewer.pal(name = "Set1", n = 6)
names(cluster_col) <- c(1,2,3,4,5,6)

# complexheatmap ----------------------
ha = HeatmapAnnotation(
  df = data.frame(cluster = samples_subgroup[colnames(lys_mat)], cancer = cancer), 
  gap = unit(c(4,2), "mm"),
  col = list(
    cancer = pcc, 
    cluster = cluster_col)
)
# dev.off()
# draw(ha, 1: 9744)

Heatmap(
  lys_mat_scale, 
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"), space = "RGB"),
  name = "expression",
  
  # annotation
  top_annotation = ha, 
  
  # show row and columns
  show_row_names = T, 
  show_column_names = FALSE, 
  show_row_dend = F,
  show_column_dend = F, 
  
  row_names_gp = gpar(fontsize = 6),
  # clustering
  clustering_distance_columns = "pearson",
  clustering_method_columns = "ward.D",
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D"
) -> ht

pdf(file.path(expr_path_a, "cluster_lys_mrna_01.pdf"), width=10, height = 10)
draw(ht, show_heatmap_legend = F, show_annotation_legend = F)
decorate_annotation(
  "cancer", 
  {grid.text("cancer", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 12))}
)
decorate_annotation(
  "cluster", 
  {grid.text("cluster", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 12))}
)
dev.off()

gene_list_expr_tumor_lys %>% 
  tibble::enframe() %>% 
  dplyr::mutate(cluster_count = purrr::map(.x = value, .f = fn_subtype, samples_subgroup = samples_subgroup)) %>% 
  dplyr::select(-value) %>% 
  tidyr::unnest() -> cancer_cluster_freq

cancer_cluster_freq %>% 
  dplyr::group_by(name) %>% 
  dplyr::mutate(freq = count / sum(count)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = as.factor(cluster), y = name, fill = freq)) +
  geom_tile() +
  geom_text(aes(label = count)) +
  scale_fill_gradient(
    name = "Percent (%)",
    limit = c(0, 1),
    high = "red",
    low = "white",
    na.value = "white"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      color = 'black'
    ),
    axis.ticks = element_blank(),
    
    panel.background = element_blank(),
    panel.border = element_rect(
      color = "black", 
      size = 0.2, 
      fill = NA
    )
  ) +
  guides(
    fill = guide_legend(
      title = "Percent (%)",
      title.position = 'left',
      title.theme = element_text(angle = 90, vjust = 2, size = 10),
      reverse = T,
      keywidth = 0.6,
      keyheight = 0.7
    )
  ) -> tile_plot

ggsave(filename = "cluster_lys_subtype.pdf", device = "pdf", plot = tile_plot, path = expr_path_a, width = 5, height = 6)

save.image(file = file.path(expr_path_a, ".rda_03_g_cluster_of_cluster"))
load(file = file.path(expr_path_a, ".rda_03_a_gene_expr.rda"))

