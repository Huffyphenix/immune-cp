library(methods)
library(magrittr)

tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")



load(file = file.path(expr_path_a, ".rda_03_h_coca.rda"))

# expr
expr_matrix %>% t() -> expr_matrix_t
factoextra::hcut(expr_matrix_t, k = 4, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> expr_hcut
cutree(expr_hcut, k = 4) -> expr_group

expr_group %>% 
  tibble::enframe(name = "sample", value = "group") %>% 
  dplyr::inner_join(sample_info, by = "sample") %>% 
  dplyr::distinct(sample, cancer_types, .keep_all = T) %>% 
  ggplot(aes(x = as.factor(group), y = cancer_types)) +
  geom_tile()

fn_encode <- function(.x){
  .d <- tibble::tibble()
  if(.x == 1) {.d <- tibble::tibble(a = 1,b = 0, c = 0, d= 0)}
  if(.x == 2) {.d <- tibble::tibble(a = 0,b = 1, c = 0, d = 0)}
  if(.x == 3) {.d <- tibble::tibble(a = 0,b = 0, c = 1, d = 0)}
  if(.x == 4) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 1)}
  .d
}

expr_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> expr_encode


cnv_matrix %>% t() -> cnv_matrix_t
factoextra::hcut(cnv_matrix_t, k = 4, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> cnv_hcut
cutree(cnv_hcut, k = 4) -> cnv_group

cnv_group %>% 
  tibble::enframe(name = "sample", value = "group") %>% 
  dplyr::inner_join(sample_info, by = "sample") %>% 
  dplyr::distinct(sample, cancer_types, .keep_all = T) %>% 
  ggplot(aes(x = as.factor(group), y = cancer_types)) +
  geom_tile()

cnv_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> cnv_encode


methy_matrix %>% t() -> methy_matrix_t
factoextra::hcut(methy_matrix_t, k = 4, hc_func = 'hclust', hc_method = 'ward.D', hc_metric = 'pearson', stand = T) -> methy_hcut
cutree(methy_hcut, k = 4) -> methy_group

methy_group %>% 
  tibble::enframe(name = "sample", value = "group") %>% 
  dplyr::inner_join(sample_info, by = "sample") %>% 
  dplyr::distinct(sample, cancer_types, .keep_all = T) %>% 
  ggplot(aes(x = as.factor(group), y = cancer_types)) +
  geom_tile()

methy_group %>% 
  tibble::enframe(name = "sample") %>% 
  dplyr::mutate(encode = purrr::map(.x = value, .f = fn_encode)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-value) -> methy_encode

list(expr_encode, cnv_encode, methy_encode) %>% 
  purrr::reduce(.f = function(x, y){x %>% dplyr::inner_join(y, by = "sample")}, .init = tibble::tibble(sample = common_names[-1])) %>% 
  tidyr::gather(key = type, value = v, -sample) %>% 
  tidyr::spread(key = sample, value = v) %>% 
  dplyr::select(-type) -> cc_d


library(ConsensusClusterPlus)
ConsensusClusterPlus(cc_d %>% as.matrix(), maxK=20, reps=500,pItem=0.8,pFeature=1, title="ConsensusClusterplus4", clusterAlg="hc",distance="pearson",seed=1262118388.71279, plot = "pdf") -> cc_res

cc_res %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_h_coca_cc_cc_res.rds.gz"), compress = "gz")

fn_best_clust <- function(k, d = d){
  cc <- d[[k]][[3]]
  l <- length(cc)
  cm <- d[[k]][[1]]
  tibble::tibble(
    sample = names(cc),
    name = paste("V", c(1:l), sep = ""),
    group = cc) %>% 
    dplyr::arrange(group) -> rank_sample
  
  cm %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(
      sample1 = paste("V", c(1:l), sep = ""), .before = 1) %>% 
    tidyr::gather(key = sample2, value = sim, -sample1) -> plot_ready
  
  plot_ready %>% 
    ggplot(aes(x= sample1, y = sample2, fill = sim)) +
    geom_tile() +
    scale_x_discrete(limits = rank_sample$name) +
    scale_y_discrete(limits = rank_sample$name) +
    scale_fill_gradient(high  = "#1DEDF2", low = "#010235") +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      
      legend.position = "none"
    ) -> p
  ggsave(filename = paste("cc4_c",k, ".tif", sep = ""), plot = p, device = "tiff", width = 8, height = 8, path = "/extraspace/liucj/projects/6.autophagy/02_autophagy_expr/03_a_gene_expr")
}
cc_res -> d

cluster <- multidplyr::create_cluster(19)
tibble::tibble(k = 2:20) %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fn_best_clust", fn_best_clust)  %>%
  multidplyr::cluster_assign_value("d", d)  %>%
  dplyr::mutate(a = purrr::walk(.x = k, .f = fn_best_clust, d = d)) 
parallel::stopCluster(cluster)

clinical <- readr::read_rds(path = file.path(tcga_path,"pancan34_clinical.rds.gz"))

clinical_simplified <- 
  clinical %>% 
  dplyr::mutate(succinct = purrr::map(.x = clinical, function(.x ){.x %>% dplyr::select(barcode, os_days, os_status)})) %>% 
  dplyr::select(-clinical) %>% 
  tidyr::unnest() %>% 
  dplyr::rename(sample = barcode)

clinical_stage <- 
  readr::read_rds(path = file.path(tcga_path,"pancan34_clinical_stage.rds.gz")) %>% 
  dplyr::filter(n >= 40) %>% 
  dplyr::select(-n) %>% 
  tidyr::unnest() %>% 
  dplyr::rename(sample =barcode)

sample_info <- readr::read_rds(path = file.path(tcga_path, "sample_info.rds.gz"))

d[[4]][[3]] %>% tibble::enframe(name = "sample", value = "group") -> .d3

.d3 %>%
  dplyr::inner_join(sample_info, by = "sample") %>% 
  dplyr::distinct(sample, cancer_types, .keep_all = T) -> .d3_sample

# --------------------- draw merged heatmap-----------------
.d3_sample %>% dplyr::arrange(group) -> .d3_sample_sort
expr_matrix[,.d3_sample_sort$sample] -> expr_matrix_sort
cnv_matrix[, .d3_sample_sort$sample]  %>% apply(1, scale) %>% t() -> cnv_matrix_sort
methy_matrix[, .d3_sample_sort$sample] %>% apply(1, scale) %>% t() -> methy_matrix_sort

cluster_col <- c("#FF0000", "#191970", "#98F5FF", "#8B008B")
names(cluster_col) <- c(1,2,3,4)

library(ComplexHeatmap)
library(circlize)

ha = HeatmapAnnotation(
  df = data.frame(cluster = .d3_sample_sort$group, cancer = .d3_sample_sort$cancer_types),
  gap = unit(c(4,2), "mm"),
  col = list(
    cancer = pcc,
    cluster = cluster_col
  )
)

fn_draw_heatmap <- function(.x, .y, ha = ha){
  Heatmap(
    .y,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"), space = "RGB"),
    name = "expression",
    
    # annotation
    top_annotation = ha,
    
    # show row and columns
    show_row_names = F, 
    show_column_names = FALSE, 
    show_row_dend = F,
    show_column_dend = F, 
    
  
    clustering_distance_rows = "pearson",
    clustering_method_rows = "ward.D",
    
    cluster_columns = F,
    
    row_title = .x,
    row_title_rot = 90
  ) -> expr_ht
  
  .filename <- stringr::str_c("coca_heatmap", .x, "tiff", sep = " ") %>% stringr::str_replace_all(pattern = " ", replacement = ".")
  .pdf_file <- file.path("/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr//03_a_gene_expr", .filename)
  
  pdf(.pdf_file, width=10, height = 5)
  
  draw(expr_ht, show_heatmap_legend = F, show_annotation_legend = F)
  decorate_annotation(
    "cancer", 
    {grid.text("cancer", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 12))}
  )
  decorate_annotation(
    "cluster", 
    {grid.text("cluster", unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left", gp = gpar(fontsize = 12))}
  )
  
  dev.off()

}

fn_draw_heatmap(.x = "mRNA Expression", .y = expr_matrix_sort, ha = ha)
fn_draw_heatmap(.x = "Copy Number Variation", .y = cnv_matrix_sort, ha = ha)
fn_draw_heatmap(.x = "Methylation", .y = methy_matrix_sort, ha = ha)



c(6, 8, 1, 1) %>% 
  tibble::enframe() %>% 
  ggplot(aes(x = 1, y = value, fill = as.factor(name))) +
  geom_bar(stat = 'identity', width = 0.02) +
  # geom_text(label = c("C4", "C3", "C2", "C1"), position = position_stack(vjust = 0.5, )) +
  annotate("text", x = 1.5, y = c(0.5, 1.5, 6, 13), label = c("C4 (492)", "C3 (543)", "C2 (4086)", "C1 (3019)")) +
  scale_y_discrete(
    position = "right"
  ) +
  scale_fill_manual(
    values = unname(cluster_col),
    guide = F) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    # axis.text.y = element_blank(),
    # axis.text.x = element_text(
    #   size = 16,
    #   face = "bold"
    # ),
    axis.ticks = element_blank(),
    
    panel.background = element_blank(),
    panel.border = element_blank(),
    
    plot.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "mm")
    
  ) +
  labs(x = NULL, y = NULL) -> p
ggsave(filename = "coca_cc_legend.pdf", device = "pdf", plot = p, path = expr_path_a)
  

.d3 %>% 
  ggplot(aes(x = as.factor(group), y = 1, fill = as.factor(group))) +
  geom_tile() +
  scale_x_discrete(
    labels = c("C1", "C2", "C3", "C4"),
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
  coord_fixed(ratio = 0.2) -> p

ggsave(filename = "coca_anno.pdf", device = "pdf", plot = p, path = expr_path_a, width = 4, height = 1)



.d3_sample %>% 
  dplyr::group_by(cancer_types, group) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::mutate(freq = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = as.factor(group), y = cancer_types, fill = freq)) +
  geom_tile() +
  geom_text(aes(label = n)) +
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
  ) -> p

ggsave(filename = "coca_tile_subtype.pdf", device = "pdf", plot = p, path = expr_path_a, width = 5, height = 6)


.d3_sample %>% 
  dplyr::inner_join(clinical_stage, by = c("cancer_types", "sample")) %>% 
  dplyr::group_by(stage, group) %>% 
  ggplot(aes(x = as.factor(group), fill = stage)) +
  geom_bar(stat = 'count')

.d3_sample %>% 
  dplyr::inner_join(clinical_stage, by = c("cancer_types", "sample")) %>% 
  dplyr::group_by(stage, group) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(key = group, value = n) %>% 
  as.data.frame() -> .d3_sample_stage
rownames(.d3_sample_stage) <- .d3_sample_stage$stage
.d3_sample_stage[-1] -> .d3_sample_stage

# gplots::balloonplot(t(as.table(as.matrix(.d3_sample_stage))))

.d3_sample_stage %>% chisq.test() 
.d3_sample_stage %>% dplyr::select(1,2) %>% chisq.test()
.d3_sample_stage %>% dplyr::select(1,3) %>% chisq.test()
.d3_sample_stage %>% dplyr::select(2, 3) %>% chisq.test()

clinical_simplified %>% 
  dplyr::inner_join(.d3_sample, by = c("cancer_types", "sample")) %>% 
  dplyr::mutate(os_status = dplyr::recode(os_status, "Dead" = 1, "Alive" = 0)) %>% 
  dplyr::filter(! is.na(os_days), os_days > 0) -> .d

.d_diff <- survival::survdiff(survival::Surv(os_days, os_status) ~ group, data = .d)

kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)

fit_x <- survival::survfit(survival::Surv(os_days, os_status) ~ as.factor(group), data = .d , na.action = na.exclude)

survminer::ggsurvplot(fit_x, data = .d, pval=T, pval.method = T, 
                      palette = unname(cluster_col),
                      xlab = "Survival in days",
                      ylab = 'Probability of survival',
                      legend.title = "COCA",
                      legend.labs = c("C1", "C2", "C3", "C4")) -> p_survival
ggsave(filename = "coca_all_survival.pdf", plot = p_survival$plot, device = 'pdf', path = expr_path_a, width = 6, height = 6)

comp_list <- list(c(1,3),c(2, 4), c(2,3), c(3,4))

.d %>% 
  ggpubr::ggboxplot(x = "group", y = "os_days", color = "group", pallete = "jco") +
  ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test") + 
  ggpubr::stat_compare_means(method = "anova", label.y = 9000)

# nodes <- 
#   .d %>% 
#   dplyr::rename(cluster = group, group = cancer_types, name = sample) %>% 
#   dplyr::mutate(size = os_days / 30) %>% 
#   dplyr::mutate(size = dplyr::case_when(size < 1 ~ 1, size > 100 ~ 100, TRUE ~ round(size))) %>% 
#   dplyr::select(name, group, size, cluster) %>% 
#   as.data.frame() %>% 
#   head(100)
# 
# inter <- d[[4]][[1]] %>% as.data.frame() %>% tibble::as_tibble()
# names(inter) <- 0: {length(.d3_sample$sample) - 1}
# 
# inter %>% 
#   tibble::add_column(source = 0: {length(.d3_sample$sample) - 1}, .before = 1) %>% 
#   tidyr::gather(key = target, value = value, - source) %>% 
#   dplyr::filter(source != target) %>% 
#   dplyr::mutate(value = ifelse(value == 1, 1, 0)) %>% 
#   dplyr::mutate(target = as.numeric(target)) -> te
#  
# 
# 
# edges <- 
#   te %>% 
#   dplyr::filter(source %in% 0:99, target %in% 0:99) %>% 
#   as.data.frame()
# 
# networkD3::forceNetwork(Links = edges, Nodes = nodes , Source = "source", Target = "target", Value = "value", NodeID = "name", Nodesize = "size", Group = "group", opacity = 1) -> res_network
# 
# res_network
# data(MisLinks, MisNodes)
# forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
#              Target = "target", Value = "value", NodeID = "name", Nodesize = "size",
#              Group = "group", opacity = 1, zoom = T) 
#
#
#
save.image(file = file.path(expr_path_a, ".rda_03_h_coca_cc.rda"))
load(file.path(expr_path_a, ".rda_03_h_coca_cc4.rda"))






