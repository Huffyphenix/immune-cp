
# library -----------------------------------------------------------------
library(magrittr)
library(GSVA)
library(ggplot2)

# path --------------------------------------------------------------------
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

# load data ---------------------------------------------------------------
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

gene_sets <- list(
  LC3 = dplyr::filter(gene_list, process == "LC3"),
  ATG5_ATG12_ATG16 = dplyr::filter(gene_list, process == "ATG5-ATG12-ATG16"),
  ULK_Complex = dplyr::filter(gene_list, process == "ULK-complex"),
  ATG9_cycling_system = dplyr::filter(gene_list, process == "ATG9-cycling-system"),
  PI3K_III_complex = dplyr::filter(gene_list, process == "PI3K-III-complex"),
  Negative = dplyr::filter(gene_list, process == "Negative"),
  Positive = dplyr::filter(gene_list, process == "Positive"),
  TRAPP = dplyr::filter(gene_list, process == "TRAPP"),
  Membrane_delivery = dplyr::filter(gene_list, process == "Membrane-delivery"),
  Mitophagy = dplyr::filter(gene_list, process == "Mitophagy"),
  Lys_comp = dplyr::filter(gene_list, process == "Lys_comp"),
  Lys_deg = dplyr::filter(gene_list, process == "Lys_deg"),
  atg = dplyr::filter(gene_list, type == "atg"),
  lys = dplyr::filter(gene_list, type == "lys"),
  atg_core = dplyr::filter(gene_list, pathway == "autophagesome formation-core")
) %>% purrr::map("symbol")

# gsva score --------------------------------------------------------------
fn_gsva <- function(.x, .y, gene_sets = gene_sets){
  # .x <- .te$cancer_types
  # .y <- .te$filter_expr[[1]]
  print(.x)
  .y %>% 
    tidyr::drop_na() %>% 
    dplyr::select( -entrez_id) -> .d
  
  .d_mat <- as.matrix(.d[,-1])
  rownames(.d_mat) <- .d$symbol
  
  .es_dif <- gsva(.d_mat, gene_sets, method = "gsva", mx.diff = TRUE, verbose = FALSE, parallel.sz = 1)
  
  .es_dif$es.obs %>% 
    as.data.frame() %>% 
    tibble::as_tibble() %>% 
    tibble::add_column(set = rownames(.es_dif$es.obs), .before = 1) -> .d_es
}

cluster <- multidplyr::create_cluster(nrow(expr))
expr %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("GSVA") %>%
  multidplyr::cluster_assign_value("fn_gsva", fn_gsva)  %>%
  multidplyr::cluster_assign_value("gene_sets", gene_sets)  %>%
  dplyr::mutate(gsva = purrr::map2(.x = cancer_types, .y = expr, .f = fn_gsva, gene_sets = gene_sets)) %>% 
  dplyr::select(-expr) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> gene_list_gsva
parallel::stopCluster(cluster)

readr::write_rds(x = gene_list_gsva, path = file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"), compress = "gz")

gene_list_gsva <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_i_gsva_gene_list_gsva.rds.gz"))

# tumor vs normal ---------------------------------------------------------
fn_tn <- function(.s){
  .s %>% stringr::str_sub(start = 14, end = 14)
}
fn_tn_pair <- function(.d){
  .d %>% 
    dplyr::mutate(
      sample = stringr::str_sub(
        string = barcode,
        start = 1,
        end = 12
      )
    ) %>%
    dplyr::mutate(type = purrr::map_chr(barcode, .f = fn_tn)) %>% 
    dplyr::filter(type %in% c("0", "1")) %>%
    dplyr::mutate(type = plyr::revalue(
      x = type,
      replace = c("0" = "Tumor", "1" = "Normal"),
      warn_missing = F
    )) %>%
    dplyr::group_by(sample) %>%
    dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
    dplyr::ungroup()
}
fn_t_test <- function(.d){
  .d %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          expr = t.test(score ~ type, data = .),
          error = function(e){1},
          warning = function(e){1}
        ))) %>% 
    dplyr::mutate(diff = ifelse(tibble::has_name(., "estimate"), estimate, 0)) %>% 
    dplyr::mutate(pval = ifelse(tibble::has_name(., "p.value"), p.value, 1)) -> .dd
  
  if (nrow(.dd) > 1) dplyr::select(.dd, 1, diff, pval) else dplyr::select(.dd, diff, pval)
  
}
fn_gsva_tn <- function(.x, .y){
  # .x <- .te$cancer_types
  # .y <- .te$gsva[[1]]
  print(.x)
  
  .y %>% 
    tidyr::gather(key = "barcode", value = "score", -set) %>% 
    fn_tn_pair() -> .d
  
  sample_type_summary <- table(.d$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) || any(sample_type_summary < c(10, 10))) return(NULL)
  
  .d %>% dplyr::group_by(set) %>% fn_t_test() %>% dplyr::ungroup()
}

gene_list_gsva %>% 
  dplyr::mutate(tn_diff = purrr::map2(.x = cancer_types, .y = gsva, .f = fn_gsva_tn)) %>% 
  dplyr::filter(!purrr::map_lgl(.x = tn_diff, is.null)) %>% 
  tidyr::unnest(tn_diff) -> gsva_tn_diff

gsva_tn_diff %>% 
  dplyr::mutate(pval = -log10(pval), diff = -diff) %>% 
  dplyr::filter(pval > -log10(0.05)) -> plot_ready

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m_diff = mean(diff)) %>% 
  dplyr::arrange(dplyr::desc(m_diff)) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(set) %>% 
  dplyr::summarise(m_diff = mean(diff)) %>% 
  dplyr::arrange(dplyr::desc(m_diff)) %>% 
  dplyr::pull(set) -> path_rank

CPCOLS <- c("#000080", "#F7F7F7", "#CD0000", "#F8F8FF")
plot_ready %>% 
  dplyr::mutate(star = dplyr::case_when(
    dplyr::between(pval, -log10(0.05), 5) ~ "*",
    dplyr::between(pval, 5, 10) ~ "**",
    pval > 10 ~ "***",
    TRUE ~ "."
  )) %>% 
  ggplot(aes(x = cancer_types, y = set)) +
  geom_tile(aes(fill = diff)) +
  geom_text(aes(label = star), color = CPCOLS[4]) +
  scale_fill_gradient2(low = CPCOLS[1], mid = CPCOLS[2], high = CPCOLS[3], name = "Diff") +
  scale_size(name = "P-value") + 
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = path_rank) +
  theme_bw() +
  labs(x = "", y = "Gene set", title = "ATG related gsva Norml vs. Tumor") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  ) -> p
ggsave(filename = file.path(expr_path_a, "fig_summary_gsva_tumor_vs_normal.pdf"), plot = p)

# tumor vs normal boxplot -------------------------------------------------
fn_boxplot_data <- function(.x, .set){
  # .x <- .te$gsva[[1]]
  # .set <- "atg"
  
  .x %>% 
    dplyr::filter(set == .set) %>% 
    tidyr::gather(barcode, gsva, -set) %>% 
    fn_tn_pair() -> .d
  
  sample_type_summary <- table(.d$type) %>% as.numeric()
  if (length(sample_type_summary) < 2 || gtools::invalid(sample_type_summary) || any(sample_type_summary < c(10, 10))) return(NULL)
  .d %>% dplyr::select(-barcode)
}
gene_list_gsva %>% 
  dplyr::mutate(gsva = purrr::map(.x = gsva, .f = fn_boxplot_data, .set = "atg")) %>% 
  dplyr::filter(!purrr::map_lgl(.x = gsva, .f = is.null)) %>% 
  tidyr::unnest() -> atg_dist

atg_dist %>% 
  dplyr::distinct(sample, type, .keep_all = T) %>%
  # dplyr::select(-sample) %>% 
  tidyr::spread(key = type, value = gsva) %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = mean(Tumor) - mean(Normal)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

atg_dist %>% 
  ggplot(aes(x = cancer_types, y = gsva, fill = type)) +
  geom_boxplot(outlier.colour = NA) +
  scale_x_discrete(limits = cancer_rank) +
  scale_fill_manual(name = "Type", values = CPCOLS[c(1,3)]) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  labs(x = "", y = "ATG score") -> atg_box


gene_list_gsva %>% 
  dplyr::mutate(gsva = purrr::map(.x = gsva, .f = fn_boxplot_data, .set = "lys")) %>% 
  dplyr::filter(!purrr::map_lgl(.x = gsva, .f = is.null)) %>% 
  tidyr::unnest() -> lys_dist


lys_dist %>% 
  ggplot(aes(x = cancer_types, y = gsva, fill = type)) +
  geom_boxplot(outlier.colour = NA) +
  scale_x_discrete(limits = cancer_rank) +
  scale_fill_manual(name = "Type", values = CPCOLS[c(1,3)]) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "", y = "Lys score") -> lys_box

ggplot_gtable(data = ggplot_build(atg_box)) -> tmp
which(sapply(tmp$grobs, function(x) x$name) == "guide-box") -> leg
tmp$grobs[[leg]] -> mylegend
library(gridExtra)
grid.arrange(arrangeGrob(atg_box + theme(legend.position = "none"),
                         lys_box + theme(legend.position = "none"),
                         ncol = 1),
             mylegend, nrow = 2,heights = c(15, 1)) -> p
ggsave(filename = "fig_summary_gsva_atg_lys_tumor_vs_normal.pdf", path = expr_path_a, plot = p, device = "pdf", width = 5, height = 4)

# gene set corelation -----------------------------------------------------
gene_list_gsva %>% 
  dplyr::mutate(gsva = purrr::map(.x = gsva, .f = function(.x){
                                    .x %>% tidyr::gather(barcode, gsva, -set) %>% 
                                      tidyr::spread(set, gsva)
                                    })) %>% 
  tidyr::unnest()  %>% 
  dplyr::select(-1, -2) %>% 
  cor() -> gene_list_gsva_corr

set_rank <- c("ULK_Complex", "PI3K_III_complex", "ATG5_ATG12_ATG16", "ATG9_cycling_system", "LC3", "Lys_comp", "Lys_deg", "atg", "lys")
gene_list_gsva_corr[set_rank, set_rank]  -> gcorr_target

xy_rank <- c("ULK_Complex" = "ULK Complex", "PI3K_III_complex" = "PI3KIII Complex", "ATG5_ATG12_ATG16" = "ATG5-ATG12", "ATG9_cycling_system" = "ATG9 Cycling", "LC3" = "LC3", "1" = "", "Lys_comp" = "LYS membrane", "Lys_deg" = "LYS Enzyme", "2" = "", "atg" = "ATG score", "lys" = "LYS score")

CPCOLS <- c("#191970", "#F8F8FF", "#FF4040")

gcorr_target %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "x") %>% 
  tibble::as_tibble() %>% 
  tidyr::gather(y, corr, -x) %>% 
  dplyr::mutate(corr = dplyr::case_when(
    corr > 0.5 ~ 0.5,
    TRUE ~ corr
  )) %>% 
  ggplot(aes(x = x, y = y, fill = corr)) +
  geom_tile() +
  scale_y_discrete(limits = rev(names(xy_rank)), labels = rev(xy_rank)) +
  scale_x_discrete(limits = names(xy_rank), labels = xy_rank) +
  scale_fill_gradient2(
    low = CPCOLS[1], 
    mid = CPCOLS[2], 
    high = CPCOLS[3], 
    limits = c(-0.1, 0.5), 
    breaks = seq(-0.1, 0.5, 0.1),
    label = c("-0.1", "0", "0.1","0.2", "0.3", "0.4",  "0.5")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
    
    axis.ticks = element_blank(),
    
    panel.grid = element_blank(),
    panel.background = element_blank()
    
  ) +
  guides(
    fill = guide_colorbar(
      title = "Pearson correlation",
      title.position = "right",
      title.theme = element_text(angle = -90),
      barwidth = 0.5, 
      barheight = 7,
      
      label.position = "left")
    ) +
  labs(x = "", y = "") +
  annotation_custom(grob = grid::textGrob("ATG core", rot = 90), xmin = -1.05, xmax = -1.05, ymin = 6.5, ymax = 11.5) +
  annotation_custom(grob = grid::linesGrob(), xmin = -0.9, xmax = -0.9, ymin = 6.5, ymax = 11.5) +
  annotation_custom(grob = grid::textGrob("Lysosome", rot = 90), xmin = -1.05, xmax = -1.05, ymin = 3.5, ymax = 5.5) +
  annotation_custom(grob = grid::linesGrob(), xmin = -0.9, xmax = -0.9, ymin = 3.5, ymax = 5.5) -> p

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
ggsave(filename = "fig_summary_gsva_atg_lys_corr.pdf", path = expr_path_a, plot = gt, device = "pdf")

# gsva survival -----------------------------------------------------------

clinical <- readr::read_rds(file.path(tcga_path, "pancan34_clinical.rds.gz"))
fn_clinical_merge <- function(.x, .y){
  .x <- .te$gsva[[1]]
  .y <- .te$clinical[[1]]
  .x %>% 
    tidyr::gather(key = barcode, value = gsva, -set) %>% 
    dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 12)) %>% 
    dplyr::mutate(type = stringr::str_sub(barcode, start = 14, end = 15)) %>% 
    dplyr::filter(type != "11") %>% 
    dplyr::distinct(set, sample, .keep_all = T) %>% 
    dplyr::select(-barcode, -type) %>% 
    dplyr::rename(barcode = sample) %>% 
    dplyr::inner_join(.y, by = "barcode") %>% 
    dplyr::select(set, gsva, barcode, time = os_days, status = os_status) %>% 
    tidyr::drop_na() %>% 
    dplyr::filter(time > 0) %>% 
    dplyr::mutate(status = plyr::revalue(status, c("Dead" = 1, "Alive" = 0)) %>% as.integer()) %>% 
    dplyr::group_by(set) %>% 
    dplyr::mutate(group = ifelse(gsva > median(gsva), "High", "Low")) %>% 
    dplyr::ungroup() -> .d
    
  .d %>% 
    dplyr::group_by(set) %>% 
    dplyr::filter(n() > 10) %>% 
    dplyr::do(
      broom::tidy(
        survival::coxph(survival::Surv(time, status) ~ gsva, data = .)
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(p.value) %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr"))
  
  .d %>% dplyr::filter(set == "lys") -> .dd
  fit_x <- survival::survfit(survival::Surv(time, status) ~ group, data = .dd , na.action = na.exclude)
  survminer::ggsurvplot(fit_x, data = .dd, pval = T, pval.method = T,
                        title = "ACC LYS suvival analysis",
                        xlab = "Survival in days",
                        ylab = 'Probability of survival')
  
}

gene_list_gsva %>% 
  dplyr::inner_join(clinical, by = "cancer_types") %>% dplyr::filter(cancer_types == "ACC") -> .te
  dplyr::mutate(survival = purrr::map2(.x = gsva, .y = clinical, .f = fn_clinical_merge)) %>% 
  tidyr::unnest(survival) %>% 
  dplyr::select(cancer_types, set, estimate, p.value, fdr) -> set_gsva_survival

set_gsva_survival %>% 
  dplyr::filter(p.value < 0.05) %>% 
  dplyr::mutate(worse = ifelse(estimate < 0, "H", "L")) %>% 
  dplyr::filter(set %in% names(xy_rank)) -> plot_ready

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::count() %>%
  dplyr::arrange(-n) %>% 
  dplyr::pull(cancer_types) -> cancer_rank
  
plot_ready %>% 
  ggplot(aes(x = cancer_types, y = set, size = -log10(p.value), color = worse)) +
  geom_point() +
  scale_y_discrete(limits = rev(names(xy_rank)), labels = xy_rank) +
  scale_x_discrete(limits = cancer_rank) +
  scale_size(name = "P-value") +
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
    # axis.text.y = element_text(color = gene_rank$color),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) -> p
ggsave(filename = file.path(expr_path_a, "fig_summary_gsva_survival.pdf"), plot = p, device = "pdf", width = 8, height  = 5)

