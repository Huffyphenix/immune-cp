library(magrittr)
library(ggplot2)
tcga_dir <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
pcbc_dir <- "/home/cliu18/liucj/projects/6.autophagy/synapse/PCBC"
csc_dir <- "/home/cliu18/liucj/projects/6.autophagy/08_cancer_stem_cell"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")

scores_gsva_sci_cor <- readr::read_rds(path = file.path(csc_dir, ".rds_04_stemness_index_scores_gsva_sci_cor.rds.gz"))

# basic function -----------------------------------
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
      )) 
  # %>%
  #   dplyr::group_by(sample) %>%
  #   dplyr::filter(n() >= 2, length(unique(type)) == 2) %>%
  #   dplyr::ungroup()
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

# sci tumor vs. normal----------------------------------
sci <- scores_gsva_sci_cor %>% dplyr::select(cancer_types, sci)
fn_sci_tn <- function(.x, .y){
  print(.x)
  # .x <- .te$cancer_types
  # .y <- .te$sci[[1]]
  
  .y %>% fn_tn_pair() -> .d
  
  sample_type_summary <- table(.d$type) %>% as.numeric()
  if (gtools::invalid(sample_type_summary) || any(sample_type_summary < c(10, 10))) return(NULL)
  
  .d %>% fn_t_test()
}

sci %>%
  dplyr::mutate(tn_diff = purrr::map2(.x = cancer_types, .y = sci, .f = fn_sci_tn)) %>%
  dplyr::filter(!purrr::map_lgl(.x = tn_diff, is.null)) %>% 
  tidyr::unnest(tn_diff) -> sci_tn_diff

CPCOLS <- c("#191970", "#EE0000")

sci_tn_diff %>% 
  dplyr::filter(pval < 0.05) %>%
  tidyr::unnest(sci) %>% 
  dplyr::mutate(or = glue::glue("{cancer_types} p-value: {signif(-log10(pval), digits = 3)}")) %>% 
  dplyr::mutate(or = reorder(or, pval)) %>% 
  fn_tn_pair() %>% 
  ggplot(aes(x = type, y = score)) +
  geom_boxplot(aes(color = type)) +
  scale_color_manual(name = "Type", values = CPCOLS) +
  facet_wrap(facets =  ~or) +
  labs(y = "Stemness score", title = "Stemness score between Tumor and Normal") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    
    strip.background = element_rect(color = "black", fill = "white"),
    
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) -> stemness_tn_p

ggsave(filename = "fig_01_stemness_tn_boxplot.pdf", plot = stemness_tn_p, device = "pdf", path = csc_dir, width = 7, height = 8)
  

# gsva tumor vs. normal---------------------------------------
gsva <- scores_gsva_sci_cor %>% dplyr::select(cancer_types, gsva)
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

gsva %>% 
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
  scale_y_discrete(limit = path_rank, labels = rev(c("Pathway", "ATG", "ATG_involved", "ATG_core", "Lys_deg", "Lys", "Lys_compo"))) +
  theme_bw() +
  labs(x = "", y = "Gene set", title = "ATG related gsva Norml vs. Tumor") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  ) -> gsva_tn_p

ggsave(filename = "fig_02_gsva_point.pdf", plot = gsva_tn_p, device = "pdf", path = csc_dir, width = 7, height = 4)

# gsea to validate gsva result --------------------------------------------
gsea <- readr::read_rds(path = file.path(expr_path_a, ".rds_02_i_gsea_cancer_es.rds.gz"))

gsea %>% dplyr::filter(q_val < 0.05) %>% dplyr::mutate(q_val = -log10(q_val)) -> plot_ready

plot_ready %>% 
  ggplot(aes(x = cancer_types, y = GS, size = q_val, color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_y_discrete(limits = c("pathway", "lys_deg", "lys_com", "lys"), labels = rev(c("Pathway", "Lys_deg", "Lys_compo", "Lys"))) +
  scale_size(name = "Q-value") + 
  ggthemes::theme_gdocs() +
  labs(x = "", y = "Gene set", title = "ATG related gsea Norml vs. Tumor") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5) 
  ) -> gsea_tn_p

ggsave(filename = "fig_03_gsea_point.pdf", plot = gsea_tn_p, device = "pdf", path = csc_dir, width = 10, height = 7)

# gsea and gsva result are not consistent.

# gsva score and stemness scores---------------------------

pathway_diff <- scores_gsva_sci_cor %>% dplyr::select(cancer_types, pathway_diff)

pathway_diff %>% 
  tidyr::unnest() %>% 
  dplyr::mutate(cor_pval = -log10(cor_pval), t_pval = -log10(t_pval)) %>% 
  dplyr::mutate(cor_pval = ifelse(cor_pval > 30, 30, cor_pval)) %>% 
  dplyr::mutate(t_pval = ifelse(t_pval > 20, 20, t_pval)) %>% 
  dplyr::mutate(pathway = plyr::revalue(pathway, replace = path_map)) %>% 
  dplyr::mutate(repel = glue::glue("{cancer_types}:{pathway}")) -> plot_ready

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = mean(cor_coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(pathway) %>% 
  dplyr::summarise(m = mean(cor_coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(pathway) -> path_rank

c("pathway" = "Pathway", "atg_core" = "ATG_core", "atg" = "ATG", "atg_form" = "ATG_involved", "lys_com" = "Lys_com", "lys" = "Lys", "lys_deg" = "Lys_deg") -> path_map

plot_ready %>%
  dplyr::filter(cor_pval > -log10(cor_pval)) %>% 
  ggplot(aes(x = cancer_types, y = pathway)) +
  geom_tile(aes(fill = cor_coef)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Coef") +
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = path_rank, labels = plyr::revalue(path_rank, path_map)) +
  labs(x = "", y = "Gene set") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) -> gsva_stemness_coef

CPCOLS <- c("#000080", "#F7F7F7", "#CD0000", "#F8F8FF")

plot_ready %>%
  dplyr::mutate(star = dplyr::case_when(
    dplyr::between(cor_pval, -log10(0.05), 10) ~ "*",
    dplyr::between(cor_pval, 10, 20) ~ "**",
    cor_pval > 20 ~ "***",
    TRUE ~ "."
  )) %>% 
  dplyr::filter(cor_pval > -log10(0.05)) %>% 
  ggplot(aes(x = cancer_types, y = pathway)) +
  geom_tile(aes(fill = cor_coef)) +
  geom_text(aes(label = star), color = CPCOLS[4]) +
  scale_fill_gradient2(low = CPCOLS[1], mid = CPCOLS[2], high = CPCOLS[3], name = "Coef") +
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = path_rank, labels = plyr::revalue(path_rank, path_map)) +
  labs(x = "", y = "Gene set", title = "GSVA vs. Stemness score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    
    # legend.position = "bottom",
    # legend.direction = "horizontal"
  ) -> gsva_stemness_coef

ggsave(filename = "fig_04_gsea_stemness_coef.pdf", plot = gsva_stemness_coef, device = "pdf", path = csc_dir, height = 4, width = 8)

plot_ready %>% 
  dplyr::mutate(sig = ifelse(cor_pval > -log10(0.05) & abs(cor_coef) > 0.3, "Sig", "Not Sig")) %>% 
  ggplot(aes(x = cor_coef, y = cor_pval)) +
  geom_point(aes(color = sig)) +
  scale_color_manual(values = c("grey", "red")) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(plot_ready, abs(cor_coef) > 0.3, cor_pval > -log10(0.05), pathway != "Pathway"), 
    aes(label = repel)) +
  theme_bw() 
  
plot_ready %>%
  dplyr::filter(t_pval > -log10(t_pval)) %>% 
  ggplot(aes(x = cancer_types, y = pathway)) +
  geom_point(aes(color = t_diff, size = t_pval)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "Diff") +
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = path_rank, labels = plyr::revalue(path_rank, path_map)) +
  scale_size(name = "P-value") +
  labs(x = "", y = "Gene set") +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) -> gsva_stemness_diff
ggsave(filename = "fig_05_gsea_stemness_diff.pdf", plot = gsva_stemness_diff, device = "pdf", path = csc_dir)

# stemness score for the autophagy gene mrna expression----------------

# diff
fn_stemness_high_low <- function(.x, .y){
  # .x <- .te$cancer_types
  # .y <- .te$sci[[1]]
  
  .y %>% 
    fn_tn_pair() %>% 
    dplyr::group_by(type) %>% 
    dplyr::mutate(hln = ifelse(score > 0, "HS", "LS")) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(hln = ifelse(type == "Normal", "Normal", hln))
}
fn_merge_expr_sci <- function(.x, .y){
  # .x <- .te$filter_expr[[1]]
  # .y <- .te$sci[[1]]
  
  .x %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(barcode, expr, -symbol) %>% 
    dplyr::filter(!is.na(expr)) %>%
    dplyr::inner_join(.y, by = "barcode") %>% 
    dplyr::mutate(expr = log2(expr + 0.01)) -> .d
  
  .d %>%
    dplyr::filter(symbol == "SNCA") %>%
    dplyr::select(expr, score) %>%
    cor(method = "spearman")
  
  .d %>% 
    dplyr::filter(hln %in% c("LS", "HS")) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          expr = t.test(expr ~ hln, data = .),
          error = function(e){1},
          warning = function(w){1}
        )
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::select(symbol, HS_expr = estimate1, LS_expr = estimate2, p.value, fdr) %>% 
    dplyr::mutate(fold = HS_expr - LS_expr) %>% 
    dplyr::mutate(sig = ifelse(fdr < 0.05 & abs(fold) > log2(1.5), "Sig", "Non-sig")) -> .dd
  .dd
}
  
sci %>% 
  dplyr::mutate(sci = purrr::map2(.x = cancer_types, .y = sci, .f = fn_stemness_high_low)) -> sci_diff

gene_list_expr <- readr::read_rds(path = file.path(expr_path_a, ".rds_03_a_gene_list_expr.rds.gz"))

gene_list_expr %>% 
  dplyr::inner_join(sci_diff, by = "cancer_types") %>% 
  # dplyr::filter(cancer_types == "TGCT") %>% 
  dplyr::mutate(hl_diff = purrr::map2(.x = filter_expr, .y = sci, .f = fn_merge_expr_sci)) ->
  gene_list_sci_diff

gene_list_sci_diff %>% 
  dplyr::select(1,4) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(sig == "Sig") %>% 
  dplyr::mutate(fdr = -log10(fdr)) %>% 
  dplyr::mutate(fdr = ifelse(fdr > 30, 30, fdr)) %>% 
  dplyr::mutate(fold = dplyr::case_when(fold > 2 ~ 2, fold < -2 ~ -2, TRUE ~ fold)) -> plot_ready

plot_ready %>% 
  ggplot(aes(x = cancer_types, y = symbol, size = fdr, color = fold)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")

# coef
fn_gene_coef <- function(.x, .y){
  # .x <- .te$filter_expr[[1]]
  # .y <- .te$sci[[1]]
  
  .x %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(barcode, expr, -symbol) %>% 
    dplyr::filter(!is.na(expr)) %>%
    dplyr::inner_join(.y, by = "barcode") %>% 
    dplyr::mutate(expr = log2(expr + 0.01)) -> .d
  
  .d %>%
    dplyr::filter(symbol == "BECN1") %>%
    dplyr::select(expr, score) %>%
    dplyr::do(cor.test(.$expr, .$score, method = "spearman") %>% broom::tidy())
  
  .d %>% 
    dplyr::filter(hln %in% c("LS", "HS")) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          expr = cor.test(.$expr, .$score, method = "spearman"),
          error = function(e){1},
          warning = function(w){1}
        )
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::select(symbol, coef = estimate, p.value, fdr) %>% 
    dplyr::mutate(sig = ifelse(fdr < 0.05 & abs(coef) > 0.5, "Sig", "Non-sig")) -> .dd
  
  # .dd %>% 
  #   ggplot(aes(y = -log10(fdr), x = coef)) +
  #   geom_point(aes(color = sig)) +
  #   ggrepel::geom_text_repel(data = dplyr::filter(.dd, sig == "Sig"), aes(label = symbol))
}

gene_list_expr %>% 
  dplyr::inner_join(sci_diff, by = "cancer_types") %>% 
  # dplyr::filter(cancer_types == "TGCT") %>% 
  dplyr::mutate(hl_diff = purrr::map2(.x = filter_expr, .y = sci, .f = fn_gene_coef)) ->
  gene_list_sci_coef

gene_list_sci_coef %>% 
  dplyr::select(1,4) %>% 
  tidyr::unnest() %>% 
  dplyr::mutate(sig = ifelse(fdr < 0.05 & abs(coef) > 0.5, "Sig", "Non-sig")) %>% 
  dplyr::filter(sig == "Sig") %>% 
  dplyr::mutate(fdr = -log10(fdr)) %>% 
  dplyr::mutate(fdr = ifelse(fdr > 30, 30, fdr)) -> plot_ready

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = sum(coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(m = sum(coef)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>%
  dplyr::arrange(color, m) %>% 
  dplyr::select(symbol, color, status, type) -> gene_rank

plot_ready %>% 
  ggplot(aes(x = coef, y = fdr, color = cancer_types)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = symbol))

CPCOLS <- c("#000080", "#F0F0F0", "#FF0000")

plot_ready %>% 
  ggplot(aes(y = cancer_types, x = symbol, color = coef, size = fdr)) +
  geom_point() +
  scale_color_gradient2(low = CPCOLS[1] , mid = CPCOLS[2], high = CPCOLS[3], name = "Coef") +
  scale_size(name = "FDR") +
  scale_x_discrete(limit = gene_rank$symbol) +
  scale_y_discrete(limit = cancer_rank) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = gene_rank$color)
  ) +
  labs(x = "", y = "", title = "Gene expression vs. Stemness score") -> expr_stemness_p
ggsave(filename = "fig_06_gene_expr_stemness.pdf", plot = expr_stemness_p, device = "pdf", path = csc_dir, width = 9, height = 4)

# gsva score correlated with ten rppa pathway--------------------
rppa <- readr::read_rds(file.path(tcga_path, "pancan32_rppa_score.rds.gz"))
fn_rppa_gsva_merge <- function(.x, .y){
  # .x <- .te$gsva[[1]]
  # .y <- .te$rppa[[1]]
  .y <- .y %>% dplyr::rename(sample = barcode)
  
  .x %>% 
    tidyr::gather(key = barcode, value = gsva, -set) %>% 
    fn_tn_pair() %>% 
    dplyr::select(set, gsva, sample) %>% 
    dplyr::distinct(set, sample, .keep_all = T) %>% 
    dplyr::inner_join(.y, by = "sample") -> .d
  
  .d %>% 
    dplyr::group_by(set, pathway) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          cor.test(.$gsva, .$score, method = "spearman"),
          error = function(e){1},
          warning = function(e){1}
        )
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(set, pathway, coef = estimate, pval = p.value)
  
}

gsva %>% 
  dplyr::inner_join(rppa, by = "cancer_types")  %>% 
  dplyr::mutate(out = purrr::map2(gsva, rppa, .f = fn_rppa_gsva_merge)) ->
  gsva_rppa_coef

gsva_rppa_coef %>% 
  tidyr::unnest(out) %>% 
  dplyr::filter(pval < 0.05, abs(coef) > 0.3) %>% 
  dplyr::mutate(pval = -log10(pval)) %>% 
  dplyr::mutate(pval = ifelse(is.infinite(pval), 20, pval)) %>% 
  dplyr::group_by(set, pathway) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() -> plot_ready

plot_ready %>% 
  ggplot(aes(x = set, y = pathway, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n)) +
  scale_fill_gradient2() +
  theme_bw() +
  labs(x = "", y = "", title = "GSVA correlaton with pathway") -> gsva_rppa_p

ggsave(filename = "fig_07_gsva_rppa.pdf", plot = gsva_rppa_p, device = "pdf", path = csc_dir, height = 4, width = 6)

# gsva stage  ------------------------------------------
fn_gsva_stage <- function(.x, .y){
  # .x <- .te$gsva[[1]]
  # .y <- .te$stage[[1]]
  .y <- .y %>% dplyr::rename(sample = barcode)
  .x %>% 
    tidyr::gather(key = barcode, value = gsva, -set) %>% 
    fn_tn_pair() %>% 
    dplyr::select(set, gsva, sample) %>% 
    dplyr::distinct(set, sample, .keep_all = T) %>% 
    dplyr::inner_join(.y, by = "sample") -> .d
  
  .d %>% 
    dplyr::group_by(set) %>% 
    dplyr::do(
      broom::tidy(
        oneway.test(gsva ~ stage, data = .)
      )
    ) %>% 
    dplyr::ungroup()  -> .dp
  
  comp_list <- list(c("Stage I", "Stage II"), c("Stage II", "Stage III"), c("Stage III", "Stage IV"))
  
  .d %>% 
    dplyr::mutate(stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))) %>% 
    dplyr::filter(set == "lys_deg") %>% 
    ggpubr::ggboxplot(x = "stage", y = "gsva", color = "stage", pallete = "jco") +
    ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test") + 
    ggpubr::stat_compare_means(method = "anova", label.y = 1) +
    labs(title = "ACC Lys_deg gsva score stage change") -> .p
  
  ggsave(filename = "fig_10_acc_lys_deg_core.pdf", device = "pdf", plot = .p, path = csc_dir, width = 6, height = 5)
  
  .d %>% 
    dplyr::group_by(set, stage) %>% 
    dplyr::summarise(m = mean(gsva)) %>% 
    dplyr::arrange(set, stage) %>% 
    dplyr::group_by(set) %>% 
    dplyr::mutate(a = m - m[1]) %>% 
    dplyr::mutate(b = m - m[2]) %>% 
    dplyr::mutate(c = m - m[3]) %>% 
    dplyr::select(set, a, b, c) %>% 
    dplyr::mutate(d = a[2], e = b[3], f = c[4]) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(set, d, e, f) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(g = dplyr::case_when(d > 0 && e > 0 && f > 0 ~ "up",
                                       d < 0 && e < 0 && f < 0 ~ "down",
                                       TRUE ~ "MISC")) %>% 
    dplyr::select(set, direction = g) %>% 
    dplyr::left_join(.dp, by = "set")
    
}

gsva %>% 
  dplyr::inner_join(clinical_stage, by = "cancer_types") %>% 
  dplyr::filter(cancer_types == "ACC") -> .te
  dplyr::mutate(out = purrr::map2(.x = gsva, .y = stage, .f = fn_gsva_stage)) -> gsva_stage

gsva_stage %>% 
  tidyr::unnest(out) %>% 
  dplyr::filter(p.value < 0.05) -> plot_ready

CPCOLS <- c("#CD0000", "#00008B", "#838B8B")

plot_ready %>% 
  dplyr::mutate(set = plyr::revalue(set, path_map)) %>% 
  ggplot(aes(x = cancer_types, y = set, color = direction, size = -log10(p.value))) +
  geom_point() +
  scale_size(name = "Anova P-value") +
  scale_color_manual(name = "Direction", values = CPCOLS, limits = c("up", "down", "MISC"), labels = c("Up", "Down", "MISC")) +
  theme_bw() +
  labs(x = "", y = "Gene set") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) -> gsva_stage_p

ggsave(filename = "fig_08_gsva_stage.pdf", plot = gsva_stage_p, device = "pdf", path = csc_dir, width = 6, height = 4)

# gsva score survival-------------------------------
clinical <- readr::read_rds(path = file.path(tcga_path,"pancan34_clinical.rds.gz"))
fn_fit <- function(.z, .d){
  # .z <- "atg"
  .d %>% dplyr::filter(set == .z) -> .df
  survival::survdiff(survival::Surv(time, status) ~ group, data = .df) -> .d_diff
  
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
  
  .fit_x <- survival::survfit(survival::Surv(time, status) ~ group, data = .df , na.action = na.exclude)
  
  survminer::ggsurvplot(.fit_x, data = .df, pval = T, pval.method = T) -> .p
  
  tibble::tibble(kmp = kmp, p = list(.p$plot))
}
fn_survival <- function(.x, .y){
  # .x <- .te$gsva[[1]]
  # .y <- .te$clinical[[1]]
  
  .y <- .y %>% dplyr::rename(sample = barcode)
  
  .x %>% 
    tidyr::gather(key = barcode, value = gsva, -set) %>% 
    fn_tn_pair() %>% 
    dplyr::filter(type == "Tumor") %>% 
    dplyr::select(set, gsva, sample) %>% 
    dplyr::distinct(set, sample, .keep_all = T) %>% 
    dplyr::inner_join(.y, by = "sample") %>% 
    dplyr::select(set, gsva, sample, status = os_status, time = os_days) %>% 
    dplyr::mutate(status = plyr::revalue(status, replace = c("Alive" = 0, "Dead" = 1))) %>% 
    dplyr::mutate(status = as.numeric(status)) %>% 
    tidyr::drop_na() %>% 
    dplyr::group_by(set) %>% 
    dplyr::mutate(group = as.factor(ifelse(gsva <= median(gsva),"Low", "High"))) %>% 
    dplyr::ungroup() -> .d

  
  .d %>% 
    dplyr::group_by(set) %>% 
    dplyr::do(
      broom::tidy(
        tryCatch(
          survival::coxph(survival::Surv(time, status) ~ gsva, data = ., na.action = na.exclude),
          error = function(e){1},
          warning = function(e){1})
      )
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(hazard_ratio = exp(estimate)) -> .dp
  
  .dp %>% 
    dplyr::mutate(out = purrr::map(.x = set, .f = fn_fit, .d)) %>% 
    tidyr::unnest() %>% 
    dplyr::select(set, estimate, coxp = p.value, kmp, p) 
}

gsva %>% 
  dplyr::inner_join(clinical, by = "cancer_types") %>% 
  # dplyr::filter(cancer_types == "LGG") %>% 
  dplyr::mutate(out = purrr::map2(.x = gsva, .y = clinical, .f = fn_survival)) %>% 
  tidyr::unnest(out) -> gsva_survival

gsva_survival %>% 
  dplyr::filter(kmp < 0.05) %>% 
  # dplyr::slice(10) %>% .$p
  dplyr::mutate(direction = ifelse(estimate > 0, "High_Worse", "Low_Worse")) %>% 
  ggplot(aes(x = cancer_types, y = set, size = -log10(coxp), color = direction)) +
  geom_point() +
  theme(
    axis.text.x = element_text(angle = 90)
  )


# save----------------------------------------------------------------------
save.image(file = file.path(csc_dir, ".rda_05_stemness_analysis.rda"))
load(file = file.path(csc_dir, ".rda_05_stemness_analysis.rda"))








