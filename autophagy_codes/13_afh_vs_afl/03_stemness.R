# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# load path ---------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
stemness_path <- file.path(afhl_path, "03_stemness")

# stemness ----------------------------------------------------------------
csc <- readr::read_rds(file.path(tcga_path, "pancan33_csc_score.rds.gz"))
p62_sample_classification <- readr::read_rds(path = file.path(afhl_class, '.rds_01_p62_sample_classification.rds.gz')) %>% 
  dplyr::select(-barcode) %>% 
  tidyr::nest(-cancer_types)

csc %>% 
  dplyr::mutate(sci = purrr::map(
    .x = sci,
    .f = function(.x){
      .x %>% 
        dplyr::filter(stringr::str_sub(barcode, start = 14, 15) == "01") %>% 
        dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 12)) %>% 
        dplyr::distinct(sample, .keep_all = T) %>% 
        dplyr::select(-barcode)
    }
  )) %>% 
  dplyr::inner_join(p62_sample_classification, by = "cancer_types") %>% 
  dplyr::mutate(merge = purrr::map2(
    .x = sci,
    .y = data,
    .f = function(.x, .y) {
      .y %>% 
        dplyr::mutate(type = ifelse(rppa > median(rppa), "L", "H")) %>% 
        dplyr::inner_join(.x, by = "sample")
    }
  )) %>% 
  dplyr::select(cancer_types, merge) -> p62_csc



p62_csc %>% 
  dplyr::mutate(pval = purrr::map(
    .x = merge,
    .f = function(.x) {
      .x %>%
        dplyr::mutate(type = dplyr::case_when(
          rppa > quantile(rppa, 0.6) ~ "L",
          rppa < quantile(rppa, 0.4) ~ "H",
          TRUE ~ "M"
        )) %>%
        dplyr::filter(type != "M") %>%
        t.test(score~type, data = .) %>%
        broom::tidy() %>% dplyr::select(estimate, p.value) -> .t_pval
      
      # t.test(score~type, data = .x) %>% 
        # broom::tidy() %>% dplyr::select(estimate, p.value) -> .t_pval
    }
  )) %>% 
  tidyr::unnest(pval) %>% 
  dplyr::mutate(cancer_types = reorder(cancer_types, estimate, sort)) -> p62_csc_pval 

CPCOLS <- c("#CD0000","#000080")

p62_csc_pval$cancer_types %>% levels() -> levs

p62_csc_pval %>% 
  dplyr::mutate(mark = dplyr::case_when(
    dplyr::between(p.value, 0.001, 0.05) ~ "*",
    dplyr::between(p.value, 0.0001, 0.001) ~ "**",
    p.value < 0.0001 ~ "***",
    TRUE ~ ""
  )) %>% 
  dplyr::select(-merge, -p.value, -estimate) %>% 
  tibble::deframe() -> markers

p62_csc_pval %>%  
  tidyr::unnest() %>% 
  ggplot(aes(x = cancer_types, y = score, fill = type)) +
  geom_boxplot(outlier.size = NA) +
  scale_fill_manual(
    name = "Group", 
    values = CPCOLS,
    labels = c("High Autophagy", "Low Autophagy")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(x = "Cancer Types", y = "Stemness score")  +
  annotate(geom = "text", x = seq_along(levs), y = 2.3, label = markers[levs]) -> stemness_pattern
  
ggsave(filename = "fig_01_stemness_pattern.pdf", 
       plot = stemness_pattern, 
       width = 8,
       height = 6, 
       path = stemness_path)
# correlation -------------------------------------------------------------

p62_csc %>% 
  dplyr::mutate(corr = purrr::map2(
    .x = cancer_types,
    .y = merge,
    .f = function(.x, .y) {
      cor.test(.y$rppa, .y$score, method = 'spearman') %>% 
        broom::tidy() %>% 
        dplyr::select(ceof = estimate, pval = p.value) -> .corr
      .label <- glue::glue("Rs = {signif(.corr[[1]], 3)}
                           p = {signif(.corr[[2]],3)}")
      .y %>% 
        ggplot(aes(x = rank(rppa), y = rank(score))) +
        geom_point() + 
        geom_smooth(se = F, method = "lm") +
        annotate(geom = "text",  x = 10, y = 10, label = .label) +
        ggthemes::theme_gdocs() +
        labs(title = .x, x = "", y = "") +
        theme(plot.title = element_text(hjust = 0.5)) -> .p
      tibble::tibble(
        coef = signif(.corr[[1]],3), 
        pval = signif(.corr[[2]],4), 
        p = list(.p))
    }
  )) %>% 
  dplyr::select(-merge) %>% 
  tidyr::unnest(corr) %>% 
  dplyr::arrange(coef) -> p62_csc_corr

p62_csc_corr %>% dplyr::select(-p) %>% dplyr::arrange(-coef) %>% knitr::kable()
  

# stemness marker genes ---------------------------------------------------
# cd44+, cd24-, cd133+, cd90+, aldh1a1+, FLOT2+, abcb5+, EpCAM+,OCT4, NANOG, and SOX2,FOXO1, FOXO3, FOXO4, and FOXO6
stemness_genes <- c("CD24", "CD44", "PROM1", "THY1", "ALDH1A1", "ABCB5", "EPCAM", "OCT4", "NANOG", "SOX2", "FOXO1", "FOXO3", "FOXO4", "FOXO6")

expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))

expr %>% 
  dplyr::inner_join(p62_sample_classification, by = "cancer_types") %>%  
  dplyr::mutate(expr = purrr::map(
    .x = expr,
    .f = function(.x) {
      .x %>% 
        dplyr::filter(symbol %in% stemness_genes) %>%
        dplyr::select(-entrez_id) %>% 
        tidyr::drop_na() %>% 
        dplyr::filter(symbol != "?") %>% 
        tidyr::gather(key = barcode, value = expr, -symbol) %>% 
        dplyr::filter(stringr::str_sub(barcode, start = 14, end = 15) == "01") %>% 
        dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 12)) %>% 
        dplyr::select(-barcode) %>% 
        dplyr::distinct(symbol, sample, .keep_all = T) 
    }
  )) %>% 
  dplyr::mutate(merge = purrr::map2(
    .x = expr,
    .y = data,
    .f = function(.x, .y) {
      .y %>% 
        dplyr::mutate(type = ifelse(rppa > median(rppa), "L", "H")) %>% 
        dplyr::select(sample, type) %>% 
        dplyr::inner_join(.x, by = "sample")
    }
  )) -> p62_expr_merge

p62_expr_merge %>% 
  dplyr::select(-expr, -data) %>% 
  dplyr::mutate(pval = purrr::map(
    .x = merge,
    .f = function(.x) {
      .x %>% 
        dplyr::mutate(expr = log2(expr + 0.1)) %>% 
        dplyr::group_by(symbol) %>% 
        dplyr::do(
          broom::tidy(
          tryCatch(
            t.test(expr ~ type, data = .) ,
            error = function(e){1}
          ))
        ) %>% 
        dplyr::select(symbol, fc = estimate, pval = p.value) %>% 
        dplyr::ungroup()
    }
  )) %>% 
  tidyr::unnest(pval) %>% 
  dplyr::filter(pval < 0.05, abs(fc) > log2(1.1)) %>% 
  dplyr::mutate(fc = ifelse(fc < -2, -2, fc)) %>% 
  dplyr::mutate(fc = ifelse(fc > 2, 2, fc)) -> plot_ready

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = mean(fc)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(m = mean(fc)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(symbol) -> gene_rank

plot_ready %>% 
  ggplot(aes(x = cancer_types, y = symbol)) +
  geom_point(aes(color = fc, size = -log10(pval))) +
  scale_color_gradient2(
    high = "red",
    mid = "white", 
    low = "blue",
    name = "log2(FC)",
    guide = F
  ) +
  scale_x_discrete(limit = levs) +
  scale_y_discrete(limit = gene_rank) +
  scale_size(name = "p-value", guide = F) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(color = c("red", "red", "red", "red", "red", "black", "black","black", "black","black", "red", "red")),
    legend.position = "bottom"
  ) +
  labs(x = "Cancer Types", y = "Symbol") -> p

ggsave(filename = "fig_02_csc_signature_gene_diff.pdf", 
       plot = p, 
       width = 8,
       height = 2.5, 
       path = stemness_path)

# save --------------------------------------------------------------------


save.image(file.path(stemness_path, ".rda_03_stemness.rda"))
