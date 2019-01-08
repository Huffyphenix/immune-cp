
# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# load path ---------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/03_a_gene_expr"

# load data ---------------------------------------------------------------

rppa_expr <- readr::read_rds(file.path(tcga_path, "pancan33_rppa_expr_l4.rds.gz")) %>% 
  dplyr::filter(cancer_types != "PANCAN19")
rppa_name <- readr::read_rds(file.path(tcga_path, "rppa_name_symbol.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path, "rds_03_a_atg_lys_gene_list.rds.gz"))
pcc <- readr::read_tsv(file = file.path(tcga_path, "PanCanAtlasTumors_color_coded_by_organ_system_20170302.tsv"))

rppa_name %>% dplyr::semi_join(gene_list, by = "symbol") -> atg_rppa

atg_rppa %>% 
  dplyr::inner_join(gene_list, by = "symbol") %>% 
  dplyr::select(1,2, process) %>% 
  dplyr::arrange(process) -> sym_func
knitr::kable(sym_func)

PI3K_AKT <- "pS473|pT308|pS9|pS21S9|pT246|pT1462" %>% stringr::str_replace_all(pattern = "\\|", replacement = "$\\|") # AKT pS473, AKT pT308, GSK3 pS9,  GSK3 pS21-pS9, PRAS40 pT246 TSC2 pT1462

# mtor activity score -----------------------------------------------------
# anotation from cancer cell paper 10.1016/j.ccell.2017.04.013
mTOR <- "pS2448|pT1135|pS65|pT37T46|pT70|pT389|pS235S236|pS240S244" %>% stringr::str_replace_all(pattern = "\\|", replacement = "$\\|") # mTOR pS2448, RICTOR pT1135, 4EBP1 pS65, 4EBP1 pT37T46, 4EBP1 pT70, S6K pT389, S6 pS235S236, S6 pS240S244

rppa_expr %>% 
  dplyr::mutate(mtor = purrr::map(
    .x = protein_expr,
    .f = function(.x) {
      .x %>% dplyr::filter(
        stringr::str_detect(protein, pattern = mTOR)
      ) %>% 
        dplyr::select(-symbol) %>% 
        dplyr::mutate(direc = as.numeric(
          plyr::revalue(protein, c("MTOR_pS2448" = 1, "P70S6K_pT389" = 1, "RICTOR_pT1135" = 1, "S6_pS235S236" = 1, "S6_pS240S244" = 1, "X4EBP1_pS65" = 1, "X4EBP1_pT37T46" = 1, "X4EBP1_pT70" = 1)))) %>% 
        tidyr::gather(key = barcode, value = rppa, -protein, -direc) %>% 
        dplyr::group_by(barcode) %>% 
        dplyr::summarise(mtor_score = sum(rppa * direc))
    }
  )) %>% 
  dplyr::mutate(pi3k_akt = purrr::map(
    .x = protein_expr,
    .f = function(.x) {
      .x %>% dplyr::filter(
        stringr::str_detect(protein, pattern = PI3K_AKT)
      ) %>% 
        dplyr::select(-symbol) %>% 
        dplyr::mutate(direc = as.numeric(
          plyr::revalue(protein, c("AKT_pS473" = 1, "AKT_pT308" = 1, "GSK3_pS9" = 1, "GSK3ALPHABETA_pS21S9" = 1, "PRAS40_pT246" = 1, "TUBERIN_pT1462" = 1)))) %>% 
        tidyr::gather(key = barcode, value = rppa, -protein, -direc) %>% 
        dplyr::group_by(barcode) %>% 
        dplyr::summarise(pik_score = sum(rppa * direc)) 
    }
  )) %>% 
  dplyr::select(-protein_expr) %>% 
  dplyr::mutate(merge = purrr::map2(
    .x = mtor, 
    .y = pi3k_akt,
    .f = function(.x, .y){
      .x %>% dplyr::inner_join(.y, by = 'barcode')
    }
  )) %>% 
  dplyr::select(cancer_types, merge) %>% 
  tidyr::unnest() -> mtor_pi3k_score 

glue::glue("{PI3K_AKT}$|{mTOR}$|P62|BECLIN|FOXO3A|RAPTOR|AMPKALPHA_pT172|AMPKALPHA|MTOR") -> PI3K_AKT_mTOR

rppa_expr %>% 
  dplyr::mutate(p62 = purrr::map(
    .x = protein_expr,
    .f = function(.x) {
      .x %>% 
        dplyr::filter(
          stringr::str_detect(protein, pattern = PI3K_AKT_mTOR)
        ) %>% 
        dplyr::select(-symbol) %>% 
        tidyr::gather(key = barcode, value = rppa, -protein) %>% 
        tidyr::spread(key = protein, value = rppa)
    }
  )) %>% 
  dplyr::select(-protein_expr) %>% 
  tidyr::unnest() %>% 
  dplyr::inner_join(mtor_pi3k_score, by = c("cancer_types", "barcode")) -> p62_mtor_pi3k_score

names(p62_mtor_pi3k_score)[-c(1,2,14)] -> ref_proteins

fn_dist <- function(.s){
  
  lev <- c("TGCT", "LGG", "GBM", "PRAD", "UCS", "PAAD", "PCPG", "THCA", "KICH",
           "MESO", "CHOL", "THYM", "SARC", "KIRP", "COAD", "BLCA", "READ", "UCEC",
           "STAD", "CESC", "BRCA", "KIRC", "HNSC", "ESCA", "SKCM", "LUAD", "OV", 
           "LUSC", "UVM", "ACC", "LIHC", "DLBC")
  
  p62_mtor_pi3k_score %>% 
    dplyr::select(1,2, .s) %>% 
    dplyr::mutate(cancer_types = factor(cancer_types, levels = lev)) %>% 
    dplyr::rename(p62 = rlang::UQ(.s)) %>% 
    dplyr::mutate(p62 = dplyr::case_when(p62 > 3 ~ 3, p62 < -2 ~ -2, TRUE ~ p62)) %>% 
    ggplot(aes(x = cancer_types, y = p62)) +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(
      aes(color = cancer_types), 
      position = position_jitter(width = 0.05), 
      alpha = 0.4,
      size = 0.8) +
    scale_color_manual(
      name = "Cancer Types",
      values = dplyr::slice(pcc, match(lev, cancer_types)) %>% dplyr::pull(color)
    ) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
      
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing.x = unit(0, "lines")
    ) +
    guides(color = F) +
    labs(
      x = "", 
      y = glue::glue("{.s} RPPA score"), 
      title = "") -> .dist
  
  ggsave(filename = file.path(afhl_class, glue::glue("fig_01_{.s}_rppa_distribution_across_cancer_types.pdf")), 
         plot = .dist,
         device = "pdf", 
         width = 7, 
         height = 3)
}

fn_p62_corr_all <- function(.s) {
  p62_mtor_pi3k_score %>% 
    dplyr::select(1,2, .s, p62 = P62LCKLIGAND) %>% 
    # dplyr::group_by(cancer_types) %>% 
    dplyr::do(
      cor.test(dplyr::pull(., p62), dplyr::pull(., .s), method = "spearman") %>% 
        broom::tidy()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(coef = estimate, pval = p.value) %>% 
    tibble::add_column(vs_type = glue::glue("p62_vs_{.s}"), .before = 1)
}
  
ref_proteins %>% purrr::walk(.f = fn_dist)
ref_proteins %>% 
  purrr::map(.f = fn_p62_corr_all) %>% 
  dplyr::bind_rows() %>% 
  dplyr::mutate(vs_type = stringr::str_replace(vs_type, "p62_vs_", "")) %>% 
  ggplot(aes(x = coef, -log10(pval))) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = vs_type)) +
  ggthemes::theme_gdocs() +
  labs(x = "Coef", y = "P-value") -> p62_corr_all_cancer_types

  ggsave(filename = "fig_01_across_cancer_type_p62_corr.pdf",
         device = "pdf",
         plot = p62_corr_all_cancer_types,
         path = afhl_class,
         width = 7,
         height = 6)

fn_p62_corr <- function(.s) {
    p62_mtor_pi3k_score %>% 
      # dplyr::mutate(ratio = exp(BECLIN) / exp(P62LCKLIGAND)) %>%
      dplyr::select(1,2, .s, p62 = P62LCKLIGAND) %>% 
      dplyr::group_by(cancer_types) %>% 
      dplyr::do(
        cor.test(dplyr::pull(., p62), dplyr::pull(., .s), method = "spearman") %>% 
          broom::tidy()
      ) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(cancer_types, coef = estimate, pval = p.value) %>% 
      tibble::add_column(vs_type = glue::glue("p62_vs_{.s}"), .before = 2)
  }
  
ref_proteins %>% purrr::map(.f = fn_p62_corr) %>% dplyr::bind_rows() -> p62_corr

p62_corr %>% 
  dplyr::filter(!vs_type %in% c("p62_vs_mtor_score", "p62_vs_pik_score")) %>%
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = mean(coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) -> cancer_rank

p62_corr %>% 
  dplyr::filter(!vs_type %in% c("p62_vs_mtor_score", "p62_vs_pik_score")) %>%
  dplyr::group_by(vs_type) %>% 
  dplyr::summarise(m = mean(coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(vs_type) -> vs_type_rank

p62_corr %>% 
  # dplyr::filter(!vs_type %in% c("p62_vs_mtor_score", "p62_vs_pik_score")) %>%
  dplyr::filter(pval < 0.05) %>% 
  dplyr::mutate(pval = -log10(pval)) %>% 
  dplyr::mutate(pval = ifelse(pval > 30, 30, pval)) %>% 
  dplyr::mutate(coef = dplyr::case_when(
    coef > 0.5 ~ 0.5,
    coef < -0.5 ~ -0.5,
    TRUE ~ coef
  )) %>%
  ggplot(aes(x = cancer_types, y = vs_type)) +
  geom_point(aes(color = coef, size = pval)) +
  scale_color_gradient2(
    name = "Coef",
    high = "red",
    mid = "white",
    low = "blue",
    breaks = seq(-0.5,0.5,length.out = 5),
    labels = c("=<-0.5", "-0.25", "0", "0.25", ">=0.5")
  ) +
  scale_size(name = "P-value") +
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = c(vs_type_rank),
    labels = c(vs_type_rank %>% 
      stringr::str_replace( pattern = "p62_vs_", "") %>% 
      stringr::str_replace(pattern = "_", replacement = " ") %>% 
        stringr::str_replace(pattern = "ALPHA", "") %>% 
        stringr::str_replace(pattern = "BETA", "") %>% 
        stringr::str_replace(pattern = "^X", ""),
      "","mTOR score", "PI3K/AKT score"),
    expand = c(0.05,0.05)) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.text.y = element_text(
      
    )
  ) +
  labs(x = "Cancer Types", y = "Proteins") -> p62_corr_protein

ggsave(filename = "fig_01_p62_corr_proteins.pdf",
       plot = p62_corr_protein, 
       path = afhl_class,
       device = "pdf", 
       width = 9,
       height = 6)

p62_mtor_pi3k_score %>% 
  dplyr::select(cancer_types, P62LCKLIGAND, BECLIN) %>% 
  dplyr::left_join(coef_pval, by = "cancer_types") -> p62_mtor_pi3k_score_mid
p62_mtor_pi3k_score_mid %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(m = mean(coef)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) -> levs

p62_mtor_pi3k_score_mid %>% 
  dplyr::mutate(cancer_types = factor(cancer_types, levels = levs)) %>% 
  dplyr::mutate(coef = paste("coef", signif(coef, 3), sep = " = ")) %>% 
  dplyr::mutate(pval = paste("p",signif(pval, 3), sep = " = ") ) %>% 
  ggplot(aes(x = rank(P62LCKLIGAND), y = rank(BECLIN))) +
  geom_point() +
  geom_smooth(se = F, method = "lm") +
  facet_wrap(~cancer_types + coef + pval, nrow = 5) +
  labs(x = "p62", y = "BECLIN") +
  theme(
    strip.background = element_rect(fill = 'white', color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black")
  ) -> p
ggsave(filename = "fig_02_lm_rppa_p62_becn1.pdf", 
       plot = p,
       device = "pdf",
       width = 20, 
       height = 20,
       path = afhl_class)


p62_mtor_pi3k_score %>% 
  dplyr::group_by(cancer_types) %>% 
  purrrlyr::by_slice(
    ..f = function(.d){
      tryCatch(
      cor.test(~P62LCKLIGAND + BECLIN, data = .d) %>% 
        broom::tidy() %>% 
        dplyr::select(coef = estimate, pval = p.value) -> .pval,
      error = function(e) return(NULL))
      .cancer <- .d$cancer_types
      glue::glue("rppa spearman correlation in {.cancer}") -> .title
      glue::glue("Rs = {signif(.pval$coef, 3)} 
                 p = {signif(.pval$pval,3)}") -> .label
      print(.d$cancer_types)
      .d %>% 
        ggplot(aes(x = rank(P62LCKLIGAND), y = rank(BECLIN))) +
        geom_point() +
        geom_smooth(se = F, method = "lm") +
        ggthemes::theme_gdocs() +
        annotate(geom = "text", label = .label, x = 2, y = 2) +
        labs(x = "p62", y = "BECLIN") -> .p
      tibble::tibble(coef = .pval$coef, pval = .pval$pval, plot = list(.p))
    }
  ) -> p62_becn_plot
p62_becn_plot %>% tidyr::unnest() %>% dplyr::select(-plot) -> coef_pval
gridExtra::arrangeGrob(grobs = p62_becn_plot$.out, ncol = 5) %>% 
  ggsave(filename = "fig_02_lm_rppa_p62_becn1.pdf", 
         plot = .,
         device = "pdf",
         width = 20, 
         height = 25,
         path = afhl_class)

# mTORC1 complex MTOR, RPTOR, MLST8, DEPTOR -----
# phosphorylated (actiavated)
p_mTOR <- c("MTOR_pS2448")

unp_mTOR <- c("MTOR", "RAPTOR")

p_ampk <- c("AMPKALPHA_pT172")

unp_ampk <- c("AMPKALPHA")

sym_names <- c("BECLIN" = "BECN1", "MTOR_pS2448" = "MTOR", "P62LCKLIGAND" = "SQSTM1",  "AMPKALPHA" = "PRKAA1", "AMPKALPHA_pT172" = "PRKAA1")

.names <- c("AMPKALPHA" = "PRKAA1", "AMPKALPHA_pT172" = "PRKAA1 pT172", "BECLIN" = "BECN1", "MTOR_pS2448" = "MTOR pS2448", "P62LCKLIGAND" = "p62", "MTOR" = "MTOR", "BCL2" = "BCL2", "RAPTOR" = "RAPTOR")

# mTOR and P62, BECN1 and BCL2---------------------------------------------------------
fn_one_corr <- function(.z, .d){
  .d %>% dplyr::select(.z) -> .dd
  cor.test(.dd[[1]], .dd[[2]]) %>% broom::tidy() -> .corr
  tibble::tibble(p1 = .z[1], p2 = .z[2], coef = .corr$estimate, pval = .corr$p.value)
}

fn_all_corr <- function(.d) {
  .protein <- colnames(.d)
  combn(.protein, 2) -> .combn
  as.list(data.frame(.combn)) %>% purrr::map(.f = fn_one_corr, .d)
}

fn_select_marker <- function(.x){
  
  .names <- c("AMPKALPHA" = "PRKAA1", "AMPKALPHA_pT172" = "PRKAA1 pT172", "BECLIN" = "BECN1", "MTOR_pS2448" = "MTOR pS2448", "P62LCKLIGAND" = "p62")
  
  .x %>% 
    dplyr::filter(protein %in% c(p_mTOR, unp_mTOR, p_ampk, unp_ampk, "P62LCKLIGAND", "BECLIN", "BCL2")) %>% 
    dplyr::mutate(protein = plyr::revalue(x = protein, replace = .names)) %>% 
    dplyr::select(-symbol)
  }

rppa_expr %>% 
  dplyr::mutate(protein_expr = purrr::map(.x = protein_expr, .f = fn_select_marker)) -> 
  atg_rppa_expr

fn_protein_dist <- function(.gene) {

  atg_rppa_expr %>% 
    dplyr::mutate(
      protein_expr = purrr::map(
        .x = protein_expr, 
        .f = function(.x) dplyr::filter(.x, protein == .gene) %>% 
          tidyr::gather(key = barcode, value = rppa, -protein) %>% 
          tidyr::spread(key = protein, value = rppa)
        )
      ) %>% 
    tidyr::unnest(protein_expr) %>% 
    dplyr::rename(p62 = rlang::UQ(.gene)) %>% 
    dplyr::group_by(cancer_types) %>% 
    dplyr::mutate(rank = rank(p62)) %>% 
    dplyr::ungroup() -> p62
  
  p62 %>% 
    dplyr::group_by(cancer_types) %>% 
    dplyr::summarise(m = median(p62)) %>% 
    dplyr::arrange(m) %>% 
    dplyr::pull(cancer_types) -> lev
  
  lev <- c("TGCT", "LGG", "GBM", "PRAD", "UCS", "PAAD", "PCPG", "THCA", "KICH",
           "MESO", "CHOL", "THYM", "SARC", "KIRP", "COAD", "BLCA", "READ", "UCEC",
           "STAD", "CESC", "BRCA", "KIRC", "HNSC", "ESCA", "SKCM", "LUAD", "OV", 
           "LUSC", "UVM", "ACC", "LIHC", "DLBC")
  
  p62 %>% 
    dplyr::mutate(cancer_types = factor(cancer_types, levels = lev)) %>% 
    dplyr::mutate(p62 = dplyr::case_when(p62 > 3 ~ 3, p62 < -2 ~ -2, TRUE ~ p62)) %>% 
    ggplot(aes(x = cancer_types, y = p62)) +
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(
      aes(color = cancer_types), 
      position = position_jitter(width = 0.05), 
      alpha = 0.4,
      size = 0.8) +
    scale_color_manual(
      name = "Cancer Types",
      values = dplyr::slice(pcc, match(lev, cancer_types)) %>% dplyr::pull(color)
    ) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
      
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.spacing.x = unit(0, "lines")
    ) +
    guides(color = F) +
    labs(
      x = "", 
      y = glue::glue("{.gene} RPPA score"), 
      title = "") -> plot_p62_rppa_dist
  
  ggsave(filename = file.path(afhl_class, glue::glue("fig_01_{.gene}_rppa_distribution_across_cancer_types.pdf")), 
         plot = plot_p62_rppa_dist,
         device = "pdf", 
         width = 7, 
         height = 3)
  
  plot_p62_rppa_dist
}

.names %>% purrr::map(fn_protein_dist) -> rppa_plot_dist
gridExtra::arrangeGrob %>% purrr::lift_dl(ncol = 1) %>% purrr::invoke(rppa_plot_dist) -> rppa_plot_dist_all
ggsave(filename = file.path(afhl_class, "fig_01_rppa_plot_dist_all.pdf"), plot = rppa_plot_dist_all, device = "pdf", width = 8, height = 15)

gridExtra::arrangeGrob(grobs = rppa_plot_dist[c(5, 3)]) %>% 
  ggsave(filename = "fig_01_rppa_p62_becn1.pdf", 
         plot = .,
         device = "pdf",
         width = 7, 
         height = 6)
#

# correlation with known proteins -------------------------------------

fn_corr_protein <- function(.x, .y){
  .corr_path <- file.path(afhl_class, "01_p62_mtor_ampk")
  if (!dir.exists(.corr_path)) dir.create(.corr_path)

  .x %>% 
    tidyr::gather(key = barcode, value = rppa, -protein) %>% 
    tidyr::spread(key = protein, value = rppa) %>% 
    dplyr::select(-barcode) -> .d
  
  .d %>% cor() -> .corr_matrix
  
  .title = glue::glue("Correlation of P62, BECLIN, AMPK and mTOR in {.y}")
  .filename = glue::glue("{.y}_p62_beclin_bcl2_mtor_ampk.pdf")
  
  ggcorrplot::ggcorrplot(.corr_matrix, hc.order = TRUE, outline.col = "white") +
    labs(title = .title) -> .p
  ggsave(filename = .filename, plot = .p, device = "pdf", path = .corr_path, width = 7, height = 7)
  
  list(
    p62_vs_BECN1 = cor.test(~p62 + BECN1, data = .d) %>% broom::tidy(),
    p62_vs_MTOR = cor.test(~p62 + MTOR, data = .d) %>% broom::tidy(),
    p62_vs_MTOR_pS2448 = cor.test(~p62 + `MTOR pS2448`, data = .d) %>% broom::tidy(),
    p62_vs_BCL2 = cor.test(~p62 + BCL2, data = .d) %>% broom::tidy(),
    BECN1_vs_BCL2 = cor.test(~BECN1 + BCL2, data = .d) %>% broom::tidy(),
    MTOR_vs_MTOR_pS2448 = cor.test(~MTOR + `MTOR pS2448`, data = .d) %>% broom::tidy()
  ) %>% 
    tibble::enframe() %>% 
    tidyr::unnest() %>% 
    dplyr::select(name, coef = estimate, pval = p.value) -> .d_corr
  
  .lm_dir <- file.path(.corr_path, "01_linear_regression")
  if (!dir.exists(.lm_dir)) dir.create(.lm_dir)
  
  .d_corr %>% 
    dplyr::filter(abs(coef) > 0.3, pval < 0.05) %>%
    dplyr::mutate(coef = signif(coef, digits = 4), pval = signif(pval, digits = 4)) %>% 
    purrr::pmap(
      .f = function(name, coef, pval) {
        
        .xy <- stringr::str_split(name, "_vs_", simplify = T) %>% 
          as.vector() %>% 
          stringr::str_replace("_", " ")
        
        .d %>% 
          ggplot(aes(x = rlang::eval_bare(rlang::sym(.xy[1])), y = rlang::eval_bare(rlang::sym(.xy[2])))) + 
          geom_point() + 
          annotate("text", x = 0, y = 0, 
                   label = glue::glue("R = {coef}
                                      p-value = {pval}")) +
          geom_smooth(se = F, method = "lm") +
          ggthemes::theme_gdocs() +
          labs(x = .xy[1], y = .xy[2]) -> .p
        
        ggsave(
            filename = glue::glue("{.y}_{paste0(.xy, collapse = '_')}_{coef}_{pval}.pdf"),
            plot = .p, device = "pdf", path = .lm_dir, 
            width = 5, height = 5)
    }) # draw linear plot

  .d_corr
}

atg_rppa_expr %>% 
  dplyr::mutate(corr = purrr::map2(.y = cancer_types, .x = protein_expr, .f = fn_corr_protein)) ->
  atg_rppa_expr


atg_rppa_expr %>% 
  tidyr::unnest(corr) %>% 
  dplyr::mutate(
    coef = 
      dplyr::case_when(
        coef < -0.4 ~ -0.4,
        coef > 0.5 ~ 0.5,
        TRUE ~ coef
      )) %>% 
  # dplyr::bind_rows(p62_mtor_pi3k_score_coef) %>% 
  dplyr::mutate(sig = ifelse(pval < 0.05, "*", "")) %>%
  dplyr::mutate(h_corr = ifelse(abs(coef) > 0.3, "H", "L")) -> 
  corr_plot_ready

cancer_rank <- 
  corr_plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(s = sum(coef)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(s) %>% 
  dplyr::pull(cancer_types)

pair_rank <- 
  corr_plot_ready %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(s = sum(coef)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(s) %>% 
  dplyr::pull(name)

CPCOLS <- c("#00008B", "#FCFCFC", "#CD0000")

corr_plot_ready %>% 
  ggplot(aes(x = cancer_types, y = name, fill = coef)) +
  geom_tile(aes(color = h_corr), width = 0.85, height = 0.95, size = 1) +
  geom_text(aes(label = sig), color = "white") + 
  # scale_color_manual(values = "black", guide = F) +
  scale_fill_gradient2(
    low = CPCOLS[1],
    mid = CPCOLS[2],
    high = CPCOLS[3],
    limit = c(-0.4, 0.5),
    breaks = c(-0.4, -0.2, 0, 0.3, 0.5)
  ) +
  scale_x_discrete(limit = cancer_rank) +
  scale_y_discrete(limit = pair_rank, labels = stringr::str_replace_all(
    pair_rank, "_vs_", " vs "
  )) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    
    panel.background = element_rect(fill = "white"),
    
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(
    title = "Coef", 
    title.position = "left",
    title.vjust = 0.2,
    title.hjust = 0.5,
    
    keywidth = 1,
    keyheight = 0.8,
    label.position = "top",
    label.hjust = 1,
    direction = "horizontal")) +
  labs(x = "", y = "") -> corr_tile_plot
ggsave(filename = "fig_01_corr_tile_plot.pdf", plot = corr_tile_plot, device = "pdf", path = afhl_class, height = 4, width = 8)

# p62 rppa and mrna correlation -------------------------------------------
mrna_expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz")) 

mrna_expr %>% 
  dplyr::mutate(expr = purrr::map(
  .x = expr,
  .f = function(.x) {
    .x %>% dplyr::filter(symbol == "SQSTM1")
  }
)) -> SQSTM1_mrna
  
atg_rppa_expr %>% 
  dplyr::select(-corr) %>% 
  dplyr::mutate(protein_expr = purrr::map(
    .x = protein_expr,
    .f = function(.x) {
      .x %>% dplyr::filter(protein == "p62")
    }
  )) %>% 
  dplyr::inner_join(SQSTM1_mrna, by = "cancer_types") ->  p62_rppa_expr

fn_corr_mrna_protein_p62 <- function(.x, .y){
  .y %>% 
    dplyr::select(-entrez_id) %>% 
    dplyr::filter(symbol == "SQSTM1") %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 16)) %>% 
    dplyr::select(symbol, sample, expr) %>% 
    dplyr::distinct(sample, .keep_all = T) %>% 
    dplyr::mutate(expr = log2(expr + 0.1)) -> .my
  
  .x %>% 
    dplyr::filter(protein == "p62") %>% 
    tidyr::gather(key = barcode, value = rppa, -protein) %>% 
    dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 16)) %>% 
    dplyr::select(protein, sample, rppa) %>% 
    dplyr::distinct(sample, .keep_all = T) -> .px
  
  .px %>% 
    dplyr::inner_join(.my, by = "sample") %>% 
    tidyr::drop_na() %>% 
    cor.test(~rppa + expr, data = ., method = "spearman") %>% 
    broom::tidy() %>% 
    dplyr::select(coef = estimate, pval = p.value)
}

p62_rppa_expr %>% 
  dplyr::mutate(
    mp_p62 = purrr::map2(
    .x = protein_expr,
    .y = expr,
    .f = fn_corr_mrna_protein_p62
  )) %>% 
  dplyr::select(-protein_expr, -expr) %>% 
  tidyr::unnest(mp_p62) -> p62_corr_rppa_expr

CPCOLS <- c("#000080", "#FF4040")

p62_corr_rppa_expr %>% 
  dplyr::mutate(pval = -log10(pval)) %>% 
  dplyr::mutate(pval = ifelse(pval > 50, 50, pval)) %>% 
  dplyr::mutate(sig = ifelse(pval > -log10(0.05) & abs(coef) > 0.3, "sig", "non-sig")) %>% 
  ggplot(aes(x = coef, y = pval)) +
  geom_point(aes(color = sig), size = 1) +
  ggrepel::geom_text_repel(aes(label = cancer_types), size = 3) +
  scale_color_manual(values = CPCOLS, name = "Significance") +
  theme_classic() +
  labs(x = "Coefficient", y = "P-value") -> p62_mrna_rppa_coef
  
ggsave(
  filename = "fig_01_p62_mrna_rppa_corr.pdf",
  plot = p62_mrna_rppa_coef,
  device = "pdf",
  path = afhl_class,
  width = 6,
  height = 5
)

# clinical data and classify ----------------------------------------------
clinical <- readr::read_rds(file.path(tcga_path, "pancan34_clinical.rds.gz"))
p62_rppa_expr %>% 
  dplyr::select(-expr) %>% 
  dplyr::mutate(
    p62 = purrr::map(
      .x = protein_expr,
      .f = function(.x){
        .x %>% 
          tidyr::gather(key = barcode, value = rppa, -protein) %>% 
          dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 12)) %>% 
          dplyr::distinct(sample, .keep_all = T) %>% 
          dplyr::mutate(type = ifelse(rppa > median(rppa), "H", "L"))
      }
    )
  ) %>% 
  tidyr::unnest(p62) -> p62_sample_classification


p62_sample_classification %>% readr::write_rds(path = file.path(afhl_class, '.rds_01_p62_sample_classification.rds.gz'), compress = "gz")

save.image(file.path(afhl_class, ".rda_01_af_h_l_classification.rda.gz"))

load(file.path(afhl_class, ".rda_01_af_h_l_classification.rda.gz"))









