
# library -----------------------------------------------------------------
library(magrittr)
library(ggplot2)

# load path --------------------------------------------------------------------
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
rppa_path <- file.path(expr_path, "03_e_rppa")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- file.path(expr_path, "03_a_gene_expr")

# load data ---------------------------------------------------------------
rppa_expr <- readr::read_rds(file.path(tcga_path, "pancan32_rppa_expr.rds.gz"))
rppa_name <- readr::read_rds(file.path(tcga_path, "rppa_name_symbol.rds.gz"))
expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))
gene_list <- readr::read_rds(file.path(expr_path, "rds_03_a_atg_lys_gene_list.rds.gz"))
rppa_name %>% dplyr::semi_join(gene_list, by = "symbol") -> atg_rppa

atg_rppa %>% 
  dplyr::inner_join(gene_list, by = "symbol") %>% 
  dplyr::select(1,2, process) %>% 
  dplyr::arrange(process) -> sym_func
knitr::kable(sym_func)

fn_merge <- function(.x, .y, atg_rppa){
  # .x <- .te$protein_expr[[1]]
  # .y <- .te$filter_expr[[1]]
  
  .x %>% 
    dplyr::select(-symbol) %>% 
    tidyr::gather(key = barcode, value = rppa, -protein) %>% 
    dplyr::inner_join(atg_rppa, by = "protein") %>% 
    dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 16)) %>% 
    dplyr::select(-barcode) -> .dx
  
  .y %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::inner_join(atg_rppa, by = "symbol") %>% 
    dplyr::mutate(sample = stringr::str_sub(barcode, start = 1, end = 16)) %>% 
    dplyr::select(-barcode) %>% 
    dplyr::mutate(expr = log2(expr + 0.01)) -> .dy
  
  .dx %>% dplyr::inner_join(.dy, by = c("symbol", "sample", "protein")) -> .d 
  
  # .d %>% dplyr::filter(symbol == "SQSTM1")
}

rppa_expr %>% 
  dplyr::inner_join(expr, by = "cancer_types") %>% 
  dplyr::mutate(merge = purrr::map2(.x = protein_expr, .y = filter_expr, .f = fn_merge, atg_rppa)) %>% 
  tidyr::unnest(merge) %>% 
  tidyr::drop_na() -> rppa_expr

rppa_expr %>% dplyr::filter(symbol == "SQSTM1") -> p62

p62 %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::mutate(rank = rank(rppa)) -> p62_d
p62_d %>% 
  dplyr::summarise(m = mean(rppa)) %>% 
  dplyr::arrange(m) %>% 
  dplyr::pull(cancer_types) -> lev

p62_d %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(cancer_types = factor(cancer_types, levels = lev)) %>% 
  ggplot(aes(x = rank, y = rppa)) +
  geom_point() +
  facet_grid(~cancer_types, scales = "free", switch = "x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5),
    
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = "black")
  ) +
  labs(x = "", y = "p62 RPPA", title = "p62 rppa distribution") -> p62_rppa
ggsave(filename = file.path(rppa_path, "BECN1_rppa_distribution.pdf"), plot = p62_rppa, device = "pdf",
       width = 14)

p62 %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::filter(n() > 10) %>% 
  dplyr::do(
    broom::tidy(
      tryCatch(
        expr = cor.test(~rppa + expr, data = ., method = "pearson")
      )
    )
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cancer_types, coef = estimate, pval = p.value) %>% 
  dplyr::arrange(dplyr::desc(coef)) %>% 
  dplyr::filter(pval < 0.05) %>% 
  dplyr::mutate(pval = -log10(pval)) %>% 
  dplyr::mutate(pval = ifelse(pval > 20, 20, pval)) %>% 
  ggplot(aes(x = coef, y = pval)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = cancer_types)) +
  labs(x = "Coef", y = "P-value", title = "p62 mRNA correlates with rppa") +
  theme_bw() -> pval_coef

ggsave(filename = file.path(rppa_path, "p62_rppa_pval_coef.pdf"), plot = pval_coef, device = "pdf")

# survival ----------------------------------------------------------------


clinical <- readr::read_rds(file.path(tcga_path, "pancan34_clinical.rds.gz"))

clinical %>% 
  # dplyr::filter(cancer_types == "BRCA") %>%
  dplyr::mutate(clinical = purrr::map(.x = clinical , dplyr::select, barcode, os_days, os_status)) %>% 
  tidyr::unnest() %>% 
  dplyr::select(cancer_types, barcode, time = os_days, status = os_status) %>% 
  tidyr::drop_na() %>% 
  dplyr::filter(time > 0) %>% 
  dplyr::mutate(status = plyr::revalue(status, c("Dead" = 1, "Alive" = 0))) %>% 
  dplyr::mutate(status = as.integer(status)) -> brca


rppa_expr %>% 
  dplyr::filter(symbol %in% c("BECN1", "SQSTM1","BCL2")) %>% 
  dplyr::select(-expr, -protein) %>% dplyr::distinct(cancer_types, symbol, sample, .keep_all = T) %>% 
  tidyr::spread(symbol, rppa) %>% 
  tidyr::drop_na() %>% 
  # dplyr::select(BCL2, BECN1, SQSTM1) %>% cor()
  dplyr::mutate(ind = SQSTM1 / BECN1) %>% 
  # dplyr::filter(cancer_types == "BRCA") %>% 
  dplyr::mutate(sample = stringr::str_sub(sample, 1, 12)) %>% 
  dplyr::rename(barcode = sample) %>% 
  dplyr::distinct(barcode, .keep_all = T) %>% 
  dplyr::inner_join(brca, by = c("cancer_types", "barcode")) %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::mutate(group = ifelse(ind > 0, "high", "low")) %>% 
  dplyr::ungroup() -> rppa_expr_d


rppa_expr_d %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::filter(n() > 10) %>% 
  dplyr::do(
    broom::tidy(
      survival::coxph(survival::Surv(time, status) ~ BECN1, data = .)
      )
    ) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(p.value < 0.05) %>% 
  dplyr::arrange(p.value) -> ind_coxp

rppa_expr_d %>% dplyr::filter(cancer_types == "LGG") -> .d
fit_x <- survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude)
survminer::ggsurvplot(fit_x, data = .d, pval = T, pval.method = T,
                      # title = paste(paste(cancer_types, gene, sep = "-"), "Coxph =", signif(p.value, 2)),
                      xlab = "Survival in days",
                      ylab = 'Probability of survival')

