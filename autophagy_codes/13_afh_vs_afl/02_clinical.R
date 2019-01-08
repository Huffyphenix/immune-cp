
# library -----------------------------------------------------------------

library(magrittr)
library(ggplot2)

# load path ---------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
clinical_path <- file.path(afhl_path, "02_clinical")

# clinical data and classify ----------------------------------------------
clinical <- readr::read_rds(file.path(tcga_path, "pancan34_clinical.rds.gz"))
p62_sample_classification <- readr::read_rds(path = file.path(afhl_class, '.rds_01_p62_sample_classification.rds.gz')) %>% 
  dplyr::select(-barcode) %>% 
  tidyr::nest(-cancer_types)

clinical %>% 
  dplyr::mutate(
    clinical = purrr::map(
      .x = clinical,
      .f = function(.x){
        .x %>% 
          dplyr::select(barcode, time = os_days, status = os_status) %>% 
          dplyr::rename(sample = barcode) %>% 
          tidyr::drop_na() %>% 
          dplyr::filter(time > 0) %>% 
          dplyr::mutate(status = plyr::revalue(status, c("Dead" = 1, "Alive" = 0))) %>% 
          dplyr::mutate(status = as.integer(status))
      }
    )
  ) %>%
  dplyr::inner_join(p62_sample_classification, by = "cancer_types") %>% 
  dplyr::mutate(
    merge = purrr::map2(
      .x = clinical,
      .y = data,
      .f = function(.x, .y) {
        .x %>% dplyr::inner_join(.y, by = "sample") %>% 
          dplyr::rename(group = type)
      }
    )
  ) %>% 
  dplyr::select(cancer_types, merge) %>% 
  tidyr::unnest() -> p62_clinical

p62_clinical %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::filter(n() > 20) %>% 
  tidyr::nest(-cancer_types) -> p62_clinical_nest

p62_clinical_nest %>% 
  dplyr::mutate(hazard_ratio = purrr::map2(
    .x = data,
    .y = cancer_types,
    .f = function(.x, .y){
      .x %>% 
        dplyr::mutate(time = time / 30) %>% 
        dplyr::mutate(time = ifelse(time > 60, 60, time)) -> .dd
      
      .dd  -> .d
      
      survival::coxph(survival::Surv(time, status) ~ rppa, data = .dd) %>% 
        broom::tidy() %>% 
        dplyr::mutate(hazard_ratio = exp(estimate)) %>% 
        dplyr::select(hazard_ratio, coxp = p.value) -> .hazard_coxp
      
      
      survival::survdiff(survival::Surv(time, status) ~ group, data = .d) -> .d_diff
      kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
      
      survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude) -> fit_x
      
      CPCOLS <- c("#CD0000", "#000080")
      survminer::ggsurvplot(fit_x, 
                            data = .d, 
                            xlab = "Survival in months",
                            ylab = 'Probability of survival',
                            palette = CPCOLS,
                            ggtheme = theme_bw()) -> .p
      
      .label <- glue::glue("{.y}
                           p = {signif(kmp, 2)}")
      
      tibble::tibble(
        hazard_ratio = .hazard_coxp$hazard_ratio,
        coxp = .hazard_coxp$coxp,
        kmp = kmp,
        p = list(.p$plot + annotate("text", x = 5, y = 0.25, label = .label) + theme(legend.position = ""))
      )
    }
  )) %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest() -> p62_survival
p62_survival %>% dplyr::filter(coxp < 0.05, kmp > 0.05) -> p62_survival_median

# use top 40% and bottom 40%
p62_clinical_nest %>% 
  dplyr::mutate(hazard_ratio = purrr::map2(
    .x = data,
    .y = cancer_types,
    .f = function(.x, .y){
      .x %>% 
        dplyr::mutate(time = time / 30) %>% 
        dplyr::mutate(time = ifelse(time > 60, 60, time)) %>% 
        dplyr::mutate(
          group = dplyr::case_when(
            rppa > quantile(rppa, 0.5) ~ "L",
            rppa < quantile(rppa, 0.5) ~ "H",
            TRUE ~ "M"
          )
        ) -> .dd
      
      .dd %>% dplyr::filter(group != "M") -> .d
      
      survival::coxph(survival::Surv(time, status) ~ rppa, data = .d) -> .cox
        broom::tidy() %>% 
        dplyr::mutate(hazard_ratio = exp(estimate)) %>% 
        dplyr::select(hazard_ratio, coxp = p.value) -> .hazard_coxp
      
      
      survival::survdiff(survival::Surv(time, status) ~ group, data = .d) -> .d_diff
      kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
      
      survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude) -> fit_x
      
      CPCOLS <- c("#CD0000", "#000080")
      survminer::ggsurvplot(fit_x, 
                            data = .d, 
                            xlab = "Survival in months",
                            ylab = 'Probability of survival',
                            palette = CPCOLS,
                            ggtheme = theme_bw()) -> .p
      
      .label <- glue::glue("{.y}
                           p = {signif(kmp, 2)}")
      
      tibble::tibble(
        hazard_ratio = .hazard_coxp$hazard_ratio,
        coxp = .hazard_coxp$coxp,
        kmp = kmp,
        p = list(.p$plot + annotate("text", x = 5, y = 0.25, label = .label) + theme(legend.position = "", panel.grid = element_blank()))
      )
    }
  )) %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest() -> p62_survival_quantile40

p62_survival_quantile40 %>% 
  dplyr::mutate(cancer_types = reorder(cancer_types, hazard_ratio, sort)) %>% 
  ggplot(aes(x = hazard_ratio, y = cancer_types)) + 
  geom_point(aes(size = -log10(coxp))) +
  geom_vline(xintercept = 1)

gridExtra::arrangeGrob(grobs = p62_survival_quantile40$p, nrow = 2) %>% 
  ggsave(filename = "fig_01_p62_survival.pdf",
         plot = .,
         device = "pdf", 
         width = 8,
         height = 6,
         path = clinical_path)

# for harzard ratio
p62_clinical_nest %>% 
  dplyr::mutate(hazard_ratio = purrr::map2(
    .x = data,
    .y = cancer_types,
    .f = function(.x, .y){
      .x %>% 
        dplyr::mutate(time = time / 30) %>% 
        dplyr::mutate(time = ifelse(time > 60, 60, time)) %>% 
        dplyr::mutate(group = dplyr::case_when(
          rppa > quantile(rppa, 0.6) ~ "L",
          rppa < quantile(rppa, 0.4) ~ "H",
          TRUE ~ "M"
        )) -> .d
        .d %>% dplyr::filter(group != "M") -> .dd
    
      survival::coxph(survival::Surv(time, status) ~ rppa, data = .d) -> .cox
      summary(.cox) -> .z
      survival::survdiff(survival::Surv(time, status) ~ group, data = .d) -> .d_diff
      kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
      
      tibble::tibble(
        n = .z$n,
        hr = .z$conf.int[1],
        hr_l = .z$conf.int[3],
        hr_h = .z$conf.int[4],
        coxp = .z$waldtest[3],
        kmp = kmp)
      
    }
  )) %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest()  -> p62_hazard_ratio

p62_hazard_ratio %>% 
  dplyr::mutate(coxp = -log10(coxp)) %>% 
  dplyr::mutate(cancer_types = glue::glue("{cancer_types} ({n})")) %>% 
  dplyr::mutate(cancer_types = reorder(cancer_types, hr, sort)) -> plot_ready

plot_ready %>% 
  ggplot(aes(y = hr, x = cancer_types)) +
  geom_point(aes(size = coxp)) +
  geom_hline(aes(yintercept = 1)) +
  scale_size(name = "p-value") +
  coord_flip() +
  ggthemes::theme_gdocs() +
  labs(y = "Hazard Ratio", x = "Cancer Types") -> p
ggsave(filename = 'fig_01_hazard_ratio.pdf', plot = p, device = 'pdf', path = clinical_path)
knitr::kable(plot_ready %>% 
               dplyr::rename(`hazard ratio` = hr) %>% 
               dplyr::mutate(coxp = 10^-coxp) %>% 
               dplyr::arrange(-`hazard ratio`) %>% 
               dplyr::mutate(hr_l = signif(hr_l, 3), hr_h = signif(hr_h,3)) %>% 
               tidyr::unite(`95%CI`, hr_l, hr_h, sep = "-"))
# p62 stage -------------------------------------------------------------

clinical_stage <- readr::read_rds(file.path(tcga_path, "pancan34_clinical_stage.rds.gz")) %>% 
  dplyr::filter(n >= 40) %>% 
  dplyr::select(-n)

p62_sample_classification %>% 
  dplyr::inner_join(clinical_stage, by = "cancer_types") %>% 
  dplyr::mutate(pval = purrr::map2(
    .x = data,
    .y = stage,
    .f = function(.x, .y){
      .y %>%
        dplyr::rename(sample = barcode) %>%
        dplyr::inner_join(.x, by = "sample") %>% 
        dplyr::mutate(stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")))-> .d
      
      if (pval < 0.05) {

        CPCOLS <- c("#191970", "#CD0000")
        .d %>%
          ggplot(aes(x = stage, y = rppa, color = stage)) +
          geom_boxplot(outlier.colour = NA, width = 0.5) +
          # stat_boxplot(geom = 'errorbar', width = 0.2) +
          scale_x_discrete(limit = c("Stage I", "Stage IV"), labels = c("Stage I & II", "Stage III & IV")) +
          scale_color_manual(values = CPCOLS) +
          annotate("text", label = glue::glue("p = {signif(pval, 3)}"), x = 1.5, y = 2) +
          theme_bw() +
          labs(x = "Stage", y = "RPPA score") +
          theme(
            panel.grid = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal"
          )

      }
    }
  )) %>% 
  dplyr::filter(purrr::map_lgl(pval, Negate(is.null))) %>% 
  dplyr::select(cancer_types, pval) %>% 
  purrrlyr::by_row(..f = function(.d) {
    ggsave(filename = glue::glue("fig_02_stage_diff_{.d$cancer_types}.pdf"), plot = .d$pval[[1]], path = clinical_path, width = 3, height = 4)
  })

p62_sample_classification %>% 
  dplyr::inner_join(clinical_stage, by = "cancer_types") %>% 
  dplyr::mutate(pval = purrr::pmap(
    .l = list(
    .x = data,
    .y = stage,
    .z = cancer_types
    ),
    .f = function(.x, .y, .z){
      .y %>%
        dplyr::rename(sample = barcode) %>%
        dplyr::inner_join(.x, by = "sample") %>% 
        dplyr::mutate(stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))) -> .d
      
      comp_list <- list(c("Stage I", "Stage II"), c("Stage II", "Stage III"), c("Stage III", "Stage IV"))
      aov(rppa ~ stage, .d) %>% broom::glance() %>% .$p.value -> anova_p
      .d %>% 
        ggpubr::ggboxplot(x = "stage", y = "rppa",  color = "stage", pallete = "jco"  ) +
        ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test") + 
        ggpubr::stat_compare_means(method = "anova") +
        labs(x = .z, y = "") +
        theme(
          legend.position = 'none',
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        ggthemes::scale_color_gdocs() -> p
      
      tibble::tibble(anova_p = signif(anova_p, 3), p = list(p))
    }
  )) %>% 
  tidyr::unnest(pval) %>% 
  dplyr::select(-data, -stage) -> plot_ready
  
plot_ready %>% dplyr::arrange(anova_p) %>% dplyr::select(-anova_p) %>% tibble::deframe() -> p
gridExtra::arrangeGrob(grobs = p, nrow = 3,
                       left = grid::textGrob("p62 rppa score", rot = 90, vjust = 1)) -> grob

ggsave(filename = 'fig_02_all_stage.pdf', 
       plot = grob, device = 'pdf', 
       path = clinical_path, width = 15, height = 10)

# p62 subtype -------------------------------------------------------------

clinical_subtype <- 
  readr::read_rds(path = file.path(tcga_path,"pancan34_clinical_subtype.rds.gz")) %>% 
  dplyr::select(-n)

p62_sample_classification %>% 
  dplyr::inner_join(clinical_subtype, by = "cancer_types") %>% 
  dplyr::mutate(pval = purrr::pmap(
    .l = list(
      .x = data,
      .y = subtype,
      .z = cancer_types),
    .f = function(.x, .y, .z){
      .y %>% 
        dplyr::rename(sample = barcode) %>%
        dplyr::inner_join(.x, by = "sample") -> .d

      aov(rppa ~ subtype, data = .d) %>% broom::glance() %>% .$p.value -> pval
      
      if (pval < 0.05) {
        .label <- glue::glue("{.z}
                             ANOVA p = {signif(pval, 3)}")
        .d %>%
          ggplot(aes(x = reorder(subtype, rppa, median), y = rppa, 
                     color = reorder(subtype, rppa, median))) +
          geom_boxplot(outlier.size = NA) +
          ggthemes::scale_color_gdocs(name = "Subtype") +
          theme_bw() +
          labs(x = "Subtype", y = "RPPA score") +
          theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
          ) +
          annotate("text", label = .label, x = 1.5, y = 1.5)
      }
    }
  )) %>% 
  dplyr::filter(purrr::map_lgl(pval, Negate(is.null))) %>% 
  dplyr::pull(pval) %>% 
  gridExtra::arrangeGrob(grobs = ., nrow = 3) %>% 
  ggsave(filename = "fig_03_subtype_all.pdf",
         plot = .,
         device = "pdf", 
         width = 10,
         height = 8,
         path = clinical_path)

