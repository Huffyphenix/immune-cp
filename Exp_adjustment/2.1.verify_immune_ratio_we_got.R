
# verify the immune ration using published data ---------------------------
library(magrittr)

# data path -----
xcell_path <- "/project/huff/huff/data/TCGA/immune_infiltration/xCell"
purity_path <- "/project/huff/huff/data/TCGA_tumor_purity_from_ncomms9971"
immune_path <- "/project/huff/huff/immune_checkpoint"
immunity_path <- "/project/huff/huff/immune_checkpoint/data/immunity"

# load data -----
# immune ratio we calculated
xCell_TCGA_RSEM.immune_stroma.ratio <- readr::read_rds(file.path(immune_path,"genelist_data","xCell_TCGA_RSEM.immune_stroma.ratio.rds.gz")) %>%
  dplyr::select(-barcode) %>%
  dplyr::rename("barcode" = "Sample ID")

# TIMER data
TIMER_immunity <- readr::read_tsv(file.path(immunity_path,"immuneEstimation.txt")) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic)

# purity data
TCGA_purity <- readr::read_tsv(file.path(purity_path,"ncomms9971-s2.txt")) %>%
  dplyr::select(`Cancer type`,`Sample ID`,`CPE`) %>%
  dplyr::mutate(barcode=substr(`Sample ID`,1,15))

# data combination ----

TIMER_immunity %>%
  dplyr::inner_join(xCell_TCGA_RSEM.immune_stroma.ratio,by="barcode") -> immunity_combine_xCell.immnue.ratio


# distribution ------------------------------------------------------------
xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::filter(!is.na(Immune_ratio)) %>%
  dplyr::group_by(`Cancer type`) %>%
  dplyr::mutate(mid = quantile(Immune_ratio,0.5)) %>%
  dplyr::select(`Cancer type`,mid) %>%
  unique() %>%
  dplyr::arrange(mid) -> cancer_rank
xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::filter(!is.na(Immune_ratio)) %>%
  ggplot(aes(x=`Cancer type`,y=Immune_ratio)) +
  geom_violin(aes(fill = `Cancer type`),alpha = 0.5) +
  scale_x_discrete(limit = cancer_rank$`Cancer type`) +
  coord_flip() +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10),
    axis.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12,angle = 90),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20),
    axis.text = element_text(colour = "black")
  )
  
ggsave(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","ImmuneRatio_distribution_in_cancers.pdf"),device = "pdf",height = 4,width = 4)
ggsave(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","ImmuneRatio_distribution_in_cancers.png"),device = "png",height = 4,width = 4)
# calculate correlation ---------------------------------------------------

fn_correlation <- function(.x){
  broom::tidy(cor.test(.x$Immune_ratio,.x$TIL)) %>%
    dplyr::mutate(x="Immune_ratio",y="Total TIL") -> .tmp1
  broom::tidy(cor.test(.x$Immune_ratio,.x$B_cell)) %>%
    dplyr::mutate(x="Immune_ratio",y="B_cell") -> .tmp2
  broom::tidy(cor.test(.x$Immune_ratio,.x$CD4_Tcell)) %>%
    dplyr::mutate(x="Immune_ratio",y="CD4_Tcell") -> .tmp3
  broom::tidy(cor.test(.x$Immune_ratio,.x$CD8_Tcell)) %>%
    dplyr::mutate(x="Immune_ratio",y="CD8_Tcell") -> .tmp4
  broom::tidy(cor.test(.x$Immune_ratio,.x$Neutrophil)) %>%
    dplyr::mutate(x="Immune_ratio",y="Neutrophil") -> .tmp5
  broom::tidy(cor.test(.x$Immune_ratio,.x$Macrophage)) %>%
    dplyr::mutate(x="Immune_ratio",y="Macrophage") -> .tmp6
  broom::tidy(cor.test(.x$Immune_ratio,.x$Dendritic)) %>%
    dplyr::mutate(x="Immune_ratio",y="Dendritic") -> .tmp7
  rbind(.tmp1,.tmp2,.tmp3,.tmp4,.tmp5,.tmp6,.tmp7)
}

immunity_combine_xCell.immnue.ratio %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::group_by(`Cancer type`) %>%
  dplyr::mutate(cor = purrr::map(data,fn_correlation)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() -> correlation_TIMER.our_results

immunity_combine_xCell.immnue.ratio %>%
  dplyr::mutate(`Cancer type` = "Total") %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::group_by(`Cancer type`) %>%
  dplyr::mutate(cor = purrr::map(data,fn_correlation)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() -> Total.correlation_TIMER.our_results
# draw plot -----
library(ggplot2)
CPCOLS <- c("red", "white", "#1C86EE")
Total.correlation_TIMER.our_results %>%
  rbind(correlation_TIMER.our_results) %>%
  ggplot(aes(y=`Cancer type`,x=y)) +
  geom_tile(aes(fill=estimate),colour = "white") +
  scale_fill_gradient2(
    name = "Correlation", #"Methylation diff (T - N)",
    low = CPCOLS[3],
    high = CPCOLS[1],
    mid = CPCOLS[2],
    breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1)
  ) +
  guides(fill=guide_colorbar(title.position="left",
                             title.theme = element_text(angle = 90))) +
  scale_y_discrete(limit = c(sort(correlation_TIMER.our_results$`Cancer type` %>% unique()),"Total")) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, angle = 90, size = 10),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12,angle = 90),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.title = element_text(size = 20),
    axis.text = element_text(colour = "black")
  ) -> p
p
ggsave(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","TIMER","heatmap_correlation.TIMER.ourImmuneRatio.pdf"),device = "pdf",height = 6,width = 4)
ggsave(file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","TIMER","heatmap_correlation.TIMER.ourImmuneRatio.png"),device = "png",height = 6,width = 4)


# survival and immune ratio -----------------------------------------------
clinical <- readr::read_rds(file.path("/project/huff/huff/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz"))
clinical %>%
  tidyr::unnest() %>%
  dplyr::select(bcr_patient_barcode,PFS,PFS.time) %>%
  tidyr::drop_na() %>%
  dplyr::rename("barcode" = "bcr_patient_barcode")-> time_status

fn_survival <-function(.cancer,.x){
  .up=75
  .low=25
  .x %>%
    dplyr::filter(!is.na(Immune_ratio)) %>%
    dplyr::mutate(group = ifelse(Immune_ratio>quantile(Immune_ratio,.up/100),"High",NA)) %>%
    dplyr::mutate(group = ifelse(Immune_ratio<quantile(Immune_ratio,.low/100),"Low",group))-> .data
  .d_diff <- try(survival::survdiff(survival::Surv(PFS.time, PFS) ~ group, data = .data),silent = TRUE)
  if ('try-error' %in% class(.d_diff)) {return(data.frame(KMP=1,Coxp=1))}else{
    # print(.d_diff)
    kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.data$group))) - 1)
    # print(kmp)
    coxp <-  broom::tidy(survival::coxph(survival::Surv(PFS.time, PFS) ~ Immune_ratio, data = .data, na.action = na.exclude))
    # print(coxp)
    
    .fit <- survival::survfit(survival::Surv(PFS.time, PFS) ~ group, data = .data, na.action = na.exclude) 
    low_n=.fit$n[2]
    high_n=.fit$n[1]
    legend <- data.frame(group=c("Low","High"),n=c(low_n,high_n))
    legend %>%
      dplyr::mutate(
        label = purrr::map2(
          .x = group,
          .y = n,
          .f = function(.x,.y){
            latex2exp::TeX(glue::glue("IR^{<<.x>>}, n = <<.y>>", .open = "<<", .close = ">>"))
          }
        )
      ) %>%
      dplyr::arrange(group)-> legend
    survminer::ggsurvplot(.fit,pval=F, #pval.method = T,
                          surv.median.line = "hv",
                          title = paste(paste(.cancer, sep = "-"), "p =", signif(kmp, 2)),
                          xlab = "Time (days)",
                          ylab = 'Survival',
                          legend.title = "Expression group:",
                          legend= c(0.8,0.8),
                          ggtheme = theme(
                            panel.border = element_blank(), panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.5), 
                            panel.background = element_rect(fill = "white"),
                            legend.key = element_blank(),
                            legend.background = element_blank(),
                            legend.text = element_text(size = 12, colour = "black"),
                            axis.text = element_text(size = 12, colour = "black"),
                            legend.title = element_blank(),
                            axis.title = element_text(size = 12,color = "black")
                          ),
                          font.main = c(16),
                          font.x = c(14),
                          font.y = c(14),
                          font.tickslab = c(12)
    ) +
      scale_color_manual(
        values = c( "#EE6363","#1C86EE"),
        labels = legend$label
      )
    # dev.off()
    fig_name <- paste(.cancer,"PFS",signif(kmp, 2),"pdf", sep = ".")
    fig_name_png <- paste(.cancer,"PFS",signif(kmp, 2),"png", sep = ".")
    ggsave(filename = fig_name, device = "pdf", path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","survival","75-25"), width = 4, height = 3)
    ggsave(filename = fig_name_png, device = "png", path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_immune_ratio","survival","75-25"), width = 4, height = 3)
    # print(.fit)
    # fn_sur_draw(.gene,kmp,.fit)
    data.frame(KMP=kmp,Coxp=coxp) -> .out
    return(.out)
  }
}

xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::mutate(sur = purrr::map2(`Cancer type`,data,fn_survival)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> PFS_50_ImmuneRatio_survival
PFS_50_ImmuneRatio_survival %>%
  dplyr::filter(KMP <= 0.05)
xCell_TCGA_RSEM.immune_stroma.ratio %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  tidyr::nest(-`Cancer type`) %>%
  dplyr::mutate(sur = purrr::map2(`Cancer type`,data,fn_survival)) %>%
  dplyr::select(-data) %>%
  dplyr::ungroup() %>%
  tidyr::unnest() -> PFS_75.25_ImmuneRatio_survival
PFS_75.25_ImmuneRatio_survival %>%
  dplyr::filter(KMP <= 0.05)
