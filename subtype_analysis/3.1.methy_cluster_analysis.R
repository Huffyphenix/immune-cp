###### do subtype analysis for each data type
#load library
library(methods)
library(magrittr)
library(CancerSubtypes)

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))

result_path <- "/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis"
# expr ----
expr_cc <- readr::read_rds(file.path(data_result_path, ".rds_PanCan28_expr_cc_10.rds.gz"))
group_cc=expr_cc$group
distanceMatrix_cc=expr_cc$distanceMatrix
p_value=survAnalysis(mainTitle="cnv_sur_result",expr_time,expr_status,group_cc,
                     distanceMatrix_cc,similarity=TRUE)

# CNV ----
cnv_cc <- readr::read_rds(file.path(data_result_path, ".rds_PanCan28_cnv_cc_10.rds.gz"))
group_cc=cnv_cc$group
group_cc %>% table()
distanceMatrix_cc=cnv_cc$distanceMatrix
p_value=survAnalysis(mainTitle="cnv_sur_result",cnv_time,cnv_status,group_cc,
                     distanceMatrix_cc,similarity=TRUE)

# combine -----
combine_cc <- readr::read_rds(file.path(data_result_path, ".rds_PanCan28_combine_cc_10.rds.gz"))
group_cc=combine_cc$group;group_cc %>% table()
distanceMatrix_cc=combine_cc$distanceMatrix
p_value=survAnalysis(mainTitle="methy_sur_result",time.combine,status.combine,group_cc,
                     distanceMatrix_cc,similarity=TRUE)

# Methylation ----
# only methylation cluster showed association with survival
methy_cc <- readr::read_rds(file.path(data_result_path, ".rds_PanCan28_methy_cc_C10.rds.gz"))
group_cc=methy_cc$group
group_cc %>% table()
# distanceMatrix_cc=methy_cc$distanceMatrix
# p_value=survAnalysis(mainTitle="methy_sur_result",methy_time,methy_status,group_cc,
                     # distanceMatrix_cc,similarity=TRUE)



# survival analysis --------
data.frame(sample = names(group_cc),group=group_cc,time = methy_time,status = methy_status) %>%
  dplyr::as.tbl() %>%
  dplyr::filter(group %in% c("1","4","8","9","10")) -> methy_group_survival_data
data.frame(sample = names(group_cc),group=group_cc,time = expr_time,status = expr_status) %>%
  dplyr::as.tbl() %>%
  dplyr::filter(group %in% c("1","2","3","4","8"))-> methy_group_survival_data
data.frame(sample = names(group_cc),group=group_cc,time = cnv_time,status = cnv_status) %>%
  dplyr::as.tbl() -> methy_group_survival_data
library(survival)
library(ggplot2)
methy_fit <- survfit(survival::Surv(time, status) ~ group, data = methy_group_survival_data, na.action = na.exclude)
methy_diff <- survdiff(survival::Surv(time, status) ~ group, data = methy_group_survival_data, na.action = na.exclude)
kmp <- 1 - pchisq(methy_diff$chisq, df = length(levels(as.factor(methy_group_survival_data$group))) - 1)
legend <- data.frame(group=paste("C",unique(methy_group_survival_data$group),sep=""),n=methy_fit$n)
legend %>%
  dplyr::mutate(
    label = purrr::map2(
      .x = group,
      .y = n,
      .f = function(.x,.y){
        latex2exp::TeX(glue::glue("<<.x>>, n = <<.y>>", .open = "<<", .close = ">>"))
      }
    )
  ) -> legend
survminer::ggsurvplot(methy_fit,pval=F, #pval.method = T,
                      surv.median.line = "hv",
                      title = paste("Methylation cluster, p =", signif(kmp, 2)), # change it when doing diff data
                      xlab = "Survival in days",
                      ylab = 'Probability of survival',
                      # legend.title = "Methyla group:",
                      legend= c(0.8,0.6),
                      # ggtheme = theme_survminer(),
                      ggtheme = theme(
                        panel.border = element_blank(), panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", 
                                                                                     size = 0.5), 
                        panel.background = element_rect(fill = "white"),
                        legend.key = element_blank(),
                        legend.background = element_blank(),
                        legend.text = element_text(size = 12, colour = "black"),
                        axis.text = element_text(size = 12, colour = "black"),
                        legend.title = element_blank(),
                        axis.title = element_text(size = 12,color = "black")
                      )
) +
  scale_color_manual(
    values = c( "#EE6363","#1C86EE","#8B8B00",c("#7B68EE", "#EE00EE", "#008B00", "#000000", "#00C5CD", "#FFA54F", "#0000CD")),
    labels = legend$label
  ) 
ggsave(filename = "methy_cluster_3_OS.pdf", device = "pdf", path = file.path(result_path,"methy"), width = 4, height = 3)
ggsave(filename = "methy_cluster_3_OS.png", device = "png", path = file.path(result_path,"methy"), width = 4, height = 3)

ggsave(filename = "methy_cluster_10_OS_main_cluster.pdf", device = "pdf", path = file.path(result_path,"methy"), width = 4, height = 3)
ggsave(filename = "methy_cluster_10_OS_main_cluster.png", device = "png", path = file.path(result_path,"methy"), width = 4, height = 3)

ggsave(filename = "expr_cluster_10_OS_main_cluster.pdf", device = "pdf", path = file.path(result_path,"expr"), width = 5, height = 4)
ggsave(filename = "expr_cluster_10_OS_main_cluster.png", device = "png", path = file.path(result_path,"expr"), width = 5, height = 4)

ggsave(filename = "CNV_cluster_10_OS.pdf", device = "pdf", path = file.path(result_path,"CNV"), width = 5, height = 4)
ggsave(filename = "CNV_cluster_10_OS.png", device = "png", path = file.path(result_path,"CNV"), width = 5, height = 4)
# muttaion burden difference ----
methy_group_survival_data %>%
  dplyr::rename("barcode" = "sample") %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") -> methy_cluster_mutation

b <- runif(nrow(methy_cluster_mutation), -0.3, 0.3)
comp_list <- list(c("1","2"),c("1","3"),c("1","4"),c("1","8"),c("2","3"),c("2","4"),c("2","8"),c("3","4"),c("3","8"),c("4","8")) # expr compare list
comp_list <- list(c("1","4"),c("1","8"),c("1","9"),c("1","10"),c("4","8"),c("4","9"),c("4","10"),c("8","9"),c("8","10"),c("9","10")) # expr compare list
methy_cluster_mutation %>%
  dplyr::mutate(sm_count = log2(sm_count)) %>%
  ggpubr::ggboxplot(x = "group", y = "sm_count",
                    color = "group" #add = "jitter",#, palette = "npg"
  ) +
  # geom_point(aes(x=as.numeric(group)+b,y=sm_count,color=group),alpha = 0.5) +
  scale_x_discrete(breaks = c(1:10),
                   labels = paste("C",c(1:10),sep="")
                   # expand = c(0.2,0.2,0.2)
                   ) +
  # facet_wrap(~ cancer_types, strip.position = "bottom", scales = "free") +
  scale_color_manual(
    values =  c( "#EE6363","#1C86EE","#8B8B00",c("#7B68EE", "#EE00EE", "#008B00", "#000000", "#00C5CD", "#FFA54F", "#0000CD"))
  )+
  # ylim(4,12) +
  ylab("log2(Mutation burden)") +
  xlab("Methylation clusters") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white",colour = "white"),
        strip.text = element_text(size = 12)) +
  # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")   -> p;p
ggsave(filename = "mutationBurdenDiff_between_methy_cluster_10.pdf", device = "pdf", path = file.path(result_path,"methy"), width = 4, height = 3)
ggsave(filename = "mutationBurdenDiff_between_methy_cluster_10.png", device = "png", path = file.path(result_path,"methy"), width = 4, height = 3)

ggsave(filename = "mutationBurdenDiff_between_methy_cluster_10_cancer_types.pdf", device = "pdf", path = file.path(result_path,"methy"), width = 6, height = 5)
ggsave(filename = "mutationBurdenDiff_between_methy_cluster_10_cancer_types.png", device = "png", path = file.path(result_path,"methy"), width = 6, height = 5)

ggsave(filename = "mutationBurdenDiff_between_expr_cluster_10.pdf", device = "pdf", path = file.path(result_path,"expr"), width = 6, height = 5)
ggsave(filename = "mutationBurdenDiff_between_expr_cluster_10.png", device = "png", path = file.path(result_path,"expr"), width = 6, height = 5)
