library(magrittr)
library(ggplot2)

out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
survival_path <- file.path(out_path, "c_3_survival")
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")

# clinical <- readr::read_rds(path = file.path(tcga_path,"pancan34_clinical.rds.gz")) 
clinical <- readr::read_rds(file.path("/project/huff/huff/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type")

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

cancer_pairs <- readr::read_tsv(file = file.path(expr_path, "tsv_02_pancan_samples_pairs.tsv"))


# merge clinical and expr
expr_clinical <- 
  gene_list_expr %>%
  # dplyr::filter(cancer_types %in% cancer_pairs$cancer_types) %>% 
  dplyr::inner_join(clinical, by = "cancer_types") %>%
  dplyr::rename("clinical"="data")

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
} # get tumor and normal info

fun_expr_survival_merge <- function(filter_expr, clinical){
  # merge clinical and expr data
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type == "01") %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(symbol, barcode, expr)  %>% 
    dplyr::rename("bcr_patient_barcode" = "barcode") %>%
    dplyr::inner_join(clinical, by = "bcr_patient_barcode") %>% 
    dplyr::select(symbol, bcr_patient_barcode, expr, time = PFS.time, status = PFS) %>%  #gender, race,
    dplyr::filter(!is.na(time), time > 0, !is.na(status)) %>% 
    dplyr::mutate(time=time/365) %>% 
    # dplyr::mutate(status = plyr::revalue(status, replace = c("Alive" = 0, "Dead" = 1))) %>%
    dplyr::mutate(status = as.numeric(status)) %>% 
    dplyr::mutate(expr = log2(expr + 1)) %>% 
    tidyr::drop_na(expr) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(group = as.factor(ifelse(expr <= median(expr),"Low", "High"))) %>% 
    dplyr::ungroup() -> expr_clinical_ready
} 

fun_draw_survival <- function(symbol, p.value, cancer_types, expr_clinical_ready){
  gene <- symbol
  p_val <- signif(p.value, digits = 3)
  fig_name <- paste(cancer_types, gene, p_val, sep = ".")
  # print(fig_name)
  .d <- 
    expr_clinical_ready %>% 
    dplyr::filter(symbol == gene)
  
  .d_diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = .d)
  
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
  
  if(kmp > 0.1) {return(NA)} else{
    fit_x <- survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude)
    survminer::ggsurvplot(fit_x, data = .d, pval=T, pval.method = T,
               title = paste(paste(cancer_types, gene, sep = "-"), "Coxph =", signif(p.value, 2)),
               surv.median.line = "hv",
               xlab = "Survival in years",
               ylab = 'Probability of survival',
               ggtheme = theme(
                 panel.border = element_blank(), 
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 axis.line = element_line(colour = "black", size = 0.5), 
                 panel.background = element_rect(fill = "white"),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 8, colour = "black"),
                 axis.text = element_text(size = 12, colour = "black"),
                 legend.title = element_blank(),
                 axis.title = element_text(size = 12,color = "black")
               ))[[1]] +
      scale_color_manual(
        values = c("red4","dodgerblue4"),
        labels = paste(c("High","Low"),",n =",fit_x$n,sep="")
      )
    ggsave(filename = paste(fig_name,"pdf",sep="."), device = "pdf", path = file.path(survival_path, "survival"), width = 4, height = 3)
    ggsave(filename = paste(fig_name,"png",sep="."), device = "png", path = file.path(survival_path, "survival"), width = 4, height = 3)
  }
}
fun_clinical_test <- function(expr_clinical_ready, cancer_types){
  if(nrow(expr_clinical_ready) < 1){return(tibble::tibble())}
  print(cancer_types)
  # cox p
  expr_clinical_ready %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::do(
     broom::tidy(
       tryCatch(
         survival::coxph(survival::Surv(time, status) ~ expr, data = ., na.action = na.exclude),
         error = function(e){1},
         warning = function(e){2}
         )
       )
      ) %>%
    dplyr::ungroup() %>% 
    # dplyr::filter(p.value < 0.05) %>% 
    dplyr::select(symbol, estimate, p.value) %>% 
    dplyr::mutate(status = ifelse(estimate > 0, "High_risk", "Low_risk"))-> d
  
  # kmp
  expr_clinical_ready %>% 
    dplyr::group_by(symbol) %>%
    dplyr::do(
      broom::tidy(
        tryCatch(
          1 - pchisq(
            survival::survdiff(survival::Surv(time, status) ~ group, data = ., na.action = na.exclude)$chisq,
            df = 1),
          error = function(e){1},
          warning = function(e){2}
        )
      )
    ) %>%
    dplyr::ungroup() %>% 
    dplyr::rename("kmp"="x") %>%
    # dplyr::filter(p.value < 0.05) %>% 
    dplyr::select(symbol, kmp) -> d.kmp
  
  d %>% 
    dplyr::inner_join(d.kmp,by="symbol") %>%
    dplyr::filter(p.value <= 0.05 | kmp <= 0.05) %>%
    dplyr::select(symbol, p.value) %>% 
    purrr::pwalk(fun_draw_survival, cancer_types = cancer_types, expr_clinical_ready = expr_clinical_ready) 
  
  return(d %>% 
           dplyr::inner_join(d.kmp,by="symbol"))
}

# expr_clinical %>%
#   dplyr::filter(cancer_types == "KIRC") %>% 
#   dplyr::mutate(merged_clean = purrr::map2(filter_expr, clinical, fun_expr_survival_merge)) %>%
#   dplyr::select(-filter_expr, -clinical) %>%
#   dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fun_clinical_test)) %>%
#   dplyr::select(-merged_clean) %>%
#   tidyr::unnest(diff_pval) -> expr_clinical_sig_pval

#>>>>>>>>>>>>>>>packages：survminer and ggpubr need R update.
# cl <- parallel::detectCores()
# cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
# expr_clinical %>%
#   multidplyr::partition(cluster = cluster) %>%
#   multidplyr::cluster_library("magrittr") %>%
#   multidplyr::cluster_library("ggplot2") %>%
#   multidplyr::cluster_library("survminer") %>%
#   multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
#   multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type) %>% 
#   multidplyr::cluster_assign_value("fun_expr_survival_merge", fun_expr_survival_merge) %>% 
#   multidplyr::cluster_assign_value("fun_clinical_test", fun_clinical_test) %>% 
#   multidplyr::cluster_assign_value("fun_draw_survival", fun_draw_survival) %>% 
#   multidplyr::cluster_assign_value("survival_path", survival_path) %>% 
  
expr_clinical %>%
  dplyr::mutate(merged_clean = purrr::map2(filter_expr, clinical, fun_expr_survival_merge)) %>%
  dplyr::select(-filter_expr, -clinical) %>%
  # head(1) %>%
  dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fun_clinical_test)) %>%
  dplyr::select(-merged_clean) %>% 
  # dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  #dplyr::select(-PARTITION_ID) %>% 
  tidyr::unnest(diff_pval) -> expr_clinical_pval
expr_clinical_pval %>%
  readr::write_csv(file.path(survival_path,"PFS_c_3_survival_genelist_pval.csv"))

expr_clinical_pval %>%
  dplyr::filter(p.value <= 0.05 | kmp <= 0.05) -> expr_clinical_sig_pval
expr_clinical_sig_pval %>%
  readr::write_csv(file.path(survival_path,"PFS_c_3_survival_genelist_sig_pval.csv"))
# on.exit(parallel::stopCluster(cluster))
#---------------------------------------------------------------------------------------------

fun_rank_cancer <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(., na.rm = T))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
} #get cancer rank
fun_rank_gene <- function(pattern){
  pattern %>% 
    dplyr::rowwise() %>%
    dplyr::do(
      symbol = .$symbol,
      rank =  unlist(.[-1], use.names = F) %>% sum(na.rm = T)
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest() %>%
    dplyr::arrange(rank)
} # get gene rank

expr_clinical_sig_pval %>% 
  dplyr::select(cancer_types, symbol) %>% 
  dplyr::mutate(n = 1) %>% 
  tidyr::spread(key = cancer_types, value = n) -> pattern

cancer_rank <- pattern %>% fun_rank_cancer()
gene_rank <- 
  pattern %>% 
  fun_rank_gene() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::filter( rank >= 5) %>% 
  dplyr::mutate(color = plyr::revalue(status, replace = c('a' = "#e41a1c", "l" = "#377eb8", "i" = "#4daf4a", "p" = "#984ea3"))) %>% 
  dplyr::arrange(color, rank)

expr_clinical_sig_pval %>% 
  ggplot(aes(x = cancer_types, y = symbol, color = status)) +
  geom_point(aes(size = -log10(p.value))) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(name = "P-value") +
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
    axis.text.y = element_text(color = gene_rank$color),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) +
  ggthemes::scale_color_gdocs(name = "Surivival Worse")-> p
ggsave(
  filename = "fig_03_d_survival_sig_genes_serminar.pdf",
  plot = p,
  device = "pdf",
  width = 14,
  height = 9,
  path = survival_path
)



save.image(file = file.path(survival_path, ".rda_03_d_survival_gene_expr.rda"))
load(file = file.path(survival_path, ".rda_03_d_survival_gene_expr.rda"))
