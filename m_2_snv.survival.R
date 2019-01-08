#snv with survival
library(magrittr)
library(ggplot2)
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")

# load cnv and gene list
snv <- readr::read_rds(file.path(tcga_path, "pancan33_snv.rds.gz"))
gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

snv %>%
  dplyr::mutate(filter_snv = purrr::map(snv, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-snv) -> gene_list_snv

#load clinical data
clinical <- readr::read_rds(path = file.path(tcga_path,"pancan34_clinical.rds.gz")) 

# merge clinical and snv
snv_clinical <- 
  gene_list_snv %>%
  # dplyr::filter(cancer_types %in% cancer_pairs$cancer_types) %>% 
  dplyr::inner_join(clinical, by = "cancer_types")

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

fun_expr_survival_merge <- function(filter_snv, clinical){
  # merge clinical and expr data
  filter_snv %>% 
    #dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = snv, -symbol) %>% 
    dplyr::select(symbol, barcode, snv)  %>% 
    dplyr::mutate(snv=as.character(snv)) %>%
    #dplyr::mutate(snv=ifelse(is.na(snv),"non mut",snv)) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(snv=plyr::revalue(snv,replace = c("0"="non mut","NA"="non mut"))) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(snv_class=ifelse(snv!="non mut","mut",snv)) %>%
    dplyr::select(symbol,barcode,snv_class) %>%
    dplyr::inner_join(clinical, by = "barcode") %>% 
    dplyr::select(symbol, barcode, snv_class, gender, race,time = os_days, status = os_status) %>% 
    dplyr::filter(!is.na(time), time > 0, !is.na(status)) %>% 
    dplyr::mutate(status = plyr::revalue(status, replace = c("Alive" = 0, "Dead" = 1))) %>%
    dplyr::mutate(status = as.numeric(status)) %>% 
    #dplyr::mutate(expr = log2(expr + 1)) %>% 
    #tidyr::drop_na(expr) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(snv_class=as.character(snv_class)) %>%
    dplyr::mutate(group = snv_class) %>% 
    dplyr::ungroup() -> expr_clinical_ready
} 

fun_draw_survival <- function(symbol, p.value, cancer_types, expr_clinical_ready){
  # symbol="CD96"
  # p.value=0.05
  # cancer_types="SKCM"
  gene <- symbol
  p_val <- signif(-log10(p.value), digits = 3)
  fig_name <- paste(cancer_types, gene, p_val, "pdf", sep = ".")
  # print(fig_name)
  .d <- 
    expr_clinical_ready %>% 
    dplyr::filter(symbol == gene)
  
  .d_diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = .d)
  
  kmp <- 1 - pchisq(.d_diff$chisq, df = length(levels(as.factor(.d$group))) - 1)
  
  if(kmp > 0.05) {return(NA)} else{
    fit_x <- survival::survfit(survival::Surv(time, status) ~ group, data = .d , na.action = na.exclude)
    #pdf(file = file.path(snv_path,"snv_survival",fig_name),width = 480,height = 480)
    survminer::ggsurvplot(fit_x, data = .d, pval=T, pval.method = T,
                          tables.height = 0.2,
                          tables.theme = theme_cleantable(),
                          title = paste(paste(cancer_types, gene, sep = "-"), "cox p =", signif(p.value, 2)),
                          xlab = "Survival in days",
                          ylab = 'Probability of survival'
                          )
    #dev.off()
    ggsave(filename = fig_name, device = "pdf", path = file.path(snv_path, "snv_survival"), width = 6, height = 6)
  }
}
fun_clinical_test <- function(expr_clinical_ready, cancer_types){
  if(nrow(expr_clinical_ready) < 1){return(tibble::tibble())}
  #if(expr_clinical_ready %>% dplyr::filter(group=="mut") %>% nrow() <nrow(expr_clinical_ready)*0.02){return(tibble::tibble())}
  
  print(cancer_types)
  expr_clinical_ready %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::summarise(sm_mut=sum(snv_class=="mut"),sm_nm=sum(snv_class=="non mut")) %>%
    dplyr::mutate(mut_per=sm_mut/(sm_mut+sm_nm)) %>%
    dplyr::filter(sm_mut>=3) %>%
    dplyr::pull(symbol) ->ready_gene
  if(ready_gene %>% length() ==0){return(tibble::tibble())}else{
    expr_clinical_ready %>%
      dplyr::filter(symbol %in% ready_gene) %>%
      dplyr::group_by(symbol) %>%
        dplyr::do(
          broom::tidy(
            tryCatch(
              survival::coxph(survival::Surv(time, status) ~ group, data = ., na.action = na.exclude),
              error = function(e){1}#,
              #warning = function(e){2}
            )
          )
       ) %>%
      dplyr::ungroup() %>% 
      #dplyr::filter(p.value < 0.05) %>% 
      dplyr::select(symbol, estimate, p.value) -> d
  d %>% 
    dplyr::select(symbol, p.value) %>% 
    purrr::pwalk(fun_draw_survival, cancer_types = cancer_types, expr_clinical_ready = expr_clinical_ready) 
  
  return(d)
  }
}

snv_clinical %>%
  dplyr::mutate(merged_clean = purrr::map2(filter_snv, clinical, fun_expr_survival_merge)) %>%
  dplyr::select(-filter_snv, -clinical) %>%
  dplyr::mutate(diff_pval = purrr::map2(merged_clean, cancer_types, fun_clinical_test)) %>%
  dplyr::select(-merged_clean) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  #dplyr::select(-PARTITION_ID) %>% 
  tidyr::unnest(diff_pval) -> expr_clinical_sig_pval
expr_clinical_sig_pval %>%
  dplyr::filter(p.value<=0.05) %>%
  readr::write_csv(file.path(snv_path,"snv_survival","c_3_snv_survival_genelist_sig_pval.csv"))

save.image(file = file.path(snv_path, "snv_survival",".rda_03_d_survival_gene_expr.rda"))
load(file = file.path(survival_path, ".rda_03_d_survival_gene_expr.rda"))



