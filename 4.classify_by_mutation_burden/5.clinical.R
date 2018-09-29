########################################
### clinical difference between mutation groups 
########################################

library(magrittr)
library(ggplot2)

# data path ---------------------------------------------------------------

immune_path <- "/project/huff/huff/immune_checkpoint"
burden_path <- "/project/huff/huff/data/TCGA"
tcga_path <- file.path(immune_path,"/data/TCGA_data")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint") 
result_path <- "/project/huff/huff/immune_checkpoint/result/mutation_load_class"

# load data ---------------------------------------------------------------

mutation_burden_class <- readr::read_rds(file.path(burden_path,"classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest()

clinical <- readr::read_rds(file.path(tcga_path,"pancan34_clinical.rds.gz"))


# merge data --------------------------------------------------------------

fn_merge <- function(cli,cancer){
  print(cancer)
  cli %>%
    dplyr::select(barcode,os_days,os_status) %>%
    dplyr::inner_join(mutation_burden_class,by="barcode") %>%
    dplyr::select(-cancer_types) %>%
    dplyr::mutate(os_status = ifelse(os_status == "Dead",1,0))
}
mutation_burden_class$cancer_types %>% unique() ->cancers_in_mutaion_burden_class
clinical %>%
  dplyr::filter(cancer_types %in% cancers_in_mutaion_burden_class) %>%
  dplyr::mutate(cli_snv_merge = purrr::map2(.x=clinical,cancer_types,fn_merge)) %>%
  dplyr::select(-clinical) -> clinical_snv_merged


# survival analysis -------------------------------------------------------
## cancers level analysis 
library(survival)
fn_survival <- function(.x){
  broom::tidy(survival::coxph(survival::Surv(os_days, os_status) ~ mutation_status, data = .x, na.action = na.exclude)) %>%
    dplyr::mutate(logrnakP = 1 - pchisq(survival::survdiff(survival::Surv(os_days, os_status) ~ mutation_status, data = .x, na.action = na.exclude)$chisq, df = length(levels(as.factor(.x$mutation_status))) - 1) )
}

cancers_with_enough_high_mutation <- c("LUAD","SKCM","BLCA","COAD","LUSC","STAD","HNSC","UCEC","BRCA","LIHC","READ","ACC")
clinical_snv_merged %>%
  dplyr::filter(cancer_types %in% cancers_with_enough_high_mutation) %>%
  dplyr::mutate(surv = purrr::map(cli_snv_merge,fn_survival)) %>%
  dplyr::select(-cli_snv_merge) %>%
  tidyr::unnest() -> clinical_snv_surv_pvalue

## pancan survival analysi
clinical_snv_merged %>%
  tidyr::unnest() -> tmp
  
cli_snv_cox <- survival::coxph(survival::Surv(os_days, os_status) ~ mutation_status, data = tmp, na.action = na.exclude)
logrnakP <- 1 - pchisq(survival::survdiff(survival::Surv(os_days, os_status) ~ mutation_status, data = tmp, na.action = na.exclude)$chisq, df = length(levels(as.factor(tmp$mutation_status))) - 1)
