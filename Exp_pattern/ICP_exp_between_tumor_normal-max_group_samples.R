########################################
# seperate the samples by ICP expression compared with normal expression
########################################

library(magrittr)
library(dplyr)

# processed path
tcga_path = c("/project/huff/huff/immune_checkpoint/data/TCGA_data")
immune_path <- "/project/huff/huff/immune_checkpoint"
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")


# load data ---------------------------------------------------------------
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz"))


# filter data -------------------------------------------------------------
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol") %>%
    dplyr::select(-entrez_id) %>%
    tidyr::gather(-symbol,key="barcode",value="exp")
}

expr %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) -> gene_list_expr.nest

gene_list_expr.nest %>%
  tidyr::unnest() %>%
  dplyr::mutate(group = ifelse(substr(barcode,14,14)==1,"Normal","Tumor")) %>%
  dplyr::mutate(Participant = substr(barcode,1,12)) -> gene_list_expr.nest.TNgrouped


# gene max expression in normal samples each cancer ----------------------

gene_list_expr.nest.TNgrouped %>%
  dplyr::filter(group == "Normal") %>%
  dplyr::group_by(cancer_types,symbol) %>%
  dplyr::mutate(max_exp = max(exp)) %>%
  dplyr::select(-barcode,-exp,-Participant) %>%
  dplyr::ungroup() %>%
  unique() ->  gene_list_normal_max_expr

# get cancers with > 10 normal samples ------------------------------------

gene_list_expr.nest.TNgrouped %>%
  dplyr::filter(group == "Normal") %>%
  dplyr::select(-symbol,-exp) %>%
  unique() %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,n) %>%
  dplyr::filter(n>=10) %>%
  dplyr::ungroup() %>%
  unique() -> normal_10samples_in_cancers

#compare expression in tumor samples with [max exp] in normal samples --------------

gene_list_expr.nest.TNgrouped %>%
  dplyr::filter(group == "Tumor") %>%
  dplyr::filter(cancer_types %in% normal_10samples_in_cancers$cancer_types) %>%
  dplyr::inner_join(gene_list_normal_max_expr,by=c("cancer_types","symbol")) %>%
  dplyr::rename("normal_max_exp" = "max_exp") %>%
  dplyr::mutate(`log2(T/N)` = log2(exp/normal_max_exp)) %>%
  dplyr::select(cancer_types,symbol,barcode,exp,normal_max_exp,`log2(T/N)`) -> gene_list_expr.T_N.by_max

# group tumor samples -----------------------------------------------------

ICP_expr_pattern <- readr::read_tsv(file.path(immune_path,"result_20171025","ICP_exp_patthern","manual_edit_2_ICP_exp_pattern_in_immune_tumor_cell.tsv"))

fn_group_sample_by_ICP_exp_in_TN_1 <- function(.name,.data){
  print(.name)
  .data %>%
    dplyr::filter(!is.na(`log2(T/N)`)) %>%
    dplyr::filter(!is.na(`Exp site`)) -> .tmp
  if(dim(.tmp)[1]>0){
    .tmp %>%
      dplyr::mutate(sample_status = purrr::map2(`log2(T/N)`,`Exp site`,fn_group_sample_by_ICP_exp_in_TN_2)) %>%
      tidyr::unnest() %>%
      dplyr::select(symbol,sample_status,description)
  }else{
    tibble::tibble(symbol=NA,sample_status=NA,description=NA)
  }
}
fn_group_sample_by_ICP_exp_in_TN_2 <- function(.x,.y){
  if(.x < 0){
    if(.y=="Mainly_Immune"){
      tibble::tibble(sample_status="Immunity_cold",description="Immune_withdraw/absent")
    }else if(.y=="Both"){
      tibble::tibble(sample_status="Immunity_cold",description="Immune_withdraw/absent & Tumor_downregulate_immuneInhibity")
    }else{
      tibble::tibble(sample_status="Immunity_cold",description="Tumor_downregulate_immuneInhibity")
    }
  }else if(.x > 0){
    if(.y=="Mainly_Immune"){
      tibble::tibble(sample_status="Immunity_hot",description="Immune_enter")
    }else if(.y=="Both"){
      tibble::tibble(sample_status="Immunity_hot",description="Immune_enter & Tumor_upregulate_immuneInhibity")
    }else{
      tibble::tibble(sample_status="Immunity_hot",description="Tumor_upregulate_immuneInhibity")
    }
  }else{
    tibble::tibble(sample_status="not clear",description="difficult to distinct")
  }
}

# using max exp value compare with every tumor samples to group tumor samples ----

gene_list_expr.T_N.by_max %>%
  # head() %>%
  dplyr::inner_join(ICP_expr_pattern) %>%
  tidyr::nest(-cancer_types,-barcode) %>%
  tidyr::unite(c_p,c("cancer_types","barcode")) %>%
  dplyr::mutate(status = purrr::map2(c_p,data,fn_group_sample_by_ICP_exp_in_TN_1)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.by_max

tumor_class_by_T_N.by_max %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_max-gene_details"))

tumor_class_by_T_N.by_max %>%
  tidyr::separate(c_p,c("cancer_types","barcode"),"_") %>%
  dplyr::group_by(barcode,sample_status) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,barcode,sample_status,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="sample_status",value="n") %>%
  dplyr::mutate(class = ifelse(Immunity_cold>Immunity_hot,"Immunity_cold","Immunity_hot")) %>%
  dplyr::mutate(class = ifelse(Immunity_cold==Immunity_hot,"not clear",class)) -> tumor_class_by_T_N.by_max.class
tumor_class_by_T_N.by_max.class %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_max"))


# score -------------------------------------------------------------------


tumor_class_by_T_N.by_max.class 
