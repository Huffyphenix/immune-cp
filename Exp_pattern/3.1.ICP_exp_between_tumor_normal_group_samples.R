########################################
# seperate the samples by ICP expression compared with normal expression
########################################

library(magrittr)
library(dplyr)

# processed path
basic_path <- "/home/huff/project"
# tcga_path = file.path(basic_path,"immune_checkpoint/data/TCGA_data")
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
expr_path <- file.path(basic_path,"immune_checkpoint/result_20171025/expr_rds")


# load data ---------------------------------------------------------------
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
expr <- readr::read_rds(file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))


# filter data -------------------------------------------------------------
# filter_gene_list <- function(.x, gene_list) {
#   gene_list %>%
#     dplyr::select(symbol) %>%
#     dplyr::left_join(.x, by = "symbol") %>%
#     dplyr::select(-entrez_id) %>%
#     tidyr::gather(-symbol,key="barcode",value="exp")
# }
# 
# expr %>%
#   dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
#   dplyr::select(-expr) -> gene_list_expr.nest

expr %>%
  dplyr::mutate(gat_exp = purrr::map(filter_expr, .f=function(.x){
    .x %>%
      tidyr::gather(-symbol,-entrez_id,key="barcode",value="exp") 
  })) %>%
  dplyr::select(-filter_expr) %>%
  tidyr::unnest() %>%
  dplyr::mutate(group = ifelse(substr(barcode,14,14)==1,"Normal","Tumor")) %>%
  dplyr::mutate(Participant = substr(barcode,1,12)) -> gene_list_expr.nest.TNgrouped


# gene mean expression in normal samples each cancer ----------------------

gene_list_expr.nest.TNgrouped %>%
  dplyr::filter(group == "Normal") %>%
  dplyr::group_by(cancer_types,symbol) %>%
  dplyr::mutate(mean_exp = mean(exp)) %>%
  dplyr::select(-barcode,-exp,-Participant) %>%
  dplyr::ungroup() %>%
  unique() ->  gene_list_normal_mean_expr

# gene peak expression in normal samples each cancer (density distribution)----------------------
fn_density_peak <- function(.name,.x){
  print(.name)
  # get peak and secondary peak
  set.seed(42)
  .x <- .x$exp
  if(length(as.vector(.x))<10 ){
    tibble::tibble(peak.x = NA, peak.y = NA)
  }else if(length(unique(as.vector(.x)))<3){
    tibble::tibble(peak.x = mean(.x), peak.y = NA)
  }else{
    d1 <- density(.x)
    d1.ym <- which.max(d1$y)
    d1.ym.x <- d1$x[d1.ym]
    d1.ym.y <- d1$y[d1.ym]
    
    tibble::tibble(peak.x = d1.ym.x, peak.y = d1.ym.y)
  }
}
gene_list_expr.nest.TNgrouped %>%
  dplyr::filter(group == "Normal") %>%
  tidyr::nest(-cancer_types,-symbol) %>%
  tidyr::unite(can_sym,c("cancer_types","symbol")) %>%
  dplyr::mutate(peak = purrr::map2(can_sym,data,fn_density_peak)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  tidyr::separate(can_sym,c("cancer_types","symbol"),"_") -> gene_list_normal_peak_expr


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


#compare expression in tumor samples with [mean exp] in normal samples --------------

gene_list_expr.nest.TNgrouped %>%
  dplyr::filter(group == "Tumor") %>%
  dplyr::filter(cancer_types %in% normal_10samples_in_cancers$cancer_types) %>%
  dplyr::inner_join(gene_list_normal_mean_expr,by=c("cancer_types","symbol")) %>%
  dplyr::rename("normal_mean_exp" = "mean_exp") %>%
  dplyr::mutate(`log2(T/N)` = log2((exp+1)/(normal_mean_exp+1))) %>%
  dplyr::select(cancer_types,symbol,barcode,exp,normal_mean_exp,`log2(T/N)`) -> gene_list_expr.T_N.by_mean

# compare expression in tumor samples with [density peak exp] in normal samples  --------------

gene_list_expr.nest.TNgrouped %>%
  dplyr::filter(group == "Tumor") %>%
  dplyr::filter(cancer_types %in% normal_10samples_in_cancers$cancer_types) %>%
  dplyr::inner_join(gene_list_normal_peak_expr,by=c("cancer_types","symbol")) %>%
  dplyr::rename("normal_peak_exp" = "peak.x") %>%
  dplyr::mutate(`log2(T/N)` = log2(exp/normal_peak_exp)) %>%
  dplyr::select(cancer_types,symbol,barcode,exp,normal_peak_exp,`log2(T/N)`) -> gene_list_expr.T_N.by_peak


# only paired T-N samples -------------------------------------------------
gene_list_expr.nest.TNgrouped %>%
  dplyr::filter(group == "Normal") %>%
  dplyr::select(Participant) %>%
  unique() -> Participant_have_normal_samples

gene_list_expr.nest.TNgrouped %>% 
  dplyr::filter(Participant %in% Participant_have_normal_samples$Participant) %>%
  dplyr::filter(substr(barcode,15,15)==1) %>%
  dplyr::select(-barcode) %>%
  tidyr::spread(key="group",value="exp") %>%
  dplyr::mutate(`log2(T/N)` = log2((Tumor+1)/(Normal+1))) -> gene_list_expr.T_N.only_paired
  


# group tumor samples by gene counts -----------------------------------------------------
ICP_expr_pattern <-
  readr::read_tsv(file.path(immune_path,"result_20171025","ICP_exp_patthern-byMeanUQ","pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(symbol,Exp_site)
# ICP_expr_pattern<- readr::read_tsv(file.path(immune_path,"result_20171025","ICP_exp_patthern-byratio","pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv"))

fn_group_sample_by_ICP_exp_in_TN_1 <- function(.name,.data){
  print(.name)
  .data %>%
    dplyr::filter(!is.na(`log2(T/N)`)) %>%
    dplyr::filter(!is.na(`Exp_site`)) -> .tmp
  if(dim(.tmp)[1]>0){
    .tmp %>%
      dplyr::mutate(sample_status = purrr::map2(`log2(T/N)`,`Exp_site`,fn_group_sample_by_ICP_exp_in_TN_2)) %>%
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

# paired tumor-normal samples grouping ----
gene_list_expr.T_N.only_paired %>%
  dplyr::inner_join(ICP_expr_pattern,by="symbol") %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Immune cell dominate","Immune and tumor cell almost"),"Useful","Useless")) %>%
  dplyr::mutate(Exp_site=Exp_site.1) %>%
  dplyr::select(-Exp_site.1) %>%
  tidyr::nest(-cancer_types,-Participant) %>%
  tidyr::unite(c_p,c("cancer_types","Participant")) %>%
  dplyr::mutate(status = purrr::map2(c_p,data,fn_group_sample_by_ICP_exp_in_TN_1)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.only_paired

tumor_class_by_T_N.only_paired %>%
  tidyr::separate(c_p,c("cancer_types","Participant"),"_") %>%
  dplyr::group_by(Participant,sample_status) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,Participant,sample_status,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="sample_status",value="n") %>%
  dplyr::mutate(class = ifelse(Immunity_cold>Immunity_hot,"Immunity_cold","Immunity_hot")) %>%
  dplyr::mutate(class = ifelse(Immunity_cold==Immunity_hot,"not clear",class)) -> tumor_class_by_T_N.only_paired.class
 
tumor_class_by_T_N.only_paired.class %>% 
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio","tumor_class_by_T_N.only_paired.by_geneCounts","tumor_class_by_T_N.only_paired"))

# scored tumor samples by mean FC product of all genes-----------------------------------------------------
# gene's FC in fantom * FC in tumor-normal
fn_score_samples_byFC <- function(.x,.y){
  FC_T_N_tcga <- 2^.x
  FC_fantom <- 2^.y
  FC_T_N_tcga*FC_fantom
}

# paired tumor-normal samples grouping ----
gene_list_expr.T_N.only_paired %>%
  dplyr::inner_join(ICP_expr_pattern) %>%
  dplyr::mutate(score = purrr::map2(`log2(T/N)`,`log2FC(I/T)`,fn_score_samples_byFC)) %>%
  dplyr::select(cancer_types,symbol,Participant,score) %>%
  tidyr::unnest() -> gene_list_expr.T_N.only_paired.gene_score

gene_list_expr.T_N.only_paired.gene_score %>%
  dplyr::group_by(Participant) %>%
  dplyr::mutate(score = ifelse(is.na(score),0,score))%>%
  dplyr::mutate(score_mean = mean(score)) %>%
  dplyr::select(cancer_types,Participant,score_mean) %>%
  unique() %>%
  dplyr::ungroup()-> gene_list_expr.T_N.only_paired.sample_score

gene_list_expr.T_N.only_paired.sample_score %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio","tumor_class_by_T_N.only_paired.by_FCProduct","tumor_class_by_T_N.only_paired"))

# all tumor samples grouping ----
gene_list_expr.T_N.by_mean %>%
  dplyr::inner_join(ICP_expr_pattern) %>%
  dplyr::mutate(score = purrr::map2(`log2(T/N)`,`log2FC(I/T)`,fn_score_samples_byFC)) %>%
  dplyr::select(cancer_types,symbol,barcode,score) %>%
  tidyr::unnest() -> gene_list_expr.T_N.by_mean.gene_score

gene_list_expr.T_N.by_mean.gene_score %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(score = ifelse(is.na(score),0,score))%>%
  dplyr::mutate(score_mean = mean(score)) %>%
  dplyr::select(cancer_types,barcode,score_mean) %>%
  unique() %>%
  dplyr::ungroup()-> gene_list_expr.T_N.by_mean.sample_score

gene_list_expr.T_N.by_mean.sample_score %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio","tumor_class_by_T_N.all-by-mean.by_FCProduct","tumor_class_by_T_N.all_tumor.by_mean"))

# scored tumor samples by log2 FC product of all genes-----------------------------------------------------
# sum(log2(FCg1.1)+log2(FCg1.2)+log2(FCg2.1)+log2(FCg2.2)+...=log2(FCg1.1*FCg1.2*FCg2.1*FCg2.2....))
# get each gene's score log2 FC in fantom + log2 FC in tumor-normal, and get each sample's score by sum(gene's score)
fn_score_samples_bylog2FC_sum <- function(.x,.y){
  FC_T_N_tcga <- .x
  FC_fantom <- .y
  FC_T_N_tcga+FC_fantom
}

# paired tumor-normal samples grouping ----
gene_list_expr.T_N.only_paired %>%
  dplyr::inner_join(ICP_expr_pattern) %>%
  dplyr::mutate(score = purrr::map2(`log2(T/N)`,`log2FC(I/T)`,fn_score_samples_bylog2FC_sum)) %>%
  dplyr::select(cancer_types,symbol,Participant,score) %>%
  tidyr::unnest() -> gene_list_expr.T_N.only_paired.gene_score.by_log2FCproduct

gene_list_expr.T_N.only_paired.gene_score.by_log2FCproduct %>%
  dplyr::filter(!is.na(score)) %>%
  dplyr::group_by(Participant) %>%
  dplyr::mutate(score_sum = sum(score)) %>%
  dplyr::select(cancer_types,Participant,score_sum) %>%
  unique() %>%
  dplyr::ungroup()-> gene_list_expr.T_N.only_paired.sample_score.by_log2FCproduct

gene_list_expr.T_N.only_paired.sample_score.by_log2FCproduct %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio","tumor_class_by_T_N.only_paired.by_logFCProduct_sum","tumor_class_by_T_N.only_paired.by_logFCProduct_sum"))

# all tumor samples grouping ----
gene_list_expr.T_N.by_mean %>%
  dplyr::inner_join(ICP_expr_pattern) %>%
  dplyr::mutate(score = purrr::map2(`log2(T/N)`,`log2FC(I/T)`,fn_score_samples_bylog2FC_sum)) %>%
  dplyr::select(cancer_types,symbol,barcode,score) %>%
  tidyr::unnest() -> gene_list_expr.T_N.gene_score.by_mean.by_log2FCproduct_sum

gene_list_expr.T_N.gene_score.by_mean.by_log2FCproduct_sum %>%
  dplyr::filter(!is.na(score)) %>%
  dplyr::group_by(barcode) %>%
  dplyr::mutate(score = ifelse(is.na(score),0,score))%>%
  dplyr::mutate(score_sum = sum(score)) %>%
  dplyr::select(cancer_types,barcode,score_sum) %>%
  unique() %>%
  dplyr::ungroup()-> gene_list_expr.T_N.sample_score.by_mean.by_log2FCproduct_sum

gene_list_expr.T_N.sample_score.by_mean.by_log2FCproduct_sum %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio","tumor_class_by_T_N.all-by-mean.by_log2FCProduct_sum","tumor_class_by_T_N.all_tumor.by_mean.by_log2FCProduct_sum"))
#### following code was not used anymore
#### not run
# using peak exp value compare with every tumor samples to group tumor samples ----
gene_list_expr.T_N.by_peak %>%
  dplyr::inner_join(ICP_expr_pattern) %>%
  tidyr::nest(-cancer_types,-barcode) %>%
  tidyr::unite(c_p,c("cancer_types","barcode")) %>%
  dplyr::mutate(status = purrr::map2(c_p,data,fn_group_sample_by_ICP_exp_in_TN_1)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.by_peak

tumor_class_by_T_N.by_peak %>%
  tidyr::separate(c_p,c("cancer_types","barcode"),"_") %>%
  dplyr::group_by(barcode,sample_status) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,barcode,sample_status,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="sample_status",value="n") %>%
  dplyr::mutate(class = ifelse(Immunity_cold>Immunity_hot,"Immunity_cold","Immunity_hot")) %>%
  dplyr::mutate(class = ifelse(Immunity_cold==Immunity_hot,"not clear",class)) -> tumor_class_by_T_N.by_peak.class

tumor_class_by_T_N.by_peak.class %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_peak"))
# using peak exp value compare with every tumor samples to group tumor samples ----

gene_list_expr.T_N.by_mean %>%
  dplyr::inner_join(ICP_expr_pattern) %>%
  tidyr::nest(-cancer_types,-barcode) %>%
  tidyr::unite(c_p,c("cancer_types","barcode")) %>%
  dplyr::mutate(status = purrr::map2(c_p,data,fn_group_sample_by_ICP_exp_in_TN_1)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> tumor_class_by_T_N.by_mean

tumor_class_by_T_N.by_mean %>%
  tidyr::separate(c_p,c("cancer_types","barcode"),"_") %>%
  dplyr::group_by(barcode,sample_status) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types,barcode,sample_status,n) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="sample_status",value="n") %>%
  dplyr::mutate(class = ifelse(Immunity_cold>Immunity_hot,"Immunity_cold","Immunity_hot")) %>%
  dplyr::mutate(class = ifelse(Immunity_cold==Immunity_hot,"not clear",class)) -> tumor_class_by_T_N.by_mean.class
tumor_class_by_T_N.by_mean.class %>%
  readr::write_tsv(file.path(immune_path,"result_20171025/ICP_exp_patthern","tumor_class_by_T_N.by_mean"))


# save image --------------------------------------------------------------

save.image(file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio","tumor_class_estimation_of_TCGAsamples.rda"))
