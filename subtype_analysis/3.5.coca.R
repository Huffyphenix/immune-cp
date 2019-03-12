###########################################################
# combined cluster analysis for scaled exp, methylation clusters, and CNV raw data cluster results
###########################################################

library(SNFtool)
library(ComplexHeatmap)
library(magrittr)

# data path ---------------------------------------------------------------

result_path <- file.path("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis")
data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

# load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))
# load(file = file.path(data_result_path, ".rda_genelist_data_survival_cancer.info.rda"))
# 
# source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")


# load cluster results ----------------------------------------------------

expr_scaled_clsuer_4 <- readr::read_tsv(file.path(result_path,"expr_scaled/Get_best_clutser_20","4_clusters_group_info.tsv"))
methy_scaled_clsuer_6 <- readr::read_tsv(file.path(result_path,"methy_scaled/Get_best_clutser_20","6_clusters_group_info.tsv"))
cnv_raw_clsuer_6 <- readr::read_tsv(file.path(result_path,"CNV/Get_best_clutser_20","6_clusters_group_info.tsv"))

# recode the cluster result -----------------------------------------
expr_scaled_clsuer_4 %>%
  dplyr::mutate(encode = purrr::map(group,.f=function(.x){
    .d <- tibble::tibble()
    if(.x == 1) {.d <- tibble::tibble(a = 1,b = 0, c = 0, d= 0)}
    if(.x == 2) {.d <- tibble::tibble(a = 0,b = 1, c = 0, d = 0)}
    if(.x == 3) {.d <- tibble::tibble(a = 0,b = 0, c = 1, d = 0)}
    if(.x == 4) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 1)}
    .d
  })) %>%
  dplyr::select(-group) %>%
  tidyr::unnest() -> expr_scaled_clsuer_4.recode

expr_scaled_clsuer_4.recode %>%
  tidyr::gather(-sample,key="type",value="value") %>%
  dplyr::mutate(type=paste("exp",type,sep="_")) %>%
  tidyr::spread(key="type",value="value") -> expr_scaled_clsuer_4.recode.g

methy_scaled_clsuer_6 %>%
  dplyr::mutate(encode = purrr::map(group,.f=function(.x){
    .d <- tibble::tibble()
    if(.x == 1) {.d <- tibble::tibble(a = 1,b = 0, c = 0, d = 0, e = 0, f =0)}
    if(.x == 2) {.d <- tibble::tibble(a = 0,b = 1, c = 0, d = 0, e = 0, f =0)}
    if(.x == 3) {.d <- tibble::tibble(a = 0,b = 0, c = 1, d = 0, e = 0, f =0)}
    if(.x == 4) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 1, e = 0, f =0)}
    if(.x == 5) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 0, e = 1, f =0)}
    if(.x == 6) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 0, e = 0, f =1)}
    .d
  })) %>%
  dplyr::select(-group) %>%
  tidyr::unnest() -> methy_scaled_clsuer_6.recode

methy_scaled_clsuer_6.recode %>%
  tidyr::gather(-sample,key="type",value="value") %>%
  dplyr::mutate(type=paste("meth",type,sep="_"))  %>%
  tidyr::spread(key="type",value="value") -> methy_scaled_clsuer_6.recode.g

cnv_raw_clsuer_6 %>%
  dplyr::mutate(encode = purrr::map(group,.f=function(.x){
    .d <- tibble::tibble()
    if(.x == 1) {.d <- tibble::tibble(a = 1,b = 0, c = 0, d = 0, e = 0, f =0)}
    if(.x == 2) {.d <- tibble::tibble(a = 0,b = 1, c = 0, d = 0, e = 0, f =0)}
    if(.x == 3) {.d <- tibble::tibble(a = 0,b = 0, c = 1, d = 0, e = 0, f =0)}
    if(.x == 4) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 1, e = 0, f =0)}
    if(.x == 5) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 0, e = 1, f =0)}
    if(.x == 6) {.d <- tibble::tibble(a = 0,b = 0, c = 0, d = 0, e = 0, f =1)}
    .d
  })) %>%
  dplyr::select(-group) %>%
  tidyr::unnest() -> cnv_raw_clsuer_6.recode

cnv_raw_clsuer_6.recode %>%
  tidyr::gather(-sample,key="type",value="value") %>%
  dplyr::mutate(type=paste("cnv",type,sep="_"))  %>%
  tidyr::spread(key="type",value="value") -> cnv_raw_clsuer_6.recode.g

expr_scaled_clsuer_4.recode.g %>%
  dplyr::inner_join(methy_scaled_clsuer_6.recode.g,by="sample") %>%
  dplyr::inner_join(cnv_raw_clsuer_6.recode.g,by="sample") -> combine_data

t(combine_data[,-1]) %>%  as.matrix() -> combine_data.m
colnames(combine_data.m) <- combine_data$sample
# 
# t(methy_scaled_clsuer_6.recode) %>% as.data.frame() %>% .[-1,] -> methy_scaled_clsuer_6.recode.m
# colnames(methy_scaled_clsuer_6.recode.m) <- methy_scaled_clsuer_6.recode$sample
# 
# t(cnv_raw_clsuer_6.recode) %>% as.data.frame() %>% .[-1,] -> cnv_raw_clsuer_6.recode.m
# colnames(cnv_raw_clsuer_6.recode.m) <- cnv_raw_clsuer_6.recode$sample
# 
# combine_data =list(GeneExp=expr_scaled_clsuer_4.recode.m,methy=methy_scaled_clsuer_6.recode.m,cnv=cnv_raw_clsuer_6.recode.m)
# results <- ExecuteSNF(combine_data,clusterNum=20,K=10)
# results %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_expr-cnv-methy_recode.combined_20.rds.gz"), compress = 'gz')
# 
# W <- results$distanceMatrix

library(ConsensusClusterPlus)
ConsensusClusterPlus(combine_data.m, maxK=20, reps=500,pItem=0.8,pFeature=1, title="recode.combined.ConsensusClusterplus", clusterAlg="hc",distance="pearson",seed=1262118388.71279, plot = "pdf") -> cc_res

cc_res %>% readr::write_rds(path = file.path(data_result_path, ".recode.combined.ConsensusClusterplus.rds.gz"), compress = "gz")