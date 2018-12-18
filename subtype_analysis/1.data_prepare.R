###### prepare data for subtype analysis

library(magrittr)
library(ggplot2)

# data path config
immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
burden_path <- "/project/huff/huff/data/TCGA"

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

mutation_burden_class <- readr::read_rds(file.path(burden_path,"classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest()

# 4. clinical data --------------------------------------------------------
clinical <- readr::read_rds(file.path(tcga_path,"pancan34_clinical.rds.gz"))

fn_merge <- function(cli,cancer){
  print(cancer)
  cli %>%
    dplyr::select(barcode,os_days,os_status) %>%
    dplyr::inner_join(mutation_burden_class,by="barcode") %>%
    dplyr::select(-cancer_types) %>%
    dplyr::mutate(os_status = ifelse(os_status == "Dead",1,0))
}
clinical %>%
  dplyr::filter(cancer_types %in% cancers_in_mutaion_burden_class) %>%
  dplyr::mutate(cli_snv_merge = purrr::map2(.x=clinical,cancer_types,fn_merge)) %>%
  dplyr::select(-clinical) %>%
  tidyr::unnest() %>%
  dplyr::select(barcode,os_days,os_status) %>%
  tidyr::drop_na() -> time_status

time_status %>%
  readr::write_tsv(file.path(data_result_path,"time_status.tsv"))

# 1. prepare expression data ------------------------------------------------------

expr_data <- readr::read_rds(file.path(tcga_path,"pancan33_expr.rds.gz"))

gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

expr_data %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) %>%
  dplyr::filter(cancer_types %in% mutation_burden_class$cancer_types) -> gene_list_expr

gene_list_expr %>%
  tidyr::unnest() %>%
  tidyr::gather(-cancer_types,-symbol,-entrez_id,key="barcode",value="expr") %>%
  dplyr::mutate(T_N = ifelse(substr(barcode,14,14) == 1,"Normal","Tumor")) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::select(-cancer_types) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(!is.na(expr)) %>%
  dplyr::mutate(mutation_status = ifelse(T_N=="Normal","Normal",mutation_status)) %>%
  dplyr::mutate(sm_count = ifelse(T_N=="Normal",0,sm_count)) %>%
  dplyr::inner_join(time_status,by="barcode") %>%   # samples of expr and clinical should be the same
  dplyr::filter(! is.na(expr)) %>% 
  dplyr::arrange(barcode) -> gene_list_expr_with_mutation_load # arrange is important for the sample corresponding between survival time and expr data. cause spread will arrange the barcode auto

gene_list_expr_with_mutation_load %>% 
  dplyr::select(barcode,os_days) %>%
  unique() %>% .$os_days -> expr_time

gene_list_expr_with_mutation_load %>% 
  dplyr::select(barcode,os_status) %>%
  unique() %>% .$os_status -> expr_status

gene_list_expr_with_mutation_load %>%
  dplyr::filter(T_N == "Tumor") %>%
  dplyr::select(symbol,barcode,expr) %>%
  dplyr::group_by(symbol,barcode) %>%
  dplyr::mutate(expr = mean(expr)) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key = barcode, value = expr) -> PanCan26_gene_list_expr_spread
  
PanCan26_gene_list_expr_matrix <- as.matrix(PanCan26_gene_list_expr_spread[,-1])
rownames(PanCan26_gene_list_expr_matrix) <- PanCan26_gene_list_expr_spread$symbol

PanCan26_gene_list_expr_matrix %>%
  readr::write_rds(file.path(data_result_path,"PanCan26_gene_list_expr_matrix.rds.gz"),compress = "gz")


# 2. prepare CNV data (raw data)------------------------------------------------------
cnv <- readr::read_rds(file.path(tcga_path,"pancan34_cnv.rds.gz"))

mutation_burden_class$cancer_types %>% unique() ->cancers_in_mutaion_burden_class

cnv %>%
  dplyr::filter(cancer_types %in% cancers_in_mutaion_burden_class) %>%    # some cancers in TCGA merged serveral cancers into one
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=cnv) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  dplyr::filter(! is.na(cnv)) %>% 
  dplyr::arrange(barcode) -> cnv_merge_snv_data

cnv_merge_snv_data %>% 
  dplyr::select(barcode,os_days) %>%
  unique() %>% .$os_days -> cnv_time

cnv_merge_snv_data %>% 
  dplyr::select(barcode,os_status) %>%
  unique() %>% .$os_status -> cnv_status

cnv_merge_snv_data %>%
  dplyr::select(symbol,barcode,cnv) %>%
  tidyr::spread(key = barcode,value=cnv) -> PanCan26_gene_list_cnv_spread

PanCan26_gene_list_cnv_matrix <- as.matrix(PanCan26_gene_list_cnv_spread[,-1])
rownames(PanCan26_gene_list_cnv_matrix) <- PanCan26_gene_list_cnv_spread$symbol

PanCan26_gene_list_cnv_matrix %>%
  readr::write_rds(file.path(data_result_path,"PanCan26_gene_list_cnv_matrix.rds.gz"))


# 3. methylation data  ----------------------------------------------------

mthy <- readr::read_rds(file.path(tcga_path,"pancan33_meth.rds.gz"))

mthy %>%
  dplyr::filter(cancer_types %in% cancers_in_mutaion_burden_class) %>%    # some cancers in TCGA merged serveral cancers into one
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  dplyr::select(-gene) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=methy) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  dplyr::filter(! is.na(methy)) %>% 
  dplyr::arrange(barcode) -> genelist_methy_mutaion_class  

genelist_methy_mutaion_class %>% 
  dplyr::select(barcode,os_days) %>%
  unique() %>% .$os_days -> methy_time

genelist_methy_mutaion_class %>% 
  dplyr::select(barcode,os_status) %>%
  unique() %>% .$os_status -> methy_status

genelist_methy_mutaion_class %>%
  dplyr::select(symbol,barcode,methy) %>%
  dplyr::group_by(symbol,barcode) %>%
  dplyr::mutate(methy = mean(methy)) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key = barcode,value = methy) -> PanCan26_gene_list_methy_spread

PanCan26_gene_list_methy_matrix <- as.matrix(PanCan26_gene_list_methy_spread[,-1])
rownames(PanCan26_gene_list_methy_matrix) <- PanCan26_gene_list_methy_spread$symbol

PanCan26_gene_list_methy_matrix %>%
  readr::write_rds(file.path(data_result_path,"PanCan26_gene_list_methy_matrix.rds.gz"),compress = "gz")


# data combine ------------------------------------------------------------
genelist_methy_mutaion_class %>% .$barcode %>% unique() -> methy_samples
cnv_merge_snv_data %>% .$barcode %>% unique() -> cnv_samples
gene_list_expr_with_mutation_load %>% .$barcode %>% unique() -> expr_samples
intersect(methy_samples,cnv_samples) %>% intersect(expr_samples) -> samples.intersect

genelist_methy_mutaion_class %>% 
  dplyr::filter(barcode %in% samples.intersect) -> genelist_methy_mutaion_class.combine
cnv_merge_snv_data %>% 
  dplyr::filter(barcode %in% samples.intersect) -> cnv_merge_snv_data.combine
gene_list_expr_with_mutation_load %>% 
  dplyr::filter(barcode %in% samples.intersect) -> gene_list_expr_with_mutation_load.combine
# save work space ---------------------------------------------------------

save.image(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))
