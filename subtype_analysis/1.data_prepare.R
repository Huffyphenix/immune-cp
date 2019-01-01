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
  tidyr::unnest() %>%
  dplyr::rename("cancer_types"="Cancer_Types")
mutation_burden_class$cancer_types %>% unique() ->cancers_in_mutaion_burden_class
cancers_except <- c("COADREAD","GBMLGG","KIPAN","STES","CESC")
# 4. clinical data --------------------------------------------------------
clinical <- readr::read_rds(file.path("/project/huff/huff/data/survival","TCGA_pancan_cancer_cell_survival_time.rds.gz"))

# fn_merge <- function(cli,cancer){
#   print(cancer)
#   cli %>%
#     dplyr::select(barcode,os_days,os_status) %>%
#     dplyr::inner_join(mutation_burden_class,by="barcode") %>%
#     dplyr::select(-cancer_types) %>%
#     dplyr::mutate(os_status = ifelse(os_status == "Dead",1,0))
# }
clinical %>%
  dplyr::filter(! type %in% cancers_except) %>%
  # dplyr::mutate(cli_snv_merge = purrr::map2(.x=clinical,type,fn_merge)) %>%
  # dplyr::select(-clinical) %>%
  tidyr::unnest() %>%
  dplyr::select(bcr_patient_barcode,PFS,PFS.time) %>%
  tidyr::drop_na() %>%
  dplyr::rename("barcode" = "bcr_patient_barcode")-> time_status

# time_status %>%
#   readr::write_tsv(file.path(data_result_path,"time_status.tsv"))

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
  dplyr::filter(cancer_types %in% mutation_burden_class$cancer_types) -> gene_list_expr.nest


gene_list_expr.nest %>%
  dplyr::filter(!cancer_types %in% cancers_except) %>%    # some cancers in TCGA merged serveral cancers into one
  tidyr::unnest() %>%
  tidyr::gather(-cancer_types,-symbol,-entrez_id,key="barcode",value="expr") %>%
  dplyr::mutate(T_N = ifelse(substr(barcode,14,14) == 1,"Normal","Tumor")) %>%
  dplyr::filter(T_N == "Tumor") %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  # dplyr::select(-cancer_types) %>%
  # dplyr::left_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(!is.na(expr)) %>%
  # dplyr::mutate(mutation_status = ifelse(T_N=="Normal","Normal",mutation_status)) %>%
  # dplyr::mutate(sm_count = ifelse(T_N=="Normal",0,sm_count)) %>%
  dplyr::inner_join(time_status,by="barcode") %>%   # samples of expr and clinical should be the same
  dplyr::filter(! is.na(expr)) %>% 
  dplyr::arrange(barcode) -> gene_list_expr # arrange is important for the sample corresponding between survival time and expr data. cause spread will arrange the barcode auto

gene_list_expr %>%
  dplyr::select(cancer_types,barcode) %>%
  unique() -> gene_list_expr.cancer_info

gene_list_expr %>% 
  dplyr::filter(T_N == "Tumor") %>%
  dplyr::select(barcode,PFS.time) %>%
  unique() %>% .$PFS.time -> expr_time

gene_list_expr %>% 
  dplyr::filter(T_N == "Tumor") %>%
  dplyr::select(barcode,PFS) %>%
  unique() %>% .$PFS -> expr_status

gene_list_expr %>%
  dplyr::filter(T_N == "Tumor") %>%
  dplyr::select(symbol,barcode,expr) %>%
  dplyr::group_by(symbol,barcode) %>%
  dplyr::mutate(expr = mean(expr)) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key = barcode, value = expr) -> PanCan26_gene_list_expr_spread
  
PanCan26_gene_list_expr_matrix <- as.matrix(PanCan26_gene_list_expr_spread[,-1])
rownames(PanCan26_gene_list_expr_matrix) <- PanCan26_gene_list_expr_spread$symbol

# PanCan26_gene_list_expr_matrix %>%
#   readr::write_rds(file.path(data_result_path,"PanCan26_gene_list_expr_matrix.rds.gz"),compress = "gz")


# 2. prepare CNV data (raw data)------------------------------------------------------
cnv <- readr::read_rds(file.path(tcga_path,"pancan34_cnv.rds.gz"))



cnv %>%
  dplyr::filter(!cancer_types %in% cancers_except) %>%    # some cancers in TCGA merged serveral cancers into one
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=cnv) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  # dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  dplyr::filter(! is.na(cnv)) %>% 
  dplyr::arrange(barcode) -> cnv_merge_snv_data

cnv_merge_snv_data %>%
  dplyr::select(cancer_types,barcode) %>%
  unique() -> cnv_merge_snv_data.cancer_info

cnv_merge_snv_data %>% 
  dplyr::select(barcode,PFS.time) %>%
  unique() %>% .$PFS.time -> cnv_time

cnv_merge_snv_data %>% 
  dplyr::select(barcode,PFS) %>%
  unique() %>% .$PFS -> cnv_status

cnv_merge_snv_data %>%
  dplyr::select(symbol,barcode,cnv) %>%
  tidyr::spread(key = barcode,value=cnv) -> PanCan26_gene_list_cnv_spread

PanCan26_gene_list_cnv_matrix <- as.matrix(PanCan26_gene_list_cnv_spread[,-1])
rownames(PanCan26_gene_list_cnv_matrix) <- PanCan26_gene_list_cnv_spread$symbol

PanCan26_gene_list_cnv_matrix %>%
  readr::write_rds(file.path(data_result_path,"PanCan26_gene_list_cnv_matrix.rds.gz"))

# 2. prepare SNV data (raw data)------------------------------------------------------
snv <- readr::read_rds(file.path("/data/shiny-data/GSCALite/TCGA/snv","pancan33_snv_from_syn7824274_gather.rds.gz"))

mutation_burden_class$cancer_types %>% unique() ->cancers_in_mutaion_burden_class

snv %>%
  # dplyr::filter(! Cancer_Types %in% cancers_except) %>%    # some cancers in TCGA merged serveral cancers into one
  # tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  dplyr::rename("barcode" = "sample","snv"="mut_n") %>%
  dplyr::select(symbol,barcode,snv) %>%
  # dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  dplyr::filter(! is.na(snv)) %>% 
  dplyr::arrange(barcode) %>%
  unique() %>%
  dplyr::mutate(snv = ifelse(snv >0, 1,snv))-> snv_merge_snv_data

snv_merge_snv_data %>% 
  dplyr::select(barcode,PFS.time) %>%
  unique() %>% .$PFS.time -> snv_time

snv_merge_snv_data %>% 
  dplyr::select(barcode,PFS) %>%
  unique() %>% .$PFS -> snv_status

snv_merge_snv_data %>%
  dplyr::select(symbol,barcode,snv) %>%
  tidyr::spread(key = barcode,value=snv) -> PanCan26_gene_list_snv_spread

PanCan26_gene_list_snv_matrix <- as.matrix(PanCan26_gene_list_snv_spread[,-1])
rownames(PanCan26_gene_list_snv_matrix) <- PanCan26_gene_list_snv_spread$symbol

# PanCan26_gene_list_snv_matrix %>%
#   readr::write_rds(file.path(data_result_path,"PanCan26_gene_list_snv_matrix.rds.gz"))

# 3. methylation data  ----------------------------------------------------

mthy <- readr::read_rds(file.path(tcga_path,"pancan33_meth.rds.gz"))

mthy %>%
  dplyr::filter(!cancer_types %in% cancers_except) %>%    # some cancers in TCGA merged serveral cancers into one
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  dplyr::select(-gene) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=methy) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  # dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::inner_join(time_status,by="barcode") %>%
  dplyr::filter(! is.na(methy)) %>% 
  dplyr::arrange(barcode) -> genelist_methy_mutaion_class  

genelist_methy_mutaion_class %>%
  dplyr::select(cancer_types,barcode) %>%
  unique() -> genelist_methy_mutaion_class.cancer_info

genelist_methy_mutaion_class %>% 
  dplyr::select(barcode,PFS.time) %>%
  unique() %>% .$PFS.time -> methy_time

genelist_methy_mutaion_class %>% 
  dplyr::select(barcode,PFS) %>%
  unique() %>% .$PFS -> methy_status

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
gene_list_expr %>% .$barcode %>% unique() -> expr_samples
# snv_merge_snv_data %>% .$barcode %>% unique() -> snv_samples
intersect(methy_samples,cnv_samples) %>% intersect(expr_samples) -> samples.intersect

genelist_methy_mutaion_class %>% 
  dplyr::filter(barcode %in% samples.intersect) %>%
  dplyr::select(cancer_types,barcode) %>%
  unique() -> combined.cancer_info

genelist_methy_mutaion_class %>% 
  dplyr::filter(barcode %in% samples.intersect) %>%
  dplyr::select(symbol,barcode,methy) %>%
  dplyr::group_by(symbol,barcode) %>%
  dplyr::mutate(methy = mean(methy)) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key = barcode,value = methy)-> genelist_methy_mutaion_class.combine.spread

PanCan26_gene_list_methy_matrix.combine <- as.matrix(genelist_methy_mutaion_class.combine.spread[,-1])
rownames(PanCan26_gene_list_methy_matrix.combine) <- genelist_methy_mutaion_class.combine.spread$symbol

cnv_merge_snv_data %>% 
  dplyr::filter(barcode %in% samples.intersect) %>%
  dplyr::select(symbol,barcode,cnv) %>%
  dplyr::group_by(symbol,barcode) %>%
  dplyr::mutate(cnv = mean(cnv)) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key = barcode,value = cnv)-> cnv_merge_snv_data.combine
cnv_merge_snv_data.matrix.combine <- as.matrix(cnv_merge_snv_data.combine[,-1])
rownames(cnv_merge_snv_data.matrix.combine) <- cnv_merge_snv_data.combine$symbol

# snv_merge_snv_data %>% 
#   dplyr::filter(barcode %in% samples.intersect) %>%
#   dplyr::select(symbol,barcode,snv) %>%
#   dplyr::group_by(symbol,barcode) %>%
#   dplyr::mutate(snv = mean(snv)) %>%
#   unique() %>%
#   dplyr::ungroup() %>%
#   tidyr::spread(key = barcode,value = snv)-> snv_merge_snv_data.combine
# snv_merge_snv_data.matrix.combine <- as.matrix(snv_merge_snv_data.combine[,-1])
# rownames(snv_merge_snv_data.matrix.combine) <- snv_merge_snv_data.combine$symbol

gene_list_expr %>% 
  dplyr::filter(barcode %in% samples.intersect) %>%
  dplyr::select(symbol,barcode,expr) %>%
  dplyr::group_by(symbol,barcode) %>%
  dplyr::mutate(expr = mean(expr)) %>%
  unique() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key = barcode,value = expr)-> gene_list_expr_with_mutation_load.combine
gene_list_expr_with_mutation_load.matrix.combine <- as.matrix(gene_list_expr_with_mutation_load.combine[,-1])
rownames(gene_list_expr_with_mutation_load.matrix.combine) <- gene_list_expr_with_mutation_load.combine$symbol

time_status %>%
  dplyr::filter(barcode %in% samples.intersect) %>%
  dplyr::select(barcode,PFS.time) %>%
  dplyr::arrange(barcode) %>%
  unique() %>% .$PFS.time ->  time.combine
time_status %>%
  dplyr::filter(barcode %in% samples.intersect) %>%
  dplyr::select(barcode,PFS) %>%
  dplyr::arrange(barcode) %>%
  unique() %>% .$PFS -> status.combine


# Tumor purity ------------------------------------------------------------

tumor_purity <- readr::read_tsv("/project/huff/huff/data/TCGA_tumor_purity_from_ABSOLUTE/syn_Purity_Ploidy_All_Samples_4-17-15_syn3582761.txt") %>%
  dplyr::rename("barcode"="individual_id")

tumor_purity %>%
  dplyr::select(barcode,tumor_type,purity,ploidy)-> tumor_purity

# save work space ---------------------------------------------------------

save.image(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))

