###############################################
# mutation load changes between groups

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

snv <- readr::read_rds(file.path(tcga_path, "pancan33_snv.rds.gz"))
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)

# snv difference calculate ------------------------------------------------
## fisher exact test the difference of snv freq of genes between mutation groups
### data prepare
snv %>%
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=snv) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(! is.na(snv)) %>%
  dplyr::group_by(symbol,mutation_status) %>%
  dplyr::mutate(snv_count=sum(snv)) %>%
  dplyr::select(symbol,mutation_status,snv_count) %>%
  dplyr::ungroup() %>%
  unique() %>%
  tidyr::spread(key=mutation_status,value=snv_count) %>%
  dplyr::mutate(high_sample = 1170-high_muation_burden,low_sample=6388-low_mutation_burden) %>%
  dplyr::select(symbol,high_muation_burden ,high_sample ,low_mutation_burden ,low_sample)-> gene_list_snv_count_between_mutation_groups

### function to do fisher exact test
fn_fisher <- function(.x){
  .x %>%
    as.data.frame() %>%
    .[1,] %>%
    unlist() %>%
    matrix(nc=2) -> .matrix
  fisher.test(.matrix) %>%
    broom::tidy() %>%
    dplyr::mutate(p.adjust=p.adjust(p.value,method = "bonferroni"))
}

### get result
gene_list_snv_count_between_mutation_groups %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(snv_diff = purrr::map(data,fn_fisher)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> gene_list_mutation_diff_result_1

gene_list_mutation_diff_result_1 %>%
  dplyr::filter(p.adjust<=0.05) %>%
  dplyr::arrange(p.value)

## snv difference in these 26 cancers
snv %>%
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=snv) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(! is.na(snv)) %>%
  dplyr::select(-sm_count) %>%
  dplyr::group_by(symbol,cancer_types.x,mutation_status) %>%
  dplyr::mutate(snv_count=sum(snv),n=n()) %>%
  dplyr::mutate(per=snv_count/n) %>%
  dplyr::select(symbol,mutation_status,snv_count,n,per) %>%
  dplyr::ungroup() %>%
  unique() -> genelist_snv_per_in_groups_cancers


cancers_with_enough_high_mutation <- c("LUAD","SKCM","BLCA","COAD","LUSC","STAD","HNSC","UCEC","BRCA","LIHC","READ","ACC")
genelist_snv_per_in_groups_cancers %>%
  dplyr::filter(cancer_types.x %in% cancers_with_enough_high_mutation) %>%
  dplyr::filter(mutation_status == "high_muation_burden") %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(sum=sum(per)) %>%
  dplyr::select(symbol,sum) %>%
  unique() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sum) -> gene_rank
genelist_snv_per_in_groups_cancers %>%
  dplyr::filter(cancer_types.x %in% cancers_with_enough_high_mutation) %>%
  dplyr::filter(mutation_status == "high_muation_burden") %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(sum=sum(per)) %>%
  dplyr::select(cancer_types.x,sum) %>%
  unique() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sum) -> cancer_rank

genelist_snv_per_in_groups_cancers %>%
  dplyr::filter(cancer_types.x %in% cancers_with_enough_high_mutation) %>%
  ggplot(aes(x=cancer_types.x,y=symbol,fill=per)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low="#F0FFF0",high="#008B00") +
  facet_grid(~mutation_status) +
  scale_x_discrete(limit=cancer_rank$cancer_types.x) +
  scale_y_discrete(limit = gene_rank$symbol)
