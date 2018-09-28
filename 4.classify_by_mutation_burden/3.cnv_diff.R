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

cnv <- readr::read_rds(file.path(tcga_path,"pancan34_cnv_threshold.rds.gz"))
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)


# cnv difference calculate ------------------------------------------------
## fisher exact test the difference of high cnv (-2 or 2) freq of genes between mutation groups
### data prepare
mutation_burden_class$cancer_types %>% unique() ->cancers_in_mutaion_burden_class
cnv %>%
  dplyr::filter(cancer_types %in% cancers_in_mutaion_burden_class) %>%    # some cancers in TCGA merged serveral cancers into one
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=cnv) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(! is.na(cnv)) -> cnv_merge_snv_data

cnv_merge_snv_data %>% dplyr::select(barcode,mutation_status) %>% unique() %>% .$mutation_status %>% table() -> cnv_merge_snv_sample_info

cnv_merge_snv_data %>%
  dplyr::mutate(cnv = ifelse(abs(cnv)<2,0,cnv/2)) %>%
  dplyr::mutate(cnv_type = ifelse(cnv < 0, "Del", "Amp")) %>%
  dplyr::group_by(symbol,mutation_status,cnv_type) %>%
  dplyr::mutate(cnv_count=sum(cnv)) %>%
  dplyr::select(symbol,mutation_status,cnv_type,cnv_count) %>%
  dplyr::ungroup() %>%
  unique() %>%
  tidyr::spread(key=mutation_status,value=cnv_count) %>%
  dplyr::mutate(high_muation_burden = ifelse(is.na(high_muation_burden),0,high_muation_burden)) %>%
  dplyr::mutate(low_mutation_burden = ifelse(is.na(low_mutation_burden),0,low_mutation_burden)) %>%
  dplyr::mutate(high_sample = cnv_merge_snv_sample_info[1]-high_muation_burden,low_sample=cnv_merge_snv_sample_info[2]-low_mutation_burden) %>%
  dplyr::select(symbol,cnv_type,high_muation_burden ,high_sample ,low_mutation_burden ,low_sample) -> gene_list_cnv_count_between_mutation_groups

### function to do fisher exact test
fn_fisher <- function(.x,.y){
  print(paste(.y))
  # print(.z)
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
gene_list_cnv_count_between_mutation_groups %>%
  dplyr::mutate(high_muation_burden=ifelse(high_muation_burden>0,high_muation_burden,-high_muation_burden)) %>%
  dplyr::mutate(low_mutation_burden=ifelse(low_mutation_burden>0,low_mutation_burden,-low_mutation_burden)) %>%
  tidyr::nest(-symbol,-cnv_type) %>%
  dplyr::mutate(cnv_diff = purrr::map2(.x=data,.y=symbol,fn_fisher)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> gene_list_cnv_diff_result

gene_list_cnv_diff_result_1 %>%
  dplyr::filter(p.adjust<=0.05) %>%
  dplyr::arrange(p.value) -> gene_list_cnv_diff_result_sig

## cnv difference in these 26 cancers
cnv_merge_snv_data %>%
  dplyr::select(cancer_types.x,barcode) %>%
  unique() %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(cancer_types.x,n) %>%
  unique() -> cancers_number
cnv_merge_snv_data %>%
  dplyr::inner_join(cancers_number,by="cancer_types.x") %>%
  dplyr::mutate(cnv = ifelse(abs(cnv)<2,0,cnv/2)) %>%
  dplyr::mutate(cnv_type = ifelse(cnv==1,"Amp","Del")) %>%
  dplyr::mutate(cnv_type = ifelse(cnv==0,"None",cnv_type)) %>%
  dplyr::filter(! cnv_type == "None") %>%
  dplyr::group_by(symbol,cancer_types.x,mutation_status,cnv_type) %>%
  dplyr::mutate(cnv_count=sum(cnv)) %>%
  dplyr::mutate(per=100*cnv_count/n) %>%
  dplyr::select(symbol,mutation_status,cnv_count,n,per) %>%
  dplyr::ungroup() %>%
  unique() -> genelist_cnv_per_in_groups_cancers

cancers_with_enough_high_mutation <- c("LUAD","SKCM","BLCA","COAD","LUSC","STAD","HNSC","UCEC","BRCA","LIHC","READ","ACC")
genelist_cnv_per_in_groups_cancers %>%
  dplyr::filter(cancer_types.x %in% cancers_with_enough_high_mutation) %>%
  dplyr::filter(mutation_status == "high_muation_burden") %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(sum=sum(per)) %>%
  dplyr::select(symbol,sum) %>%
  unique() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sum) -> gene_rank
genelist_cnv_per_in_groups_cancers %>%
  dplyr::filter(cancer_types.x %in% cancers_with_enough_high_mutation) %>%
  dplyr::filter(mutation_status == "high_muation_burden") %>%
  dplyr::group_by(cancer_types.x) %>%
  dplyr::mutate(sum=sum(per)) %>%
  dplyr::select(cancer_types.x,sum) %>%
  unique() %>%
  dplyr::ungroup() %>%
  dplyr::arrange(sum) -> cancer_rank

genelist_cnv_per_in_groups_cancers %>%
  dplyr::filter(cancer_types.x %in% cancers_with_enough_high_mutation) %>%
  ggplot(aes(x=cancer_types.x,y=symbol,fill=per)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(mid="white",low="#008B00",high = "red",
                       na.value = 0,
                       name = "CNV Frequency (%)") +
  facet_grid(cnv_type~mutation_status) +
  scale_x_discrete(limit=cancer_rank$cancer_types.x) +
  scale_y_discrete(limit = gene_rank$symbol) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = grid::unit(c(0, 0, 0, 0), "mm"),
    panel.border  = element_rect(fill=NA),
    axis.text.x = element_text(angle = 45,hjust = 1)
) -> p


    