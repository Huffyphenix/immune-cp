###############################################
# Methylation changes between groups

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

mthy <- readr::read_rds(file.path(tcga_path,"pancan33_meth.rds.gz"))
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_type <- read.table(file.path(gene_list_path, "checkpoint.type"),header=T)


# methy difference calculate ------------------------------------------------
## t test to get the difference of methylation levels of genes between mutation groups
### data prepare
mutation_burden_class$cancer_types %>% unique() ->cancers_in_mutaion_burden_class
mthy %>%
  dplyr::filter(cancer_types %in% cancers_in_mutaion_burden_class) %>%    # some cancers in TCGA merged serveral cancers into one
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  dplyr::select(-gene) %>%
  tidyr::gather(-cancer_types,-symbol,key=barcode,value=methy) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::inner_join(mutation_burden_class,by="barcode") %>%
  dplyr::filter(! is.na(methy)) %>%
  tidyr::nest(-symbol) -> genelist_methy_mutaion_class

### function to calculate difference

calculate_fc_pvalue <- function(.x) {
  # pvalue & fdr
  .x %>%
    dplyr::filter(! mutation_status=="Normal") -> .x
  broom::tidy(t.test(methy ~ mutation_status, data = .x)) %>%
    dplyr::mutate(`log2FC(h/l)`=log2(estimate1/estimate2)) %>%
    dplyr::mutate(high_mutation_burden_methy = estimate1,low_mutation_burden_methy = estimate2) %>%
    dplyr::select(high_mutation_burden_methy,low_mutation_burden_methy,`log2FC(h/l)`,p.value) -> df_pvalue
  return(df_pvalue)
}

### get result
genelist_methy_mutaion_class %>%
  dplyr::mutate(diff = purrr::map(data,calculate_fc_pvalue)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::filter(p.value <= 0.05) -> genelist_methy_diff_in_mutaion_groups.sig

genelist_methy_diff_in_mutaion_groups.sig %>%
  readr::write_tsv(file.path(result_path,"4.methy","genelist_methy_diff_in_mutaion_groups.sig"))

## specific gene violin plot 
genelist_methy_mutaion_class %>%
  dplyr::filter(symbol %in% c("BTNL9","ICOSLG","PDCD1LG2","BTN2A1")) %>%
  tidyr::unnest() %>%
  dplyr:::mutate(mutation_status=ifelse(mutation_status=="low_mutation_burden","Low","High")) %>%
  dplyr::rename("Mutation Load" = "mutation_status") %>%
  ggplot(aes(x= `Mutation Load`,y=methy)) +
  geom_violin(trim=F,aes(fill= `Mutation Load`)) +
  geom_boxplot(width=0.1) +
  theme(
    panel.background = element_blank(),
    panel.border =  element_rect(fill = NA)
  ) +
  facet_grid(~symbol)


