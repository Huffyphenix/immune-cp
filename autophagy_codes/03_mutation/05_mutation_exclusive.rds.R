library(magrittr)
library(ggplot2)
out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")
tcga_path <- "/home/huff/project/immune_checkpoint/data/TCGA_data"
expr_path <-c("/home/huff/project/immune_checkpoint/result_20171025/expr_rds")

# load cnv and gene list
snv_syn7824274 <- readr::read_rds(file.path("/data/shiny-data/GSCALite/TCGA/snv","pancan33_snv_from_syn7824274_gather.rds.gz"))
snv <- readr::read_rds(file.path(tcga_path, "pancan33_snv.rds.gz"))

basic_path <- "/home/huff/project/immune_checkpoint"
gene_list_path <- file.path(basic_path,"checkpoint/20171021_checkpoint")
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

fn_merge_data <- function(.data){
  .data %>%
    tidyr::gather(-symbol,key = sample, value = count) %>%
    dplyr::left_join(snv_syn7824274,by=c("sample","symbol")) %>%
    dplyr::select(symbol,sample,mut_n) %>%
    tidyr::spread(key = sample, value = mut_n) 
}
snv %>%
  dplyr::mutate(snv_syn7824274 = purrr::map(snv,fn_merge_data)) %>%
  dplyr::select(-snv) -> snv_syn7824274.with_full_samples


# filter out hypermutation samples -----
burden_path <- "/project/huff/huff/data/TCGA"
mutation_burden_class <- readr::read_rds(file.path(burden_path,"classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest()

# gene_list_snv %>%
#   tidyr::unnest() %>%
#   tidyr::gather(-cancer_types,-symbol,key="barcode",value="count") %>%
#   dplyr::left_join(mutation_burden_class,by="barcode") %>%
#   dplyr::filter(!sm_count > 1000) %>%
#   dplyr::select(-Cancer_Types) -> gene_list_snv.filter_hypermutation 

fn_get_percent <- function(cancer_types, filter_snv){
  print(cancer_types)
  n <- length(filter_snv) - 1
  filter_snv %>%
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .))) -> .d
  .d %>% 
    tidyr::gather(-symbol,key = barcode, value = count) %>% 
    dplyr::left_join(mutation_burden_class,by="barcode") %>%
    dplyr::filter(!sm_count > 1000) %>%
    dplyr::mutate(samples = ifelse(count > 0, 1, 0)) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::summarise(sm_count = sum(count), sm_sample = sum(samples)) %>% 
    dplyr::mutate(per = sm_sample / n) %>%
    dplyr::mutate(n=n)-> .d_count
  
  return(.d_count)
  # tibble::tibble( mut_count = list(.d_count))
}
fn_get_hyper_mut_data_class <- function(cancer_types, filter_snv){
  print(cancer_types)
  n <- length(filter_snv) - 1
  filter_snv %>%
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .))) -> .d
  .d %>% 
    tidyr::gather(key = barcode, value = count, -symbol) %>% 
    dplyr::left_join(mutation_burden_class,by="barcode") %>%
    dplyr::mutate(class = ifelse(sm_count > 1000,"hyper","non-hyper")) 
}

# filter out hypermutation samples -----

snv_syn7824274.with_full_samples %>% 
  # head(1) %>% .$filter_snv_syn7824274 %>% .[[1]] -> filter_snv
  plyr::mutate(res = purrr::map2(cancer_types, snv_syn7824274, fn_get_percent)) %>%
  dplyr::select(-snv_syn7824274) -> gene_list_snv.count_per.filter_hypermutation

snv_syn7824274.with_full_samples %>% 
  # head(1) %>% .$filter_snv %>% .[[1]] -> filter_snv
  plyr::mutate(res = purrr::map2(cancer_types, snv_syn7824274, fn_get_hyper_mut_data_class)) %>%
  dplyr::select(-snv_syn7824274) -> syn7824274.snv.hypermutation_class

# get most mutated genes in cancers -----
get_top_mutated_genes <- function(.x){
  .x %>%
    dplyr::arrange(desc(sm_count)) %>%
    dplyr::filter(per > 0.1)
}

gene_list_snv.count_per.filter_hypermutation %>%
  dplyr::mutate(high_mutated = purrr::map(res,get_top_mutated_genes)) %>%
  dplyr::select(-res) %>%
  tidyr::unnest() -> high_mutated_genes_in_cancers

# GET MUTATION DATA TO DO EXCLUSIVE MUTATION ANALYSIS ----------------------------
snv_syn7824274.with_full_samples %>%
  dplyr::mutate(ICP_SNV = purrr::map(snv_syn7824274,.f=function(.x){
    .x %>%
      dplyr::filter(symbol %in% gene_list$symbol)
  })) %>%
  dplyr::mutate(highGene_SNV = purrr::map2(snv_syn7824274,cancer_types,.f=function(.x,.y){
    high_mutated_genes_in_cancers %>%
      dplyr::filter(cancer_types == .y) %>%
      dplyr::select(symbol) %>%
      dplyr::inner_join(.x,by="symbol")
  })) %>%
  dplyr::select(-snv_syn7824274) -> ICP_highGene_snv

ICP_highGene_snv %>%
  readr::write_rds(file.path(snv_path,"ICP_highGene_snv.rds.gz"))

save.image(file.path(snv_path,".rda_snv_exclusive.rda"))
load(file.path(snv_path,".rda_snv_exclusive.rda"))
