### mean expression of each gene features in TCGA
library(magrittr)
library(tidyverse)

# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/home/huff/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")

res_path <- file.path(immune_res_path,"ICP_score.new")

# load data ---------------------------------------------------------------

# expression data
expr <- readr::read_rds(file.path(TCGA_path, "pancan33_expr.rds.gz"))

gene_list <- readr::read_tsv(file.path(gene_list_path, "ICPs_all_info_class-new.tsv")) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("N"),"Not_sure",Exp_site)) 

# 1.get mean expression of all possible features of ICP ----------------------------------------------
# 1.1. get gene feature from gene list ----
genelist <- list()
#### gene list feature by gene exp site #####
for(expsite in c("Tumor cell dominate","Immune and tumor cell almost","Immune cell dominate")){
  genelist[[expsite]] <- gene_list %>%
    dplyr::filter(Exp_site.1 == expsite) %>%
    .$symbol
}

#### gene list feature by gene function in immune system #####
for(fun_class in unique(gene_list$functionWithImmune)){
  genelist[[fun_class]] <- gene_list %>%
    dplyr::filter(functionWithImmune == fun_class) %>%
    .$symbol
}
#### gene list feature by ligand-receptor pairs #####
for(p in unique(gene_list$Recepter_pairs)){
  if(!is.na(p)){
    genelist[[p]] <- gene_list %>%
      dplyr::filter(Recepter_pairs == p) %>%
      .$symbol
  }
  
}
#### gene list feature by gene family #####
for(f in unique(gene_list$family)){
  if(f != "Other"){
    genelist[[f]] <- gene_list %>%
      dplyr::filter(family == f) %>%
      .$symbol
  }
}
genelist[["All_gene"]] <- gene_list$symbol

# 1.2 funciton to get mean expression ----
# fn_mean <- function(names,genelist,exp){
#   genes <- genelist[[names]]
#   exp %>%
#     dplyr::filter(symbol %in% genes) %>%
#     dplyr::group_by(barcode) %>%
#     dplyr::mutate(value = mean(exp)) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(Features = paste("Mean.",names,sep=""))%>%
#     dplyr::select(barcode,value,Features) %>%
#     unique()
# }
# 
expr %>%
  dplyr::mutate(exp_deal = purrr::map(expr,.f=function(.x){
    .x %>%
      filter(symbol %in% gene_list$symbol) %>%
      dplyr::select(-entrez_id) %>%
      tidyr::gather(-symbol,key="barcode",value="exp")
  })) %>%
  dplyr::select(-expr) %>%
  tidyr::unnest() -> expr_ready

# 1.3 calculation of mean exp-----
# tibble::tibble(names = names(genelist)) %>%
#   dplyr::mutate(mean = purrr::map(names,fn_mean,genelist=genelist,exp=expr_ready)) -> TCGA_mean_exp

# 2.get fold expression of all possible features of ICP ----------------------------------------------
# 2.1 function to get ratio of (geneA+geneB+..)/geneA ------
fn_ratio_A <- function(names,genelist,exp){
  genes <- genelist[[names]]
  exp %>%
    dplyr::filter(symbol %in% genes) %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(value = sum(exp)) %>%
    dplyr::ungroup() %>%
    dplyr::select(barcode,value) %>%
    unique() -> exp.sum
  exp %>%
    dplyr::filter(symbol %in% genes) %>%
    dplyr::inner_join(exp.sum,by="barcode") %>%
    dplyr::mutate(value = (value+1)/(exp+1))  %>%
    dplyr::mutate(Features = paste("Fold.",names,".",symbol,sep="")) %>%
    dplyr::select(barcode,value,Features) %>%
    unique()
}

# 2.2 calculation -------
tibble::tibble(names = names(genelist)) %>%
  dplyr::mutate(fold = purrr::map(names,fn_ratio_A,genelist=genelist,exp=expr_ready)) -> TCGA_fold_exp

# 3.get ratio expression of all possible features of ICP ----------------------------------------------
# 2.1 function to get ratio of geneA/geneB ------
fn_ratio_B <- function(names,genelist,exp){
  genes <- genelist[[names]]
  combn(genes,2) %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    tidyr::gather(key="combn",value="symbol") %>% 
    tidyr::nest(-combn) -> .combn
  
  .combn %>%
    dplyr::mutate(ratio = purrr::map(data,.f=function(.x){
      exp %>%
        dplyr::filter(symbol %in% .x$symbol) %>%
        tidyr::spread(key="symbol",value="exp") %>%
        dplyr::rename("G1"=.x$symbol[1],"G2"=.x$symbol[2]) %>%
        dplyr::mutate(value = (G1+1)/(G2+1)) %>%
        dplyr::mutate(Features = paste("Ratio.",.x$symbol[1],".",.x$symbol[2],sep=""))%>%
        dplyr::select(barcode,value,Features) %>%
        unique()
    })) %>%
    dplyr::select(-data,-combn) %>%
    tidyr::unnest() 
}

# 2.2 calculation -------
lapply(genelist, function(x){length(x)}) %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::as.tbl() %>% 
  tidyr::gather(key="names",value="count") %>% 
  tidyr::unnest() %>%
  dplyr::filter(count < 5) %>%
  .$names -> names 

tibble::tibble(names = names) %>%
  dplyr::mutate(fold = purrr::map(names,fn_ratio_B,genelist=genelist,exp=expr_ready)) -> TCGA_ratio_exp


# 4.combine all score -----------------------------------------------------

TCGA_ratio_exp %>%
  tidyr::unnest() %>%
  rbind(TCGA_fold_exp %>% tidyr::unnest()) %>%
  # rbind(TCGA_mean_exp %>% tidyr::unnest()) %>%
  dplyr::select(-names) %>%
  tidyr::spread(key="Features",value="value") -> TCGA_mean_fold_ratio_features_value

TCGA_mean_fold_ratio_features_value %>%
  readr::write_rds(file.path(res_path,"new-TCGA_mean_fold_ratio_features_value.rds.gz"),compress = "gz")

