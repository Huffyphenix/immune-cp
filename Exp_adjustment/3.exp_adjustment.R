
# Ultilize of FANTOM 5 data and tumor purity data to adjust ICP expression in TCGA ---------------------------------

library(magrittr)
# data path ---------------------------------------------------------------

immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")
purity_path <- "/project/huff/huff/data/TCGA_tumor_purity_from_ncomms9971"
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")

# load data ---------------------------------------------------------------
TCGA_purity <- readr::read_tsv(file.path(purity_path,"ncomms9971-s2.txt")) 
# Fantom expression
ICP_fantom.gene_exp.cell_line.Immune_cell.combine <- readr::read_rds(file.path(immune_path,"genelist_data","ICP_fantom.gene_exp.cell_line.Immune_cell.combine.rds.gz"))
# classification of tcga cancers
TCGA_tissue <- readr::read_tsv("/project/huff/huff/data/TCGA/TCGA_cancer_tissue_classification.txt")
# gene list
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
# xCell data
xCell_TCGA_RSEM.immune_stroma.ratio <- readr::read_rds(file.path(immune_path,"genelist_data","Pancan21.tcga.xCell.immune_stroma.ratio.rds.gz")) %>%
  dplyr::select(-barcode) %>%
  dplyr::rename("barcode" = "Sample ID")
# tcga expression data
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

gene_list_expr <- readr::read_rds(file.path(tcga_path,"pancan33_expr.rds.gz")) %>%
  dplyr::mutate(filter_expr = purrr::map(expr, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-expr) 

ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
  dplyr::filter(`Characteristics[Tissue]` %in% TCGA_tissue$Tissues) %>%
  .$`Characteristics[Tissue]` %>%
  unique() -> TCGA_Fantom_share_tissue

# function to adjust expression ----------------------------------------------------------------
fn_adjust_exp <- function(.cancer,.x){
  TCGA_tissue %>%
    dplyr::filter(TCGA_Cancer %in% .cancer) %>%
    unique() -> .tissue
  print(.cancer)
  print(.tissue$Tissues)
  if(is.na(grep(.tissue$Tissues,TCGA_Fantom_share_tissue))){
    return(tibble::tibble())
  }else{
    ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
      dplyr::filter(`Characteristics[Tissue]` %in% c("PrimaryCell",.tissue$Tissues)) %>%
      dplyr::select(symbol,Group,gene_mean_exp) %>%
      tidyr::spread(key="Group",value="gene_mean_exp") %>%
      dplyr::rename("I_F"="Immune Cell","S_F" = "Stromal Cell","T_F"="Tumor Cell") -> .fantom.exp
    
    .x %>%
      dplyr::select(-entrez_id) %>%
      tidyr::gather(-symbol,key="barcode",value=exp) %>%
      dplyr::mutate(barcode=substr(barcode,1,15)) %>%
      dplyr::inner_join(.fantom.exp,by="symbol") %>%
      dplyr::inner_join(xCell_TCGA_RSEM.immune_stroma.ratio,by="barcode") %>%
      dplyr::rename("T_X"="CPE","I_X"="Immune_ratio","S_X"="Stromal_ratio") %>%
      dplyr::mutate(S_X=S_X+other_Stromal_ratio) %>%   # hold the assumption that little difference od genes' expression in stromal and other stroma cell.
      dplyr::select(symbol,barcode,exp,"I_F","S_F","T_F","T_X","I_X","S_X")-> .tmp_data
    
    .tmp_data %>%
      tidyr::nest(-symbol,-barcode) %>%
      dplyr::group_by(symbol,barcode) %>%
      dplyr::mutate(exp_adjust = purrr::map(data,fn_adjust_calculate)) %>%
      tidyr::unnest()
  }
}
fn_adjust_calculate <- function(.x){
    exp <- .x$exp
    IF<-.x$I_F # Gene expression in immune cell of FANTOM 5 sample.
    IX<-.x$I_X # Ratio of immune cell in xCell predict TCGA sample.
    TF<-.x$T_F # Gene expression in tumor cell of FANTOM 5 sample.
    TX<-.x$T_X # Ratio of tumor cell in xCell predict TCGA sample.
    SF<-.x$S_F # Gene expression in stromal cell of FANTOM 5 sample.
    SX<-.x$S_X # Ratio of stromal cell in xCell predict TCGA sample.
    
    IT<-exp*IF*IX/(TF*TX*TX+IF*IX*IX+SF*SX*SX) # Gene expression in immune cell of TCGA sample.
    TT<-exp*TF*TX/(TF*TX*TX+IF*IX*IX+SF*SX*SX) # Gene expression in tumor cell of TCGA sample.
    ST<-exp*SF*SX/(TF*TX*TX+IF*IX*IX+SF*SX*SX) # Gene expression in stromal cell of TCGA sample.
    
    predic_exp_from_celldata <- IF*IX+TF*TX+SF*SX # to compare with TCGA bulk exp
    tibble::tibble(T_T=TT,I_T=IT,S_T=ST,predic_exp_from_celldata=predic_exp_from_celldata)
}


# Running exp adjustment --------------------------------------------------

library(multidplyr)
library(parallel)
cl<-21
cluster <- create_cluster(cores = cl)

gene_list_expr %>%
  dplyr::filter(cancer_types %in% unique(TCGA_purity$`Cancer type`))  %>%
  partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("tidyverse") %>%
  multidplyr::cluster_library("stringr") %>%
  multidplyr::cluster_assign_value("fn_adjust_exp",fn_adjust_exp) %>%
  multidplyr::cluster_assign_value("fn_adjust_calculate",fn_adjust_calculate) %>%
  multidplyr::cluster_assign_value("TCGA_tissue",TCGA_tissue) %>%
  multidplyr::cluster_assign_value("TCGA_Fantom_share_tissue",TCGA_Fantom_share_tissue) %>%
  multidplyr::cluster_assign_value("ICP_fantom.gene_exp.cell_line.Immune_cell.combine",ICP_fantom.gene_exp.cell_line.Immune_cell.combine) %>%
  multidplyr::cluster_assign_value("xCell_TCGA_RSEM.immune_stroma.ratio",xCell_TCGA_RSEM.immune_stroma.ratio) %ta>%
  dplyr::mutate(adjust_exp = purrr::map2(cancer_types,filter_expr,fn_adjust_exp)) %>%
  collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -filter_expr) -> gene_list_adjust_expr
on.exit(parallel::stopCluster(cluster))

  
gene_list_adjust_expr %>%
  readr::write_tsv(file.path(gene_list_path,"pancan21.ICP.exp_adjust.by_cell_ratio.rds.gz"),compress="gz")
