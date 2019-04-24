############################
# GSVA score
# all 79 ICP genes as one feature
############################
library(GSVA)
library(magrittr)
library(ggplot2)
library(survival)

# data path ---------------------------------------------------------------
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/2.GSVA-ICPs_exp_site_5_feature")

# load data ---------------------------------------------------------------
exp_data <- readr::read_rds(file.path(TCGA_path,"pancan33_expr.rds.gz"))
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header = T)
gene_list$symbol <- as.character(gene_list$symbol)
survival_data <- readr::read_rds(file.path("/home/huff/project/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz")) %>%
  dplyr::rename("cancer_types" = "type")

ICPs_exp_in_TI.cancer.site_groups <- 
  readr::read_tsv(file.path(res_path, "ICPs_exp_in_TI.cancer.site_groups.tsv")) %>%
  dplyr::inner_join(gene_list,by = "symbol") %>%
  dplyr::select(symbol, Tissues, `log2FC(I/T)`, TCGA_Cancer, ICP_exp_group, GeneID, type, functionWithImmune)

source("/home/huff/project/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# calculation of GSVA score -----------------------------------------------
fn_get_gene_list_feature   <- function(cancer_types){
  ICPs_exp_in_TI.cancer.site_groups %>%
    dplyr::filter(TCGA_Cancer %in% cancer_types) -> .data
  
  .group <- unique(.data$ICP_exp_group)
  
  .genelist <- list()
  for (group in .group){
    .data %>%
      dplyr::filter(ICP_exp_group %in% group) %>%
      .$GeneID -> .geneid
    .genelist[[group]] <- .geneid
  }
  .genelist
}

fn_GSVA <- function(cancer_types, exp){
  print(paste("GSVA",cancer_types))
  
  genelist <- fn_get_gene_list_feature(cancer_types)
  if (length(genelist)>0) {
    index <- which(substr(colnames(exp),14,14) == "0")
    exp <- exp[,c(2,index)]
    
    exp <- as.matrix(exp)
    rownames.exp <- exp[,1]
    exp <- exp[,-1]
    exp <- apply(exp, 2, as.numeric)
    rownames(exp) <- rownames.exp
    res.gsva <- gsva(exp,genelist, mx.diff = FALSE, 
                     method = c("gsva"),
                     kcdf = c("Gaussian"), 
                     verbose = FALSE, parallel.sz = 1)
    res.gsva %>%
      as.data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::mutate(feature = rownames(res.gsva)) %>%
      tidyr::gather(-feature, key = "barcode",value = "GSVA_score") %>%
      tidyr::spread(key = "feature", value = "GSVA_score") -> gsva.score
  } else{
    gsva.score <- tibble::tibble()
  }
  
  
  # fn_compare_TN(gsva.score, cancer_types, result_path = file.path(res_path,"TN_compare"))
  
  # fn_survival.calculate(gsva.score, cancer_types)
  
  # fn_heatmap(res.gsva)
  gsva.score
}



exp_data %>%
  # head(1) %>%
  dplyr::mutate(GSVA = purrr::map2(cancer_types, expr, fn_GSVA)) %>%
  dplyr::select(-expr) -> GSVA.score

GSVA.score %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score.rds.gz"), compress = "gz")
