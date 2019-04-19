############################
# GSVA score
############################
library(GSVA)
# data path ---------------------------------------------------------------
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score")

# load data ---------------------------------------------------------------
exp_data <- readr::read_rds(file.path(TCGA_path,"pancan33_expr.rds.gz"))
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header = T)
gene_list$symbol <- as.character(gene_list$symbol)


# calculation of GSVA score -----------------------------------------------

fn_GSVA <- function(exp){
  genelist <- list(ICPs = gene_list$GeneID)
  exp <- as.matrix(exp)
  rownames.exp <- exp[,2]
  exp <- exp[,-c(1:2)]
  exp <- apply(exp, 2, as.numeric)
  rownames(exp) <- rownames.exp
  res.gsva <- gsva(exp,genelist, mx.diff = FALSE, 
                   method = c("gsva"),
                   kcdf = c("Gaussian"), 
                   verbose = FALSE, parallel.sz = 1)
  res.gsva %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    tidyr::gather(key = "barcode",value = "GSVA score")
}

exp_data %>%
  dplyr::mutate(GSVA = purrr::map(data,fn_GSVA)) %>%
  dplyr::select(-data) -> GSVA.score

GSVA.score %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score.rds.gz"))