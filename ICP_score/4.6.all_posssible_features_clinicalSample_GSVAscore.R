############################
# GSVA score
# loggic regression to select efficient features to do GSVA score
# all posible feature I can imagine 
# features constructed by 79 ICP genes classified into groups by exp site, and these groups are used as features to do GSVA
# loggic regression
############################

library(GSVA)
library(magrittr)
library(ggplot2)
library(survival)

# data path ---------------------------------------------------------------
# server 1 basic path
basic_path <- file.path("/home/huff/project")

# E zhou basi path
basic_path <- file.path("F:/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/2.1.Clinical_validation-GSVA-ICPs_exp_site_5_feature")

# load data ---------------------------------------------------------------
exp_data <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data/mRNA_exp","all_FPKM_expression_2.txt"))
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header = T)
gene_list$symbol <- as.character(gene_list$symbol)


ICP_family <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/checkpoint/ICP_gene_family.txt"))%>%
  dplyr::inner_join(gene_list, by = "symbol")
ICP_ligand_receptor_pair <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/checkpoint/ICP_gene_ligand_receptor_pairs.txt")) %>%
  dplyr::mutate(pairs = paste("pair",pairs)) %>%
  dplyr::inner_join(gene_list, by = "symbol")

# source("/home/huff/project/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# exp data classfication --------------------------------------------------

sample_info <-
  readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv"))%>%
  dplyr::filter(Library_strategy == "RNA-Seq") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment")

library(clusterProfiler)
library(org.Hs.eg.db)
bitr(gene_list$symbol,"SYMBOL", "ENSEMBL",OrgDb = org.Hs.eg.db)

sample_info %>%
  dplyr::select(Run, Cancer.y, blockade) %>%
  unique() -> Run_pubmed.id
exp_data %>%
  tidyr::gather(-gene_id, key="Run", value = exp) %>%
  dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
  dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
  tidyr::nest(-Cancer.y, -blockade) -> exp_data.gather

exp_data.gather %>%
  dplyr::mutate(data_spread = purrr::map(data,.f = function(.x){
    .x %>%
      tidyr::spread(key="Run",value="exp")
  })) %>%
  dplyr::select(-data) -> exp_data.nest

# gene id transfer --------------------------------------------------------
ensembl_gene <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","Homo_sapiens.gene_info.geneid.symbol.ensembl"))

gene_list %>%
  dplyr::inner_join(ensembl_gene, by = "GeneID") %>%
  dplyr::filter(!is.na(ens_id)) -> gene_list

ICPs_exp_in_TI.cancer.site_groups <- 
  readr::read_tsv(file.path(immune_res_path, "ICP_score/2.GSVA-ICPs_exp_site_5_feature","ICPs_exp_in_TI.cancer.site_groups.tsv")) %>%
  dplyr::inner_join(gene_list,by = "symbol") %>%
  dplyr::select(symbol, Tissues, `log2FC(I/T)`, TCGA_Cancer, ICP_exp_group, GeneID, ens_id, type, functionWithImmune)

# calculation of GSVA score -----------------------------------------------
fn_get_gene_list_feature_bysite  <- function(tissue){
  ICPs_exp_in_TI.cancer.site_groups %>%
    dplyr::filter(Tissues %in% tissue) -> .data
  
  .group <- unique(.data$ICP_exp_group)
  
  .genelist <- list()
  for (group in .group){
    .data %>%
      dplyr::filter(ICP_exp_group %in% group) %>%
      .$ens_id -> .geneid
    .genelist[[group]] <- .geneid
  }
  .genelist
}

fn_GSVA <- function(tissue, exp){
  print(paste("GSVA",tissue))
  
  #### gene list feature by gene exp site #####
  genelist <- fn_get_gene_list_feature_bysite(tissue)
  #### gene list feature by gene function in immune system #####
  for(fun_class in unique(gene_list$functionWithImmune)){
    genelist[[fun_class]] <- gene_list %>%
      dplyr::filter(functionWithImmune == fun_class) %>%
      .$ens_id
  }
  #### gene list feature by ligand-receptor pairs #####
  for(p in unique(ICP_ligand_receptor_pair$pairs)){
    genelist[[p]] <- ICP_ligand_receptor_pair %>%
      dplyr::filter(pairs == p) %>%
      .$ens_id
  }
  #### gene list feature by gene family #####
  for(f in unique(ICP_family$family)){
    genelist[[f]] <- ICP_family %>%
      dplyr::filter(f == f) %>%
      .$ens_id
  }
  #### do GSVA ####
  if (length(genelist)>0) {
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
      tidyr::gather(-feature, key = "Run",value = "GSVA_score") %>%
      tidyr::spread(key = "feature", value = "GSVA_score") -> gsva.score
  } else{
    gsva.score <- tibble::tibble()
  }
  
  
  # fn_compare_TN(gsva.score, cancer_types, result_path = file.path(res_path,"TN_compare"))
  
  # fn_survival.calculate(gsva.score, cancer_types)
  
  # fn_heatmap(res.gsva)
  gsva.score
}



exp_data.nest %>%
  dplyr::mutate(tissue = ifelse(Cancer.y=="gastric cancer", "stomach", "skin")) %>%
  # head(1) %>%
  dplyr::mutate(GSVA = purrr::map2(tissue, data_spread, fn_GSVA)) %>%
  dplyr::select(-data_spread) -> GSVA.score

GSVA.score %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score-by_tissue_target-all-possible-features.rds.gz"), compress = "gz")
