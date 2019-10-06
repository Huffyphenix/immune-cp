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
# basic_path <- file.path("F:/我的坚果云")

# Hust basi path
# basic_path <- file.path("S:/坚果云/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
# res_path <- file.path(immune_res_path,"ICP_score/2.2.Clinical_validation-GSVA-ICPs_exp_site_5_feature") # old result
res_path <- file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")

load(file.path(res_path,"get_clinical_gsvascore.rds"))
# load data ---------------------------------------------------------------
# exp_data.fpkm <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data/mRNA_exp","all_FPKM_expression_2.txt"))
exp_data <- readr::read_rds(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp","Batch_corrected_CPM_exp.rds.gz"))  
exp_data %>%
  as.data.frame() %>%
  dplyr::mutate(gene_id = rownames(exp_data)) -> exp_data
gene_list <- readr::read_tsv(file.path(gene_list_path, "ICPs_all_info_class.tsv")) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Immune","Mainly_exp_on_Immune"),"Mainly_exp_on_Immune",Exp_site)) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("Only_exp_on_Tumor","Mainly_exp_on_Tumor" ),"Mainly_exp_on_Tumor",Exp_site.1)) %>%
  dplyr::mutate(Exp_site.1 = ifelse(Exp_site %in% c("N"),"Not_sure",Exp_site.1)) 

# ICP_family <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/checkpoint/ICP_gene_family.txt"))%>%
#   dplyr::inner_join(gene_list, by = "symbol")
# ICP_ligand_receptor_pair <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/checkpoint/ICP_gene_ligand_receptor_pairs.txt")) %>%
#   dplyr::mutate(pairs = paste("pair",pairs)) %>%
#   dplyr::inner_join(gene_list, by = "symbol")

# source("/home/huff/project/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# exp data classfication --------------------------------------------------

sample_info <-
  readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv"))%>%
  dplyr::filter(Library_strategy == "RNA-Seq") %>%
  dplyr::inner_join(readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","anti_ICB_study_info.txt")) %>%
                      dplyr::select(`data ID`,Author) %>%
                      unique(), by="data ID")

library(clusterProfiler)
library(org.Hs.eg.db)
bitr(gene_list$symbol,"SYMBOL", "ENSEMBL",OrgDb = org.Hs.eg.db)

sample_info %>%
  dplyr::select(Run, Cancer.y, Cancer_type,  blockade,Biopsy_Time,Author) %>%
  unique() -> Run_pubmed.id

sample_info %>%
  dplyr::select(Run,Cancer.y,blockade,blockade,Biopsy_Time,Author,Response) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "Response", "non-Response")) %>%
  dplyr::group_by(Cancer.y,blockade,blockade,Biopsy_Time,Author,Response) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(-Run) %>%
  unique() %>%
  tidyr::spread(key="Response",value="n") %>%
  dplyr::mutate(`Response_percentage(%)` = Response*100/(`non-Response`+Response)) %>%
  readr::write_tsv(file.path(res_path,"Response_statistic_each_dataset.tsv"))

# gene id transfer --------------------------------------------------------
ensembl_gene <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","Homo_sapiens.gene_info.geneid.symbol.ensembl"))

gene_list %>%
  dplyr::inner_join(ensembl_gene, by = "GeneID") %>%
  dplyr::filter(!is.na(ens_id)) -> gene_list

# 1.get GSVA score of all possible features of ICP ----------------------------------------------
# 1.1. get gene feature from gene list ----
genelist <- list()
#### gene list feature by gene exp site #####
for(expsite in c("Both_exp_on_Tumor_Immune","Mainly_exp_on_Immune","Mainly_exp_on_Tumor")){
  genelist[[expsite]] <- gene_list %>%
    dplyr::filter(Exp_site.1 == expsite) %>%
    .$ens_id
}

#### gene list feature by gene function in immune system #####
for(fun_class in unique(gene_list$functionWithImmune)){
  genelist[[fun_class]] <- gene_list %>%
    dplyr::filter(functionWithImmune == fun_class) %>%
    .$ens_id
}
#### gene list feature by ligand-receptor pairs #####
for(p in unique(gene_list$Recepter_pairs)){
  if(!is.na(p)){
    genelist[[p]] <- gene_list %>%
      dplyr::filter(Recepter_pairs == p) %>%
      .$ens_id
  }
  
}
#### gene list feature by gene family #####
for(f in unique(gene_list$family)){
  if(f != "Other"){
    genelist[[f]] <- gene_list %>%
      dplyr::filter(family == f) %>%
      .$ens_id
  }
}
genelist[["All_gene"]] <- gene_list$ens_id

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
  #### do GSVA ####
  if (length(genelist)>0) {
    exp <- as.matrix(exp)
    rownames.exp <- exp[,1]
    exp <- exp[,-1]
    exp <- apply(exp, 2, as.numeric)
    rownames(exp) <- rownames.exp
    exp[is.na(exp)] <- 0
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

# cancer and type specific ----
exp_data %>%
  tidyr::gather(-gene_id, key="Run", value = exp) %>%
  dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
  dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
  dplyr::select(-Biopsy_Time,-Author, -blockade) %>%
  tidyr::nest(-Cancer.y, -Cancer_type) -> exp_data.gather

exp_data.gather %>%
  dplyr::mutate(data_spread = purrr::map(data,.f = function(.x){
    .x %>%
      tidyr::spread(key="Run",value="exp")
  })) %>%
  dplyr::select(-data) -> exp_data.nest

exp_data.nest %>%
  dplyr::mutate(tissue = ifelse(Cancer.y=="gastric cancer", "stomach", "skin")) %>%
  # head(1) %>%
  dplyr::mutate(GSVA = purrr::map2(tissue, data_spread, fn_GSVA)) %>%
  dplyr::select(-data_spread) -> GSVA.score

GSVA.score %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score_all-possible-features_all-cancer_metastic_specific.rds.gz"), compress = "gz")


# cancer specific ----
exp_data %>%
  tidyr::gather(-gene_id, key="Run", value = exp) %>%
  dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
  dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
  dplyr::select(-Cancer_type,-Biopsy_Time,-blockade,-Author) %>%
  tidyr::nest(-Cancer.y) %>%
  dplyr::mutate(data_spread = purrr::map(data,.f = function(.x){
    .x %>%
      tidyr::spread(key="Run",value="exp")
  })) %>%
  dplyr::select(-data) -> exp_data.nest.cancer_specific

exp_data.nest.cancer_specific %>%
  dplyr::mutate(tissue = ifelse(Cancer.y=="gastric cancer", "stomach", "skin")) %>%
  # head(1) %>%
  dplyr::mutate(GSVA = purrr::map2(tissue, data_spread, fn_GSVA)) %>%
  dplyr::select(-data_spread) -> GSVA.score.cancer_specific

GSVA.score.cancer_specific %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score_all-possible-features_all-cancer_specific.rds.gz"), compress = "gz")

# cancer, study, blockage,  pre/on-treatment specific ----
exp_data %>%
  tidyr::gather(-gene_id, key="Run", value = exp) %>%
  dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
  dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
  dplyr::select(-Cancer_type) %>%
  tidyr::nest(-Cancer.y,-Biopsy_Time,-blockade, -Author) %>%
  dplyr::mutate(data_spread = purrr::map(data,.f = function(.x){
    .x %>%
      tidyr::spread(key="Run",value="exp")
  })) %>%
  dplyr::select(-data) -> exp_data.nest.all_specific

exp_data.nest.all_specific[-4,] %>%
  dplyr::mutate(tissue = ifelse(Cancer.y=="gastric cancer", "stomach", "skin")) %>%
  # head(1) %>%
  dplyr::mutate(GSVA = purrr::map2(tissue, data_spread, fn_GSVA)) %>%
  dplyr::select(-data_spread) -> GSVA.score.all_specific

GSVA.score.all_specific %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score_all-possible-features_all-specific.rds.gz"), compress = "gz")

# cancer, study, blockage,  pre/on-treatment specific ----
exp_data %>%
  dplyr::mutate(tissue = "All") %>%
  tidyr::nest(-tissue,.key="data_spread") %>%
  dplyr::mutate(GSVA = purrr::map2(tissue, data_spread, fn_GSVA)) %>%
  dplyr::select(-data_spread) -> GSVA.score.all

GSVA.score.all %>%
  readr::write_rds(file.path(res_path, "ICP_GSVA_score_all-possible-features_all-togather.rds.gz"), compress = "gz")


save.image(file.path(res_path,"get_clinical_gsvascore.rds"))
