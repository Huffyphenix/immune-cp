# other expression predictors value in clincal data
library(magrittr)


# path config -------------------------------------------------------------
basic_path <- file.path("/home/huff/project")
immune_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
data_path <- file.path(basic_path,"immune_checkpoint/clinical_response_data")


exp_data<- readr::read_rds(file.path(data_path,"clinical_exp_response_combine.rds.gz"))


# cytolytic activity ------------------------------------------------------
# Ref: Molecular and genetic properties of tumors associated with local immune cytolytic activity
# (log-average (geometric mean) of GZMA and PRF1 expression in transcripts per million (TPM)) 
fn_CYT <- function(exp){
  exp %>%
    dplyr::filter(SYMBOL %in% c("GZMA","PRF1")) %>%
    dplyr::select(-ENSEMBL) %>%
    tidyr::spread(key="SYMBOL",value="exp") %>%
    dplyr::mutate(CYT_score = log(sqrt(GZMA*PRF1))) %>%
    dplyr::select(Run,CYT_score)
}

exp_data %>%
  dplyr::mutate(CYT_score = purrr::map(exp, fn_CYT)) %>%
  dplyr::select(-exp) -> clinical_CYT

# Expanded immune ------------------------------------------------------
# Ref: IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade
# After performance of quantile normalization, a log10 transformation was applied, and signature scores were calculated by
# averaging of the included genes for the IFN-γ (6-gene) and expanded immune (18-gene) signatures.
library(preprocessCore)

fn_expanded_immune <- function(exp){
  exp %>%
    dplyr::mutate(exp = ifelse(is.na(exp),0.001,exp)) %>%
    dplyr::mutate(exp = ifelse(exp==0,0.001,exp)) %>%
    dplyr::filter(SYMBOL %in% c("CD3D","IDO1","CIITA","CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13",
                                "IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")) %>%
    dplyr::select(-ENSEMBL) %>%
    tidyr::spread(key="Run",value="exp")-> exp.filter
  
  # normalize.quantiles(as.matrix(as.matrix(exp.filter[,-1]))) %>%
  #   as.data.frame() %>%
  #   dplyr::as.tbl() %>%
  #   dplyr::mutate(SYMBOL=exp.filter$SYMBOL) -> exp.filter.QN
  # colnames(exp.filter.QN) <- c(colnames(exp.filter[,-1]),"SYMBOL")
  exp.filter %>%
    tidyr::gather(-SYMBOL,key="Run",value="exp") %>%
    dplyr::mutate(exp=log10(as.numeric(exp))) %>%
    dplyr::group_by(Run) %>%
    dplyr::mutate(expanded_immune_score = mean(exp)) %>%
    dplyr::select(Run,expanded_immune_score) %>%
    unique()
}

exp_data %>%
  dplyr::mutate(expanded_immune_score = purrr::map(exp, fn_expanded_immune)) %>%
  dplyr::select(-exp) -> clinical_expanded_immune.QN

# IFN-γ  ------------------------------------------------------
# Ref: IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade
# After performance of quantile normalization, a log10 transformation was applied, and signature scores were calculated by
# averaging of the included genes for the IFN-γ (6-gene) and expanded immune (18-gene) signatures.
library(preprocessCore)

fn_IFN <- function(exp){
  exp %>%
    dplyr::mutate(exp = ifelse(is.na(exp),0.001,exp)) %>%
    dplyr::mutate(exp = ifelse(exp==0,0.001,exp)) %>%
    dplyr::filter(SYMBOL %in% c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")) %>%
    dplyr::select(-ENSEMBL) %>%
    tidyr::spread(key="Run",value="exp")-> exp.filter
  
  # normalize.quantiles(as.matrix(as.matrix(exp.filter[,-1]))) %>%
  #   as.data.frame() %>%
  #   dplyr::as.tbl() %>%
  #   dplyr::mutate(SYMBOL=exp.filter$SYMBOL) -> exp.filter.QN
  # colnames(exp.filter.QN) <- c(colnames(exp.filter[,-1]),"SYMBOL")
  exp.filter %>%
    tidyr::gather(-SYMBOL,key="Run",value="exp") %>%
    dplyr::mutate(exp=log10(as.numeric(exp))) %>%
    dplyr::group_by(Run) %>%
    dplyr::mutate(expanded_immune_score = mean(exp)) %>%
    dplyr::select(Run,expanded_immune_score) %>%
    unique()
  # exp %>%
  #   dplyr::mutate(exp = ifelse(is.na(exp),0.001,exp)) %>%
  #   dplyr::mutate(exp = ifelse(exp==0,0.001,exp)) %>%
  #   dplyr::filter(SYMBOL %in% c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")) %>%
  #   dplyr::mutate(exp=log10(as.numeric(exp))) %>%
  #   dplyr::group_by(Run) %>%
  #   dplyr::mutate(IFN_score = mean(exp)) %>%
  #   dplyr::select(Run,IFN_score) %>%
  #   unique()
}

exp_data %>%
  dplyr::mutate(IFN_score = purrr::map(exp, fn_IFN)) %>%
  dplyr::select(-exp) -> clinical_IFN_score.QN

# IMPRESS  ------------------------------------------------------
# Ref: Robust prediction of response to immune checkpoint blockade therapy in metastatic melanoma
# After performance of quantile normalization of the 28 immune checkpoint genes, a log10 transformation was applied, 
# and signature scores were calculated by f(i,j) = ifelse(i>j,1,0)
# threshold=8
library(preprocessCore)

fn_IMPRESS <- function(exp){
  exp %>%
    dplyr::filter(SYMBOL %in% c("BTLA","VSIR","CD200","CD200R1",
                                "CD27","CD274","CD276","CD28","CD40",
                                "CD80","CD86","CEACAM1","CTLA4",
                                "HAVCR2","IDO1","IL2RB","LAG3","TNFRSF18",
                                "PDCD1","PDCD1LG2","PVR","NECTIN2",
                                "TIGIT","TNFRSF14","TNFSF4","TNFRSF4","TNFRSF9",
                                "OX40L","TNFSF9")) %>%
    dplyr::select(-ENSEMBL) %>%
    tidyr::nest(-SYMBOL) %>%
    dplyr::mutate(exp_normalize = purrr::map(data,.f=function(.x){
      .x %>%
        dplyr::mutate(na_x = ifelse(is.na(exp),1,0)) %>%
        .$na_x %>%
        sum() -> na_one
      if(na_one >=1){
        .x %>%
          dplyr::mutate(exp = 0)
      }else{
        .x
      }
    })) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    tidyr::spread(key="Run",value="exp")-> exp.filter
  
  # normalize.quantiles(as.matrix(as.matrix(exp.filter[,-1]))) %>%
  #   as.data.frame() %>%
  #   dplyr::as.tbl() %>%
  #   dplyr::mutate(SYMBOL=exp.filter$SYMBOL) -> exp.filter.QN
  # colnames(exp.filter.QN) <- c(colnames(exp.filter[,-1]),"SYMBOL")
  exp.filter %>%
    tidyr::gather(-SYMBOL,key="Run",value="exp") %>%
    tidyr::nest(-Run) %>%
    dplyr::mutate(IMPRESS = purrr::map(data,fn_IMPRESS_score)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() 
}
fn_IMPRESS_score <- function(data){
  gene_pairs <- list(pair1 = c("CD274","VSIR"),
                     pair2 = c("CD28","CD276"),
                     pair3 = c("CD86","TNFSF4"),
                     pair4 = c("CD86","CD200"),
                     pair5 = c("CTLA4","TNFSF4"),
                     pair6 = c("PDCD1","TNFSF4"),
                     pair7 = c("CD80","TNFSF9"),
                     pair8 = c("CD86","HAVCR2"),
                     pair9 = c("CD28","CD86"),
                     pair10 = c("CD27","PDCD1"),
                     pair11 = c("CD40","CD274"),
                     pair12 = c("CD40","CD80"),
                     pair13 = c("CD40","CD28"),
                     pair14 = c("CD40","PDCD1"),
                     pair15 = c("TNFRSF14","CD86"))
  res <- tibble::tibble()
  for (i in names(gene_pairs)) {
    tmp <- c(gene_pairs[[i]])
    data %>%
      dplyr::filter(SYMBOL %in% tmp) %>%
      tidyr::spread(key = "SYMBOL", value ="exp") %>%
      dplyr::rename("g1"=tmp[1],"g2"=tmp[2]) %>%
      dplyr::mutate(score = ifelse(g1<g2,1,0)) %>%
      dplyr::select(score) -> res1
    rbind(res,res1) -> res
  }
  res %>%
    dplyr::mutate(IMPRESS = sum(score)) %>%
    dplyr::select(-score) %>%
    unique()
}

exp_data %>%
  dplyr::mutate(IMPRESS = purrr::map(exp, fn_IMPRESS)) %>%
  dplyr::select(-exp) -> clinical_IMPRESS_score.QN

# combine all data -------------------------------------------------------
clinical_IFN_score.QN %>%
  dplyr::inner_join(clinical_expanded_immune.QN %>% dplyr::select(-response), by=c("Cancer.y","blockade","Author","Drug")) %>%
  dplyr::inner_join(clinical_CYT %>% dplyr::select(-response), by=c("Cancer.y","blockade","Author","Drug"))%>%
  dplyr::inner_join(clinical_IMPRESS_score.QN %>% dplyr::select(-response), by=c("Cancer.y","blockade","Author","Drug")) -> combined_score

combined_score %>%
  readr::write_rds(file.path("/home/huff/project/immune_checkpoint/result_20171025/ICP_score.new/logistic_model_predict_Response/use_filtered_signatures_permutation_and_combination-from_GSVA_add_exp_ratio_cancerSpecific/select_best_and_compare","CYT_IFN_IMPRESS_score.rds.gz"),compress = "gz")

