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
    dplyr::mutate(exp=log10(as.numeric(exp))) %>%
    dplyr::group_by(Run) %>%
    dplyr::mutate(expanded_immune_score = mean(exp)) %>%
    dplyr::select(Run,expanded_immune_score) %>%
    unique()
}

exp_data %>%
  dplyr::mutate(expanded_immune_score = purrr::map(exp, fn_expanded_immune)) %>%
  dplyr::select(-exp) -> clinical_expanded_immune

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
    dplyr::mutate(exp=log10(as.numeric(exp))) %>%
    dplyr::group_by(Run) %>%
    dplyr::mutate(IFN_score = mean(exp)) %>%
    dplyr::select(Run,IFN_score) %>%
    unique()
}

exp_data %>%
  dplyr::mutate(IFN_score = purrr::map(exp, fn_IFN)) %>%
  dplyr::select(-exp) -> clinical_IFN_score



# combine all data -------------------------------------------------------
clinical_IFN_score %>%
  dplyr::inner_join(clinical_expanded_immune %>% dplyr::select(-response), by=c("Cancer.y","blockade","Author","Drug")) %>%
  dplyr::inner_join(clinical_CYT %>% dplyr::select(-response), by=c("Cancer.y","blockade","Author","Drug")) -> combined_score

combined_score %>%
  readr::write_rds(file.path("/home/huff/project/immune_checkpoint/result_20171025/ICP_score.new/logistic_model_predict_Response/use_filtered_signatures_permutation_and_combination-from_GSVA_add_exp_ratio_cancerSpecific/select_best_and_compare","CYT_IFN_score.rds.gz"),compress = "gz")

