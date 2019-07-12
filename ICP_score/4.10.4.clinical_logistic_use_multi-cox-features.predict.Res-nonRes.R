############ Features filtered by multi-cox analysis from TCGA data
############ add thses features to logistic model to predict clinical sample response to immunotherapy
library(magrittr)
library(tidyverse)
library(caret)
library(glmnet)
library(broom)
library(pROC)
library(MASS)

# data path ---------------------------------------------------------------
# server 1 basic path
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
score_path <- file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")
res_path <- file.path(score_path,"logistic_model_predict_Response/use_only_features_fit_multi-cox-model")

# load data ---------------------------------------------------------------
# GSVA score
gsva.score <- readr::read_rds(file.path(score_path,"ICP_GSVA_score-by-matastatic-or-not_all-possible-features_class-cancer_specific.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(-Cancer.y, -tissue)

# feature filtered by tcga multi-cox analysis
multi_cox_PFS <- readr::read_tsv(file.path(immune_res_path,"TCGA_GSVAScore/cancer_specific/3.survival_with_GSVA_score.new","GSVA.score.multi-varite.surv.PFS.tsv")) %>%
  dplyr::filter(cancer_types %in% c("SKCM","STAD")) %>%
  tidyr::nest(-cancer_types,.key="Features")

multi_cox_OS <- readr::read_tsv(file.path(immune_res_path,"TCGA_GSVAScore/cancer_specific/3.survival_with_GSVA_score.new","2.GSVA.score.multi-varite.surv.OS.tsv")) %>%
  dplyr::filter(cancer_types %in% c("SKCM","STAD")) %>%
  tidyr::nest(-cancer_types,.key="Features")

# exp data classfication --------------------------------------------------
sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment")


sample_info %>%
  dplyr::select(Run, Cancer.y, Cancer_type, blockade) %>%
  unique() -> Run_pubmed.id
gsva.score %>%
  dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
  # dplyr::filter(symbol %in% gene_list$symbol) %>%
  dplyr::select(-blockade,-Cancer_type) %>%
  tidyr::nest(-Cancer.y) -> exp_data.gather.genelist

exp_data.gather.genelist %>%
  dplyr::rename("data_spread" = "data") -> exp_data.nest.genelist

sample_info %>%
  dplyr::select(Run, Cancer.y, Response) %>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no")) %>%
  tidyr::nest(-Cancer.y, .key = "response") %>%
  dplyr::inner_join(exp_data.nest.genelist, by = c("Cancer.y")) -> data_for_logistic

# functions ---------------------------------------------------------------

fn_select_train_test <- function(response,data_spread,percent = 0.7){
  train_number <- response$Response %>%
    table() %>%
    as.data.frame() %>%
    dplyr::mutate(train_n = round(Freq*percent))
  train_number %>%
    dplyr::filter(. == "no") %>%
    .$train_n -> no_n
  train_number %>%
    dplyr::filter(. == "yes") %>%
    .$train_n -> yes_n
  set.seed(123)
  rbind(sample_n(response %>% dplyr::filter(Response == "no"), no_n),
        sample_n(response %>% dplyr::filter(Response == "yes"), yes_n)) %>%
    .$Run -> train_samples
  response %>%
    dplyr::filter(! Run %in% train_samples) %>%
    .$Run -> test_samples
  tibble::tibble(Run = c(train_samples,test_samples), 
                 usage = c(rep("train",length(train_samples)),rep("test",length(test_samples))))
}

fn_logistic <- function(Cancer_type,exp, response, sample_group,Features, filename,test_n=7, train_n=7, iteration = 500, times = 15 ){
  # apply(exp[,-1],2,mad) %>% # mad() Compute the Median Absolute Deviation
  #   sort() %>%
  #   rev() %>%
  #   .[1:40] %>%
  #   names() -> top20_mad
  print(Cancer_type)
  exp %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::inner_join(sample_group, by = "Run") %>%
    dplyr::filter(usage == "train") %>%
    dplyr::select(-Run, -usage) %>%
    dplyr::mutate(Response = as.factor(Response)) -> train_set
  train_set[is.na(train_set)] <- 0
  colnames(train_set) <- gsub(" ","_",colnames(train_set))
  
  exp %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::inner_join(sample_group, by = "Run") %>%
    dplyr::filter(usage == "test") %>%
    dplyr::select(-Run, -usage) %>%
    dplyr::mutate(Response = as.factor(Response)) -> test_set
  test_set[is.na(test_set)] <- 0
  colnames(test_set) <- gsub(" ","_",colnames(test_set))
  
  Features %>%
    dplyr::filter(Features %in% colnames(train_set)) %>%
    dplyr::filter(coxp<0.1) %>%
    .$Features -> Features
  formula <- as.formula(paste("Response ~ ", paste(Features,collapse = "+")))
  
  auc_test <- fn_auc(test_set = test_set, train_set = train_set, formula = formula,filename=paste(filename,"AUC_on_test_set",sep="."))
  auc_train <- fn_auc(test_set = train_set, train_set = train_set, formula = formula,filename=paste(filename,"AUC_on_train_set",sep="."))
  
  tibble::tibble(auc_train = auc_train, auc_test = auc_test, formula = Reduce(paste, deparse(formula)))

  # iteration_n <- iteration
  # res <- apply(array(1:iteration_n,dim = iteration_n),1,FUN=fn_feature_filter, data = data.ready, times = times, test_n = test_n, train_n = train_n)
  # print("res get done")
  
  # res
  # print("res output done")
}

fn_auc <- function(test_set, train_set, formula,filename){
  # Fit the model
  model <- glm( formula, data = train_set, family = binomial)
  # Make predictions
  probabilities <- model %>% predict(test_set, type = "response")
  predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
  observed.classes <- test_set$Response
  # Model accuracy
  accuracy <- mean(predicted.classes == test_set$Response)
  error <- mean(predicted.classes != test_set$Response)
  
  # table(observed.classes, predicted.classes)
  # confusionMatrix(as.factor(predicted.classes), observed.classes,
  #                 positive = "yes")
  # ROC
  res.roc <- roc(observed.classes, probabilities) 

  auc <- res.roc$auc[[1]]
  
  ## plot
  png(file.path(res_path,paste(filename,auc,"ROC.png", sep="_")),width = 400, height = 300)
  plot.roc(res.roc, print.auc = TRUE, print.thres = "best")
  dev.off()
  
  auc
}

# select train set ----------
data_for_logistic %>%
  dplyr::mutate(sample_group = purrr::map2(response,data_spread,fn_select_train_test,percent = 0.7)) -> data_for_logistic

# calculation --------
data_for_logistic %>%
  dplyr::mutate(cancer_types = ifelse(Cancer.y=="melanoma","SKCM","STAD")) %>%
  dplyr::inner_join(multi_cox_OS,by="cancer_types") %>% 
  dplyr::mutate(filename = paste(cancer_types,"OS",sep="_")) %>%
  dplyr::mutate(AUC = purrr::pmap(list(Cancer_type=Cancer.y,exp=data_spread,response=response,sample_group=sample_group,Features=Features,filename=filename),fn_logistic)) %>%
  dplyr::select(Cancer.y,AUC) %>%
  tidyr::unnest()


data_for_logistic %>%
  dplyr::mutate(cancer_types = ifelse(Cancer.y=="melanoma","SKCM","STAD")) %>%
  dplyr::inner_join(multi_cox_PFS,by="cancer_types") %>%
  dplyr::mutate(filename = paste(cancer_types,"PFS",sep="_")) %>%
  dplyr::mutate(AUC = purrr::pmap(list(Cancer_type=Cancer.y,exp=data_spread,response=response,sample_group=sample_group,Features=Features,filename=filename),fn_logistic)) %>%
  dplyr::select(Cancer.y,AUC) %>%
  tidyr::unnest()

# GSVA * coef of cox, the same as uppon result.
# data_for_logistic %>%
#   dplyr::mutate(cancer_types = ifelse(Cancer.y=="melanoma","SKCM","STAD")) %>%
#   dplyr::inner_join(multi_cox_OS,by="cancer_types") %>% 
#   dplyr::mutate(data_spread = purrr::map2(data_spread,Features, .f=function(.x,.y){
#     .x %>%
#       tidyr::gather(-Run,key="Features",value="GSVA") %>%
#       dplyr::mutate(Features=gsub(" ","_",Features)) %>%
#       dplyr::inner_join(.y %>%
#                           dplyr::select(Features,coef),by="Features") %>%
#       dplyr::mutate(GSVA=GSVA*coef) %>%
#       dplyr::select(-coef) %>%
#       tidyr::spread(key = "Features",value="GSVA")
#   })) %>%
#   dplyr::mutate(filename = paste(cancer_types,"OS","GSVAxCoef",sep="_")) %>%
#   dplyr::mutate(AUC = purrr::pmap(list(Cancer_type=Cancer.y,exp=data_spread,response=response,sample_group=sample_group,Features=Features,filename=filename),fn_logistic)) %>%
#   dplyr::select(Cancer.y,AUC) %>%
#   tidyr::unnest()
# 
# 
# data_for_logistic %>%
#   dplyr::mutate(cancer_types = ifelse(Cancer.y=="melanoma","SKCM","STAD")) %>%
#   dplyr::inner_join(multi_cox_PFS,by="cancer_types") %>%
#   dplyr::mutate(filename = paste(cancer_types,"PFS",sep="_")) %>%
#   dplyr::mutate(AUC = purrr::pmap(list(Cancer_type=Cancer.y,exp=data_spread,response=response,sample_group=sample_group,Features=Features,filename=filename),fn_logistic)) %>%
#   dplyr::select(Cancer.y,AUC) %>%
#   tidyr::unnest()

save.image(file.path(res_path,"logistic.rda"))
