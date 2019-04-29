##################################
# logistic regression
# clinical data of ICB blockage
# ICP gene expression data as feature
#################################
library(magrittr)
library(tidyverse)
library(caret)
library(glmnet)
library(broom)
library(pROC)

theme_set(theme_classic())

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
# gene id transfer --------------------------------------------------------
ensembl_gene <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","Homo_sapiens.gene_info.geneid.symbol.ensembl"))

gene_list %>%
  dplyr::inner_join(ensembl_gene, by = "GeneID") %>%
  dplyr::filter(!is.na(ens_id)) -> gene_list

exp_data %>%
  dplyr::rename("ens_id" = "gene_id") %>%
  dplyr::inner_join(ensembl_gene, by = "ens_id") %>%
  dplyr::select(-GeneID, -ens_id) %>%
  dplyr::rename("symbol" = "Symbol") -> exp_data
# source("/home/huff/project/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

# exp data classfication --------------------------------------------------
sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq") %>%
  dplyr::filter(Biopsy_Time == "pre-treatment")


sample_info %>%
  dplyr::select(Run, Cancer.y, blockade) %>%
  unique() -> Run_pubmed.id
exp_data %>%
  tidyr::gather(-symbol, key="Run", value = exp) %>%
  dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
  dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  tidyr::nest(-Cancer.y, -blockade) -> exp_data.gather.genelist

exp_data.gather.genelist %>%
  dplyr::mutate(data_spread = purrr::map(data,.f = function(.x){
    .x %>%
      tidyr::spread(key="symbol",value="exp")
  })) %>%
  dplyr::select(-data) -> exp_data.nest.genelist

sample_info %>%
  dplyr::select(Run, Cancer.y, blockade, Response) %>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no")) %>%
  tidyr::nest(-Cancer.y, -blockade, .key = "response") %>%
  dplyr::inner_join(exp_data.nest.genelist, by = c("Cancer.y", "blockade")) -> data_for_logistic

# logistic regression iteration result -----------------------------------------------------
# 1. by targets and cancers
data_for_logistic %>%
  head(2) %>% # last two data sets of melanoma dont have enough response observation to train
  dplyr::mutate(itertion_res = purrr::map2(data_spread, response, .f=fn_logistic,iteration = 500)) %>%
  dplyr::select(-response,-data_spread)-> logistic_feature_res

# 2. only by cancers
exp_data %>%
  tidyr::gather(-symbol, key="Run", value = exp) %>%
  dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
  dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  dplyr::select(-blockade) %>%
  tidyr::nest(-Cancer.y) -> exp_data.gather.genelist.bycancer

exp_data.gather.genelist.bycancer %>%
  dplyr::mutate(data_spread = purrr::map(data,.f = function(.x){
    .x %>%
      tidyr::spread(key="symbol",value="exp")
  })) %>%
  dplyr::select(-data) -> exp_data.nest.genelist.bycancer

sample_info %>%
  dplyr::select(Run, Cancer.y, Response) %>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no")) %>%
  tidyr::nest(-Cancer.y, .key = "response") %>%
  dplyr::inner_join(exp_data.nest.genelist.bycancer, by = c("Cancer.y")) -> data_for_logistic.bycancer

data_for_logistic.bycancer %>%
  dplyr::mutate(itertion_res = purrr::map2(data_spread, response, .f=fn_logistic,iteration = 500)) %>%
  dplyr::select(-response,-data_spread)-> logistic_feature_res.bycancer

############## fnctions to do logistic regression iteration ##################
fn_logistic <- function(exp, response, iteration){
  apply(exp[,-1],2,mad) %>% # mad() Compute the Median Absolute Deviation
    sort() %>%
    rev() %>%
    .[1:30] %>%
    names() -> top20_mad
  exp[,c("Run",top20_mad)] %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::select(-Run) %>%
    dplyr::mutate(Response = as.factor(Response)) -> data.ready
  data.ready <- na.omit(data.ready)
  colnames(data.ready) <- gsub("-",".",colnames(data.ready))
  iteration_n <- iteration
  res <- apply(array(1:iteration_n,dim = iteration_n),1,FUN=fn_feature_filter, data = data.ready, times = 15, test_n = 3, train_n = 10)
  # print("res get done")

  res
  # print("res output done")
}

fn_feature_filter <- function(x, data, times = 15, test_n = 3, train_n = 10){
  # data prepare
  test_set <- sample_n(data %>% dplyr::filter(Response == "yes"),test_n) %>%
    rbind(sample_n(data %>% dplyr::filter(Response == "no"),test_n))
  train_set <- sample_n(data %>% dplyr::filter(Response == "yes"),train_n) %>%
    rbind(sample_n(data %>% dplyr::filter(Response == "no"),train_n))
  feature <- tibble::tibble(feature = colnames(train_set)[-ncol(train_set)] )
  
  for(i in 1:times){
    if(i == 1){
      f <- sample_n(feature,1) %>%
        .$feature
      feature <- feature %>% dplyr::filter(! feature %in% f)
      formula1 <- as.formula(paste("Response ~ ", f))
      auc <- 0
    } else{
      f1 <- sample_n(feature,1) %>%
        .$feature
      feature <- feature %>% dplyr::filter(! feature %in% f1)
      formula1 <- as.formula(paste(Reduce(paste, deparse(formula)), f1, sep = "+"))
    }
    print(Reduce(paste, deparse(formula1)))
    auc1 <- fn_auc(test_set = train_set, train_set = train_set, formula = formula1)
    if(i < times){
      if(auc1 < 1){
        if(auc1 > auc){
          auc <- auc1
          formula <- formula1
          }else{
            formula <- formula
            auc <- auc
            }
        }else if(auc1 ==1){
          formula <- formula1
          auc <- auc
          auc_test <- fn_auc(test_set = test_set, train_set = train_set, formula = formula)
          break()
        }
    } else{
      if(auc1 <= 1){
        if(auc1 > auc){
          auc <- auc1
          formula <- formula1
        }else{
          formula <- formula
          auc <- auc
        }
      }
      auc_test <- fn_auc(test_set = test_set, train_set = train_set, formula = formula)
    }
  }
  
  tibble::tibble(iteration = x, auc_train = auc, auc_test = auc_test, formula = Reduce(paste, deparse(formula)))
  # print("iteration done")
}

fn_auc <- function(test_set, train_set, formula){
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
  ## best threshold
  
  # plot.roc(res.roc, print.auc = TRUE)
  auc
}


# seletion of features from iteration result ------------------------------
logistic_feature_res %>%
  dplyr::mutate(itertion_res_datafrme = purrr::map(itertion_res,.f=function(.x){
    unlist(.x) %>%
      matrix(nrow = length(.x), byrow =T) %>%
      data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::rename("iteration" = "X1", "train_auc" = "X2", "test_auc" = "X3", "formula" = "X4")
  })) %>%
  dplyr::select(-itertion_res) -> logistic_feature_res.forFeatureSelect

logistic_feature_res.bycancer %>%
  dplyr::mutate(itertion_res_datafrme = purrr::map(itertion_res,.f=function(.x){
    unlist(.x) %>%
      matrix(nrow = length(.x), byrow =T) %>%
      data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::rename("iteration" = "X1", "train_auc" = "X2", "test_auc" = "X3", "formula" = "X4")
  })) %>%
  dplyr::select(-itertion_res) -> logistic_feature_res.bycancer.forFeatureSelect

############## fnctions to do feature score and selection ##################

#### calculate the feature score ####
fn_feature_score <- function(.x, pos_auc = 0.6, neg_auc = 0.4){
  .x %>%
    dplyr::mutate(features = purrr::map(formula, fn_formula_split)) %>%             # merge features into a column
    dplyr::select(-formula) %>%
    tidyr::unnest() %>%
    dplyr::mutate(test_auc = as.numeric(as.character(test_auc))) -> for_score
  
  for_score %>%
    dplyr::mutate(score_individual = ifelse(test_auc >= pos_auc, 1, 0)) %>%            # successful iterations (iterations with test-AUC >=0.6)
    dplyr::mutate(score_individual = ifelse(test_auc <= neg_auc, -1, score_individual)) %>%            # unsuccessful iterations (iterations with test-AUC <=0.4)
    dplyr::group_by(features) %>%
    dplyr::mutate(score_all = sum(score_individual)) %>%
    dplyr::arrange(desc(score_all)) %>%
    dplyr::ungroup() %>%
    dplyr::select(features, score_all) %>%
    unique() -> .score
    
  for_score %>%
    dplyr::mutate(score_individual.class = ifelse(test_auc >= pos_auc, "success", "not_sure")) %>% 
    dplyr::mutate(score_individual.class = ifelse(test_auc <= neg_auc, "failures", score_individual.class)) %>%
    dplyr::group_by(features,score_individual.class) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(features,score_individual.class,n) %>%
    unique() %>%
    tidyr::spread(key="score_individual.class", value = "n") -> for_binomtest
  
  # binom test
  for_binomtest %>%
    dplyr::mutate(fail_notsure = failures+not_sure) %>%
    dplyr::mutate(binom.p = purrr::pmap(fn_feature_selection))
  
}
#### split the formula string and extract the features ####
fn_formula_split <- function(formula){
  as.character(formula) %>%
    gsub(pattern = " ",replacement = "") %>%
    strsplit(split = "\\+") %>%
    .[[1]] %>%
    strsplit(split = "\\~") %>%
    unlist() %>%
    .[-1] -> .features
  tibble::tibble(features = .features)
}

#### binom test to test pvalue ####
fn_feature_selection <- function(features, failures, not_sure, success){
  print(features)
  binom.test(c(success,failures))$p.value
}
