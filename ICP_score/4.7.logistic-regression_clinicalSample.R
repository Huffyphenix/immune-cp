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

# home basic path
basic_path <- file.path("F:/胡斐斐/我的坚果云")

# Hust basi path
basic_path <- file.path("S:/坚果云/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical/class_metastic_type")

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
  dplyr::select(Run, Cancer.y, Cancer_type, blockade) %>%
  unique() -> Run_pubmed.id
exp_data %>%
  tidyr::gather(-symbol, key="Run", value = exp) %>%
  dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
  dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
  dplyr::filter(symbol %in% gene_list$symbol) %>%
  tidyr::nest(-Cancer.y,-Cancer_type, -blockade) -> exp_data.gather.genelist

exp_data.gather.genelist %>%
  dplyr::mutate(data_spread = purrr::map(data,.f = function(.x){
    .x %>%
      tidyr::spread(key="symbol",value="exp")
  })) %>%
  dplyr::select(-data) -> exp_data.nest.genelist

sample_info %>%
  dplyr::select(Run, Cancer.y, Cancer_type,blockade, Response) %>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no")) %>%
  tidyr::nest(-Cancer.y,-Cancer_type, -blockade, .key = "response") %>%
  dplyr::inner_join(exp_data.nest.genelist, by = c("Cancer.y", "Cancer_type","blockade")) -> data_for_logistic


############## fnctions to do logistic regression iteration ##################
fn_logistic <- function(exp, response, sample_group, test_n, train_n, iteration, times = 15 ){
  apply(exp[,-1],2,mad) %>% # mad() Compute the Median Absolute Deviation
    sort() %>%
    rev() %>%
    .[1:40] %>%
    names() -> top20_mad
  exp[,c("Run",top20_mad)] %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::inner_join(sample_group, by = "Run") %>%
    dplyr::filter(usage == "train") %>%
    dplyr::select(-Run, -usage) %>%
    dplyr::mutate(Response = as.factor(Response)) -> data.ready
  data.ready <- na.omit(data.ready)
  colnames(data.ready) <- gsub("-",".",colnames(data.ready))
  iteration_n <- iteration
  res <- apply(array(1:iteration_n,dim = iteration_n),1,FUN=fn_feature_filter, data = data.ready, times = times, test_n = test_n, train_n = train_n)
  # print("res get done")

  res
  # print("res output done")
}

fn_feature_filter <- function(x, data, times = 15, test_n = 3, train_n = 6){
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
    
  # for_score %>%
  #   dplyr::mutate(score_individual.class = ifelse(test_auc >= pos_auc, "success", "not_sure")) %>% 
  #   dplyr::mutate(score_individual.class = ifelse(test_auc <= neg_auc, "failures", score_individual.class)) %>%
  #   dplyr::group_by(features,score_individual.class) %>%
  #   dplyr::mutate(n = n()) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(features,score_individual.class,n) %>%
  #   unique() %>%
  #   tidyr::spread(key="score_individual.class", value = "n") %>%
  #   dplyr::mutate(failures = ifelse(is.na(failures), 0, failures)) %>%
  #   dplyr::mutate(not_sure = ifelse(is.na(not_sure), 0, not_sure)) %>%
  #   dplyr::mutate(success = ifelse(is.na(success), 0, success))-> for_binomtest
  
  # binom test
  # for_binomtest %>%
  #   dplyr::mutate(binom.p = purrr::pmap(.l = list(features, failures, not_sure, success),.f=fn_feature_selection)) %>%
  #   dplyr::select(features,binom.p) %>%
  #   tidyr::unnest() -> binomtest.p
  
  # binomtest.p %>%
  #   dplyr::inner_join(.score,by="features")
  .score
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

############## fnctions to do select top n features overlap ##################
fn_overlap_select <- function(.data, top_n = 10, overlpa_n = 5){
  topn <- data.frame()
  for(i in 1:length(.data)){
    # .data[[i]] %>%
    #   dplyr::mutate(topn = purrr::map(selection,function(.x){
    #     .x %>%
    #       top_n(top_n, score_all) %>%
    #       dplyr::select(features)
    #   })) %>%
    #   dplyr::select(-selection) %>%
    #   tidyr::unnest() -> topn_tmp
    .data[[i]] %>%
      tidyr::unnest() %>%
      dplyr::select(Cancer.y, Cancer_type,blockade, features, score_all) %>%
      dplyr::arrange(score_all) %>%
      dplyr::mutate(rank = nrow(.):1)-> topn_tmp 
    rbind(topn, topn_tmp) -> topn
  }
  
  # topn %>%
  #   tidyr::nest(features) %>%
  #   # dplyr::group_by(Cancer.y) %>%
  #   dplyr::mutate(topn_count = purrr::map(data,function(.x){
  #     .x %>%
  #       table() %>%
  #       as.data.frame() %>%
  #       dplyr::filter(Freq > overlpa_n) %>%
  #       dplyr::select(".") %>%
  #       dplyr::rename("features" = ".")
  #   })) %>%
  #   dplyr::select(-data)
  topn %>%
    dplyr::group_by(Cancer.y,blockade, features) %>%
    dplyr::mutate(rank_sum = sum( rank)) %>%
    dplyr::select(Cancer.y,Cancer_type,blockade,features,rank_sum) %>%
    dplyr::ungroup() %>%
    unique() %>%
    dplyr::arrange(rank_sum) %>%
    # dplyr::filter(rank_sum < top_n*10) %>%
    dplyr::select(Cancer.y, Cancer_type,blockade, features, rank_sum) %>%
    tidyr::nest(features,rank_sum) %>%
    dplyr::mutate(topn_count = purrr::map(data,function(.x){
          .x %>%
            top_n(top_n, rank_sum) %>%
            dplyr::select(features)
        })) %>%
        dplyr::select(-data) 
}
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
########################################## logistic regression iteration result -----------------------------------------------------
########################################## 1. by targets and cancers
data_for_logistic %>%
  dplyr::mutate(sample_group = purrr::map2(response,data_spread,fn_select_train_test,percent = 0.7)) -> data_for_logistic
data_for_logistic %>%
  readr::write_rds(file.path(res_path,"genelist_exp_for_logistic.rds.gz"),compress = "gz")
############## logistic regression
n=10 # do n fold corss validation to get more reliable results
logistic_feature_res_nrepeat <- list()
for(i in 1:n){
  data_for_logistic %>%
    # dplyr::mutate(test_n = ifelse(blockade == "anti–PD-1", 3, 2)) %>%
    # dplyr::mutate(train_n = ifelse(blockade == "anti–PD-1", 6, 3)) %>%
    dplyr::mutate(test_n = 3, train_n = 4) %>%
    dplyr::filter(! blockade == "anti–PD-1、CTLA-4") %>% # last two data sets of melanoma dont have enough response observation to train
    dplyr::mutate(itertion_res = purrr::pmap(list(data_spread, response, sample_group, test_n, train_n),.f=fn_logistic,iteration = 500, times = 15)) %>%
    dplyr::select(-response,-data_spread)-> logistic_feature_res_nrepeat[[i]]
}
############## score for regression iteration features
logistic_feature_res.forFeatureSelect_nrepeat <- list()
logistic_feature_Selected_nrepeat <- list()
for(i in 1:n){
  logistic_feature_res_nrepeat[[i]] %>%
    dplyr::mutate(itertion_res_datafrme = purrr::map(itertion_res,.f=function(.x){
      unlist(.x) %>%
        matrix(nrow = length(.x), byrow =T) %>%
        data.frame() %>%
        dplyr::as.tbl() %>%
        dplyr::rename("iteration" = "X1", "train_auc" = "X2", "test_auc" = "X3", "formula" = "X4")
    })) %>%
    dplyr::select(-itertion_res) -> logistic_feature_res.forFeatureSelect_nrepeat[[i]]
  
  logistic_feature_res.forFeatureSelect_nrepeat[[i]] %>%
    dplyr::mutate(selection = purrr::map(itertion_res_datafrme,fn_feature_score)) %>%
    dplyr::select(-itertion_res_datafrme, -sample_group) -> logistic_feature_Selected_nrepeat[[i]]
}

############## features selection
fn_overlap_select(logistic_feature_Selected_nrepeat,top_n = 10, overlpa_n = 5) -> final_feature.6
# 30 most variant genes, iteration = 500, times = 15, 10fold cross validation,top_n = 10, overlpa_n = 5
final_feature.1 <- final_feature
final_feature.2
final_feature.3
# 30 most variant genes, iteration = 1000, times = 20, 10fold cross validation,,top_n = 10, overlpa_n = 5
final_feature.4 %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"by_cancer_targets_30test.4","final_feature.4.tsv"))

# 30 most variant genes, iteration = 500, times = 15, 10fold cross validation, sum of the rank, top_n sum
final_feature.5 %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"by_cancer_targets_30test.5","final_feature.5.tsv"))
final_feature.6 %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"by_cancer_targets_30test.6","final_feature.6.tsv"))

##################### repeat 100 times get the best features
for(j in 54:100){
  n=10
  logistic_feature_res_nrepeat <- list()
  for(i in 1:n){
    data_for_logistic %>%
      dplyr::mutate(test_n = ifelse(Cancer.y == "gastric cancer", 3, 2)) %>%
      dplyr::mutate(train_n = ifelse(Cancer.y == "gastric cancer", 6, 3)) %>%
      # dplyr::mutate(test_n = 3, train_n = 4) %>%
      dplyr::filter(! blockade == "anti–PD-1、CTLA-4") %>% # last two data sets of melanoma dont have enough response observation to train
      dplyr::mutate(itertion_res = purrr::pmap(list(data_spread, response, sample_group, test_n, train_n),.f=fn_logistic,iteration = 500, times = 15)) %>%
      dplyr::select(-response,-data_spread)-> logistic_feature_res_nrepeat[[i]]
  }
  ############## score for regression iteration features
  logistic_feature_res.forFeatureSelect_nrepeat <- list()
  logistic_feature_Selected_nrepeat <- list()
  for(i in 1:n){
    logistic_feature_res_nrepeat[[i]] %>%
      dplyr::mutate(itertion_res_datafrme = purrr::map(itertion_res,.f=function(.x){
        unlist(.x) %>%
          matrix(nrow = length(.x), byrow =T) %>%
          data.frame() %>%
          dplyr::as.tbl() %>%
          dplyr::rename("iteration" = "X1", "train_auc" = "X2", "test_auc" = "X3", "formula" = "X4")
      })) %>%
      dplyr::select(-itertion_res) -> logistic_feature_res.forFeatureSelect_nrepeat[[i]]
    
    logistic_feature_res.forFeatureSelect_nrepeat[[i]] %>%
      dplyr::mutate(selection = purrr::map(itertion_res_datafrme,fn_feature_score)) %>%
      dplyr::select(-itertion_res_datafrme, -sample_group) -> logistic_feature_Selected_nrepeat[[i]]
  }
  ############## features selection
  fn_overlap_select(logistic_feature_Selected_nrepeat,top_n = 10, overlpa_n = 5) -> final_feature
  
  final_feature %>%
    tidyr::unnest() %>%
    readr::write_tsv(file.path(res_path,"repeat_100times_feature",paste(j,"final_feature.tsv",sep="_")))
}
########################################################## 2. only by cancers

############## logistic regression
# exp_data %>%
#   tidyr::gather(-symbol, key="Run", value = exp) %>%
#   dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
#   dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
#   dplyr::filter(symbol %in% gene_list$symbol) %>%
#   dplyr::select(-blockade) %>%
#   tidyr::nest(-Cancer.y) -> exp_data.gather.genelist.bycancer
# 
# exp_data.gather.genelist.bycancer %>%
#   dplyr::mutate(data_spread = purrr::map(data,.f = function(.x){
#     .x %>%
#       tidyr::spread(key="symbol",value="exp")
#   })) %>%
#   dplyr::select(-data) -> exp_data.nest.genelist.bycancer
# 
# sample_info %>%
#   dplyr::select(Run, Cancer.y, Response) %>%
#   dplyr::filter(! Response %in% c("NE", "X")) %>%
#   dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no")) %>%
#   tidyr::nest(-Cancer.y, .key = "response") %>%
#   dplyr::inner_join(exp_data.nest.genelist.bycancer, by = c("Cancer.y")) -> data_for_logistic.bycancer
# 
# data_for_logistic.bycancer  %>%
#   dplyr::mutate(sample_group = purrr::map2(response,data_spread,fn_select_train_test,percent = 0.7)) -> data_for_logistic.bycancer
# 
# logistic_feature_res.bycancer_nrepeat <- list()
# for(i in 1:n){
#   data_for_logistic.bycancer %>%
#     dplyr::mutate(itertion_res = purrr::pmap(list(data_spread, response, sample_group), .f=fn_logistic,iteration = 500)) %>%
#     dplyr::select(-response,-data_spread)-> logistic_feature_res.bycancer_nrepeat[[i]]
# }

############## score for regression iteration features
# logistic_feature_res.forFeatureSelect_nrepeat.bycancer <- list()
# logistic_feature_Selected_nrepeat.bycancer <- list()
# for(i in 1:n){
#   logistic_feature_res.bycancer_nrepeat[[i]] %>%
#     dplyr::mutate(itertion_res_datafrme = purrr::map(itertion_res,.f=function(.x){
#       unlist(.x) %>%
#         matrix(nrow = length(.x), byrow =T) %>%
#         data.frame() %>%
#         dplyr::as.tbl() %>%
#         dplyr::rename("iteration" = "X1", "train_auc" = "X2", "test_auc" = "X3", "formula" = "X4")
#     })) %>%
#     dplyr::select(-itertion_res, -sample_group) -> logistic_feature_res.forFeatureSelect_nrepeat.bycancer[[i]]
#   
#   logistic_feature_res.forFeatureSelect_nrepeat.bycancer[[i]] %>%
#     dplyr::mutate(selection = purrr::map(itertion_res_datafrme,fn_feature_score)) %>%
#     dplyr::select(-itertion_res_datafrme) -> logistic_feature_Selected_nrepeat.bycancer[[i]]
# }

############## features selection
# fn_overlap_select(logistic_feature_Selected_nrepeat.bycancer) -> final_feature.bycancer

# save the results
save.image(file = file.path(res_path,"logistic_feature_selection.rdata"))
load(file = file.path(res_path,"logistic_feature_selection.rdata"))
##################################### performing of the feature selected #####################################
fn_logistc_draw <- function(Cancer.y,Cancer_type,blockade.x,response,data_spread,sample_group,blockade,topn_count,dir,stepAIC=T){
  # if(blockade.x != "anti–PD-1"){
  #   sample_group %>%
  #     dplyr::mutate(usage = "test") -> sample_group
  # }
  ## data prepare
  data_spread %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::inner_join(sample_group, by = "Run") %>%
    dplyr::filter(usage == "test") %>%
    dplyr::select(-Run, -usage) %>%
    dplyr::mutate(Response = as.factor(Response)) -> data.ready
  data.ready <- na.omit(data.ready)
  colnames(data.ready) <- gsub("-",".",colnames(data.ready))
  
  # fomula prepare
  for(i in 1:nrow(topn_count)){
    if(i == 1){
      formula <- paste("Response", topn_count$features[i], sep = "~")
    } else{
      formula <- paste(formula, topn_count$features[i], sep = "+")
    }
  }
  # for(i in 1:nrow(topn_count)){
  #   if(i == 1){
  #     formula <- paste("Response", topn_count$.[i], sep = "~")
  #   } else{
  #     formula <- paste(formula, topn_count$.[i], sep = "+")
  #   }
  # }
  ## do logistic regression
  model <- glm( as.formula(formula), data = data.ready, family = binomial)
  
  ## stepwise regression model
  if(stepAIC){
    step.model <- stepAIC(model, direction = "backward",trace = F)
    model <- step.model
  } else{
    model <-model
  }
  
  
  # Make predictions
  probabilities <- model %>% predict(type = "response")
  predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
  observed.classes <- data.ready$Response
  # Model accuracy
  accuracy <- mean(predicted.classes == data.ready$Response)
  error <- mean(predicted.classes != data.ready$Response)
  
  # table(observed.classes, predicted.classes)
  # confusionMatrix(as.factor(predicted.classes), observed.classes,
  #                 positive = "yes")
  ## ROC
  res.roc <- roc(observed.classes, probabilities)
  png(file.path(res_path,dir,paste(Cancer.y,blockade.x,"by-feature-from",Cancer_type,blockade,"ROC.png", sep="_")),width = 400, height = 300)
  plot.roc(res.roc, print.auc = TRUE, print.thres = "best")
  dev.off()
  
  formula <- model$formula
  fn_formula_split(Reduce(paste, deparse(formula))) %>%
    readr::write_tsv(file.path(res_path,dir,paste(Cancer.y,blockade.x,"by-feature-from",blockade,"features.tsv", sep="_")))
}
########################################## 1. features selected by targets and cancers
############## used on data classified only by cancer
data_for_logistic.bycancer %>%
  dplyr::mutate(sample_group = purrr::map2(Cancer.y,sample_group, .f=function(.x, .y){
    if(.x == "melanoma"){
      .y %>%
        dplyr::mutate(usage = "test") -> .y
    } else{
      .y <- .y
    }
    .y
  })) %>%
  dplyr::inner_join(final_feature.6, by = c("Cancer.y")) %>%
  dplyr::mutate(blockade.x = "ALL30") %>%
  # dplyr::select(-test_n,-train_n) %>%
  purrr::pwalk(.f=fn_logistc_draw,dir="by_cancer_targets_30test.6",stepAIC = F) 

############## used on data classified by cancer and targets
data_for_logistic %>%
  dplyr::left_join(final_feature.6, by = c("Cancer.y","Cancer_type")) %>%
  dplyr::mutate(sample_group = purrr::pmap(list(blockade.x,sample_group,blockade.y), .f=function(.x, .y, .z){
    if(.x == "anti–PD-1、CTLA-4"){
      .y %>%
        dplyr::mutate(usage = "test") -> .y
    } else if (.x!=.z){  # predict all melanoma anti-CTLA4 by using features from  melanoma anti-PD1
      .y %>%
        dplyr::mutate(usage = "test") -> .y
    }else{
      .y <- .y
    }
    .y
  })) %>%
  dplyr::rename("blockade" = "blockade.y") %>%
  # dplyr::select(-test_n,-train_n) %>%
  purrr::pwalk(.f=fn_logistc_draw,dir="by_cancer_targets_30test.6",stepAIC = F)

########################################## 2. features selected by only cancers
############## used on data classified only by cancer
# data_for_logistic.bycancer %>%
#   dplyr::inner_join(final_feature.bycancer, by = c("Cancer.y")) %>%
#   dplyr::mutate(blockade.x = "ALL") %>%
#   purrr::pwalk(.f=fn_logistc_draw,dir="by_only_cancer_20test")
# 
# ############## used on data classified by cancer and targets
# data_for_logistic %>%
#   dplyr::mutate(sample_group = purrr::map2(Cancer.y,sample_group, .f=function(.x, .y){
#     if(.x == "melanoma"){
#       .y %>%
#         dplyr::mutate(usage = "test") -> .y
#     } else{
#       .y <- .y
#     }
#     .y
#   })) %>%
#   dplyr::inner_join(final_feature.bycancer, by = c("Cancer.y")) %>%
#   dplyr::rename("blockade.x" = "blockade") %>%
#   purrr::pwalk(.f=fn_logistc_draw,dir="by_only_cancer_20test")

########################################## 2. combine features selected from melenoma anti-pd1 and anti-ctla4 to predict all melanoma samples
final_feature.2 %>%
  tidyr::unnest() %>%
  dplyr::select(Cancer.y, ".") %>%
  unique() %>%
  dplyr::filter(Cancer.y=="melanoma")-> final_feature.2.combine.melanoma

# prepare fomula
for(i in 1:nrow(final_feature.2.combine.melanoma)){
  if(i == 1){
    formula <- paste("Response", final_feature.2.combine.melanoma$.[i], sep = "~")
  } else{
    formula <- paste(formula, final_feature.2.combine.melanoma$.[i], sep = "+")
  }
}

# prepare data
data_for_logistic.bycancer$response[[2]] %>%
  dplyr::inner_join(data_for_logistic.bycancer$data_spread[[2]], by = "Run") %>%
  dplyr::select(-Run) %>%
  dplyr::mutate(Response = as.factor(Response)) -> data.ready
colnames(data.ready) <- gsub("-",".",colnames(data.ready))
## do logistic regression
model <- glm( as.formula(formula), data = data.ready, family = binomial)


## stepwise regression model
library(MASS)
step.model <- stepAIC(model, direction = "backward")
model <- step.model

## regsubsets
library(leaps)
regsubset.model <- regsubsets(as.formula(formula), data = data.ready,method="backward",nvmax = 5)

# Make predictions
probabilities <- model %>% predict(type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
observed.classes <- data.ready$Response
# Model accuracy
accuracy <- mean(predicted.classes == data.ready$Response)
error <- mean(predicted.classes != data.ready$Response)

# table(observed.classes, predicted.classes)
# confusionMatrix(as.factor(predicted.classes), observed.classes,
#                 positive = "yes")
## ROC
res.roc <- roc(observed.classes, probabilities)
png(file.path(res_path,dir,paste(Cancer.y,blockade.x,"by-feature-from",blockade,"ROC.png", sep="_")),width = 400, height = 300)
plot.roc(res.roc, print.auc = TRUE, print.thres = "best")
dev.off()