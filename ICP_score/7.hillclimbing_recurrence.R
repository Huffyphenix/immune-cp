
# reproducibility of hill-climbing filter features ------------------------

################################
library(magrittr)
library(tidyverse)
library(caret)
library(glmnet)
library(broom)
library(pROC)
library(MASS)
library(ggplot2)
library(multidplyr)
library(parallel)

# data path ---------------------------------------------------------------
# server 1 basic path
basic_path <- file.path("/home/huff/project")
# hust basic path
# basic_path <- file.path("S:/坚果云/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
exp_score_path <- file.path(immune_res_path,"ICP_score.new")
gsva_score_path <-  file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")
res_path <- file.path(exp_score_path,"logistic_model_predict_Response/hill_climbing_191223")

# load data ---------------------------------------------------------------
readr::read_rds(file.path(exp_score_path,"clinical_data_for_logistic.rds.gz")) -> data_for_logistic
filtered_features <- readr::read_tsv(file.path("/home/huff/project/immune_checkpoint/result_20171025/TCGA_GSVAScore/GSVA_add_exp_ratio_cancerSpecific-new/ALL_Feature_ranked_by_FC.tsv")) %>%
  dplyr::filter(top100_counts>=10) %>%
  .$Features
data_for_logistic %>%
  dplyr::mutate(GSVA = purrr::map(GSVA, .f = function(.x){
    colnames(.x) <- gsub(colnames(.x),pattern = " ",replacement = "_")
    colnames(.x) <- gsub(colnames(.x),pattern = "-",replacement = ".")
    filtered_features <- intersect(colnames(.x),filtered_features)
    .x[,c("Run",filtered_features)]
  })) -> data_for_logistic

Response_statistic <- readr::read_tsv(file.path(gsva_score_path,"Response_statistic_each_dataset.tsv")) %>%
  dplyr::select(-Cancer.y) %>%
  dplyr::mutate(Response = ifelse(is.na(Response),0,Response)) %>%
  dplyr::group_by(Author,Biopsy_Time) %>%
  dplyr::mutate(`non-Response`=sum(`non-Response`),Response=sum(Response)) %>%
  dplyr::select(-`Response_percentage(%)`) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::mutate(`Response_percentage(%)` = 100*Response/(Response+`non-Response`))

my_theme <-   theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_line(colour=NA),
  axis.text.y = element_text(size = 10,colour = "black"),
  axis.text.x = element_text(size = 10,colour = "black"),
  # legend.position = "none",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.background = element_blank(),
  legend.key = element_rect(fill = "white", colour = "black"),
  plot.title = element_text(size = 20),
  axis.text = element_text(colour = "black"),
  strip.background = element_rect(fill = "white",colour = "black"),
  strip.text = element_text(size = 10),
  text = element_text(color = "black")
)

############## fnctions to do logistic regression iteration ##################
fn_logistic <- function(Cancer_type,exp, response, test_n, train_n, iteration, CV, times = 10, fold_info){
  # apply(exp[,-1],2,mad) %>% # mad() Compute the Median Absolute Deviation
  #   sort() %>%
  #   rev() %>%
  #   .[1:40] %>%
  #   names() -> top20_mad
  # print(Cancer_type)
  exp %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::inner_join(fold_info, by = "Run") %>%
    dplyr::filter(! fold %in% CV) %>%
    dplyr::select(-fold) %>%
    dplyr::mutate(Response = as.factor(Response)) -> data.ready
  data.ready[is.na(data.ready)] <- 0
  # colnames(data.ready) <- gsub(" ",".",colnames(data.ready))
  iteration_n <- iteration
  res <- apply(array(1:iteration_n,dim = iteration_n),1,FUN=fn_feature_filter, data = data.ready, times = times, test_n = test_n, train_n = train_n, CV = CV)
  # print("res get done")
  
  res
  # print("res output done")
}

fn_feature_filter <- function(x, data, times, test_n = 3, train_n = 6, CV){
  # data prepare
  train_test_set <- readr::read_tsv(file.path(repeat_res_path,paste("I",x,"CV",CV,"train_test_set.tsv",sep=".")))
  data %>%
    dplyr::inner_join(train_test_set, by="Run") %>%
    dplyr::filter(class == "train") %>%
    dplyr::select(-Run,-class) -> train_set
  data %>%
    dplyr::inner_join(train_test_set, by="Run") %>%
    dplyr::filter(class == "test") %>%
    dplyr::select(-Run,-class) -> test_set
  
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
    # print(Reduce(paste, deparse(formula1)))
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
        auc <- auc1
        auc_test <- fn_auc(test_set = test_set, train_set = train_set, formula = formula)
        break()
      }
    } else{
      if(auc1 > auc){
        auc <- auc1
        formula <- formula1
      }else{
        formula <- formula
        auc <- auc
      }
      auc_test <- fn_auc(test_set = test_set, train_set = train_set, formula = formula)
    }
  }
  
  tibble::tibble(CV = CV, iteration = x, auc_train = auc, auc_test = auc_test, formula = Reduce(paste, deparse(formula)))
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
  res.roc <- roc(observed.classes, probabilities,quiet=TRUE)
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
    dplyr::group_by(CV,features) %>%
    dplyr::mutate(score_all = sum(score_individual)) %>%
    dplyr::arrange(desc(score_all)) %>%
    dplyr::ungroup() %>%
    dplyr::select(CV,features, score_all) %>%
    unique() -> .score
  
  for_score %>%
    dplyr::mutate(score_individual.class = ifelse(test_auc >= pos_auc, "success", "not_sure")) %>%
    dplyr::mutate(score_individual.class = ifelse(test_auc <= neg_auc, "failures", score_individual.class)) %>%
    dplyr::group_by(CV, features,score_individual.class) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(CV,features,score_individual.class,n) %>%
    unique() %>%
    tidyr::spread(key="score_individual.class", value = "n") -> for_binomtest
  if(!"failures" %in% colnames(for_binomtest)){
    for_binomtest %>%
      dplyr::mutate(failures=0) -> for_binomtest
  }else{
    for_binomtest %>%
      dplyr::mutate(failures = ifelse(is.na(failures), 0, failures)) %>%
      dplyr::mutate(not_sure = ifelse(is.na(not_sure), 0, not_sure)) %>%
      dplyr::mutate(success = ifelse(is.na(success), 0, success))-> for_binomtest
  }
  
  # binom test
  for_binomtest %>%
    dplyr::mutate(binom.p = purrr::pmap(.l = list(features, failures, success),.f=fn_feature_selection)) %>%
    dplyr::select(CV,features,binom.p) %>%
    tidyr::unnest() -> binomtest.p
  
  # success frequency
  for_binomtest %>%
    dplyr::mutate(success_freq = success/(failures+not_sure+success)) -> success_freq
  
  binomtest.p %>%
    dplyr::inner_join(.score,by=c("features","CV")) %>%
    dplyr::inner_join(success_freq,by=c("features","CV"))
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
fn_feature_selection <- function(features, failures, success){
  # print(features)
  binom.test(c(success,failures),p=0.8,alternative = c("greater"),conf.level = 0.95)$p.value
}

############## fnctions to do select top n features overlap ##################
fn_feature_select_accross_CV <- function(.data, top_n = 10, overlpa_n = 5){
  score <- tibble::tibble()
  binom <- tibble::tibble()
  success <- tibble::tibble()
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
      dplyr::select(CV, features, binom.p) %>% 
      tidyr::spread(key="CV",value="binom.p") -> binom_tmp
    
    .data[[i]] %>%
      tidyr::unnest() %>%
      dplyr::select(CV, features, success_freq) %>% 
      tidyr::spread(key="CV",value="success_freq") -> success_tmp
    
    .data[[i]] %>%
      tidyr::unnest() %>%
      dplyr::select(CV, features, score_all) %>%  
      dplyr::mutate(rank = rank(-score_all)) %>%
      dplyr::select(-score_all) %>%
      tidyr::spread(key="CV",value="rank") -> score_tmp
    if (i == 1) {
      binom_tmp -> binom
      score_tmp -> score
      success_tmp -> success
    } else{
      binom %>%
        dplyr::inner_join(binom_tmp, by="features") -> binom
      score %>%
        dplyr::inner_join(score_tmp, by="features") -> score
      success %>%
        dplyr::inner_join(success_tmp, by="features") -> success
    }
  }
  
  binom %>%
    dplyr::filter(`1`<=0.05 & `2`<=0.05 & `3` <=0.05 & `4` <=0.05 &`5`<=0.05) -> features_all_binomP_sig
  binom %>%
    dplyr::mutate(mean = (`1`+`2`+`3`+`4`+`5`)/5) %>%
    dplyr::arrange(mean) %>%
    readr::write_tsv(file.path(repeat_res_path,"recurrence", "binomP_accross_CV.tsv"))
  success %>%
    dplyr::filter(`1`>=0.6 & `2`>=0.6 & `3` >=0.6 & `4` >=0.6 &`5`>=0.6) -> features_all_success_0.6
  success %>%
    dplyr::mutate(mean = (`1`+`2`+`3`+`4`+`5`)/5) %>%
    dplyr::arrange(desc(mean)) %>%
    readr::write_tsv(file.path(repeat_res_path,"recurrence", "successFreq_accross_CV.tsv"))
  # list(features_all_binomP_sig = features_all_binomP_sig$features,
  #      features_all_success = features_all_success_0.6$features)
  binom %>%
    dplyr::mutate(binom_rank = rank((`1`+`2`+`3`+`4`+`5`)/5)) %>%
    dplyr::inner_join(success%>%
                        dplyr::mutate(success_rank = rank(-(`1`+`2`+`3`+`4`+`5`)/5)), by ="features")
    
}

fn_auc_on_validation <- function(response,data_spread, model){
  ## data prepare
  data_spread %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::select(-Run) %>%
    dplyr::mutate(Response = as.factor(Response)) -> data.ready
  # data.ready <- na.omit(data.ready)
  data.ready[is.na(data.ready)] <- 0
  colnames(data.ready) <- gsub(" ",".",colnames(data.ready))
  
  # Make predictions
  probabilities <- model %>% predict(data.ready, type = "response")
  predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
  observed.classes <- data.ready$Response
  
  # ROC
  res.roc <- roc(observed.classes, probabilities,quiet=TRUE)
  
  auc <- res.roc$auc[1]
  tibble::tibble(Sensitivity=res.roc$sensitivities,
                 Specificity=res.roc$specificities,
                 thresholds = res.roc$thresholds) %>%
    dplyr::mutate(n=1:nrow(.)) %>%
    dplyr::group_by(n) %>%
    dplyr::mutate(sum = sum(Sensitivity+Specificity)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(best = ifelse(sum == max(sum),"best","no")) %>%
    dplyr::select(Sensitivity, Specificity,  thresholds, best) %>%
    dplyr::mutate(auc=auc) %>%
    tidyr::nest(-auc,.key="roc_data")
}


# recurrence --------------------------------------------------------------

test_n=4
train_n=8
times=858

repeat_res_path <- file.path(res_path,paste("repeat",times,sep="_"))
fold_res <- 
  readr::read_tsv(file.path(repeat_res_path,"cross_fold_info.tsv"))

data_for_logistic %>%
  dplyr::mutate(n=purrr::map(response,.f=function(.x){
    nrow(.x)
  })) %>%
  tidyr::unnest(n) %>%
  dplyr::filter(n>5) %>%
  dplyr::select(-n) -> .data_for_logistic  
##################### repeat 100 times get the best features

feature_rank <- list()

for(x in 1:100){
  n=5 # n fold cross validation
  logistic_feature_res_nrepeat <- list()
  for(i in 1:n){
    .data_for_logistic %>%
      dplyr::filter(Biopsy_Time=="on-treatment" & Author=="Riaz") %>%
      dplyr::mutate(test_n = test_n, train_n = train_n) %>%
      dplyr::mutate(itertion_res = purrr::pmap(list(Cancer.y,GSVA, response, test_n, train_n),.f=fn_logistic,iteration = 500, CV = i, times = 10, fold_info=fold_res)) %>%
      dplyr::select(-response,-GSVA)-> logistic_feature_res_nrepeat[[i]]
  }
  
  ############## score for regression iteration features
  print("# get hill climbing score")
  logistic_feature_res.forFeatureSelect_nrepeat <- list()
  logistic_feature_Selected_nrepeat <- list()
  for(i in 1:n){
    logistic_feature_res_nrepeat[[i]] %>%
      dplyr::mutate(itertion_res_datafrme = purrr::map(itertion_res,.f=function(.x){
        unlist(.x) %>%
          matrix(nrow = length(.x), byrow =T) %>%
          data.frame() %>%
          dplyr::as.tbl() %>%
          dplyr::rename("CV" = "X1", "iteration" = "X2", "train_auc" = "X3", "test_auc" = "X4", "formula" = "X5")
      })) %>%
      dplyr::select(-itertion_res) -> logistic_feature_res.forFeatureSelect_nrepeat[[i]]
    
    logistic_feature_res.forFeatureSelect_nrepeat[[i]] %>%
      dplyr::mutate(selection = purrr::map(itertion_res_datafrme,fn_feature_score,pos_auc = 0.6, neg_auc = 0.4)) %>%
      dplyr::select(-itertion_res_datafrme) -> logistic_feature_Selected_nrepeat[[i]]
  }
  ############## features selection
  print("# features selection")
  fn_feature_select_accross_CV(logistic_feature_Selected_nrepeat,top_n = 10, overlpa_n = 5) -> feature_rank[[x]]
  # x = x + 1
}
feature_rank %>%
  readr::write_rds(file.path(repeat_res_path,"recurrence","feature_rank.rds.gz"),compress = "gz")

  ############## draw picture, using stepAI
  # get model from 70% train and use it on 30% Test and validation data
#   print("# perform of features selected among 5 CV test data")
#   TEST_AUC <- tibble::tibble()
#   for (i in 1:n) {
#     .data_for_logistic %>%
#       dplyr::filter(Biopsy_Time=="on-treatment" & Author=="Riaz")  %>% 
#       dplyr::mutate(data.ready = purrr::map2(GSVA,response,.f=function(.x,.y){
#         .x %>%
#           dplyr::inner_join(.y, by = "Run") %>%
#           dplyr::inner_join(fold_res, by= "Run") %>%
#           dplyr::filter(! fold %in% i) %>%
#           dplyr::select(-fold) %>%
#           dplyr::select(-Run) %>%
#           dplyr::mutate(Response = as.factor(Response))
#       })) %>%
#       dplyr::select(data.ready) %>%
#       tidyr::unnest() -> data.ready 
#     data.ready[is.na(data.ready)] <- 0
#     
#     # Make predictions
#     .data_for_logistic %>%
#       dplyr::filter(Biopsy_Time=="on-treatment" & Author=="Riaz")  %>% 
#       dplyr::mutate(GSVA = purrr::map(GSVA,.f=function(.x){
#         .x %>%
#           dplyr::inner_join(fold_res, by= "Run") %>%
#           dplyr::filter(fold %in% i) %>%
#           dplyr::select(-fold) 
#       })) %>%
#       dplyr::mutate(response = purrr::map(response,.f=function(.x){
#         .x %>%
#           dplyr::inner_join(fold_res, by= "Run") %>%
#           dplyr::filter(fold %in% i) %>%
#           dplyr::select(-fold) 
#       }))  -> data.ready.test
#     
#     # select p all significant
#     if(length(final_feature$features_all_binomP_sig)>0){
#       formula.train_binom <- paste("Response~", paste(final_feature$features_all_binomP_sig,collapse = "+"))
#       model.train_binom <- glm( as.formula(formula.train_binom), data = data.ready, family = binomial)
#       step.model.train_binom <- stepAIC(model.train_binom, direction = "backward",trace = F)
#       final_feature_binom <- rownames(as.data.frame(summary(step.model.train_binom)$coefficients))[-1]
#       
#       data.ready.test %>%
#         dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model.train_binom)) %>%
#         dplyr::select(-response,-GSVA) %>%
#         tidyr::unnest() %>%
#         dplyr::mutate(CV = i, model = "from_binom") %>%
#         dplyr::select(-roc_data) -> TEST_auc_binom
#     }else{
#       TEST_auc_binom <- tibble::tibble(Cancer.y="melanoma",blockade="anti–PD-1",Biopsy_Time="on-treatment",Author="Riaz",auc=0,CV= i,model= "from_binom")
#     }
#     
#     
#     # select frequency all > 0.6
#     if(length(final_feature$features_all_success)>0){
#       formula.train_success <- paste("Response~", paste(final_feature$features_all_success,collapse = "+"))
#       model.train_success <- glm( as.formula(formula.train_success), data = data.ready, family = binomial)
#       step.model.train_success <- stepAIC(model.train_success, direction = "backward",trace = F)
#       final_features_success <- rownames(as.data.frame(summary(step.model.train_success)$coefficients))[-1]
#       
#       data.ready.test %>%
#         dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model.train_success)) %>%
#         dplyr::select(-response,-GSVA) %>%
#         tidyr::unnest() %>%
#         dplyr::mutate(CV = i, model = "from_successFreq")%>%
#         dplyr::select(-roc_data)  -> TEST_auc_success
#     }else{
#       TEST_auc_success <- tibble::tibble(Cancer.y="melanoma",blockade="anti–PD-1",Biopsy_Time="on-treatment",Author="Riaz",auc=0,CV= i,model= "from_successFreq")
#     }
#     
#     TEST_AUC <- rbind(TEST_AUC, TEST_auc_binom, TEST_auc_success)
#   }
#   if (all(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc>=0.6)) {
#     success <- "yes"
#   } else{
#     success <- "no"
#   }
# } 
# 
# if(mean(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc) > mean(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc)){
#   model_type <- "from_binom"
#   formula.train <- paste("Response~", paste(final_feature$features_all_binomP_sig,collapse = "+"))
#   mean_test_auc <- mean(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc)
#   final_features <- final_feature$features_all_binomP_sig
#   if (all(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc>=0.6)) {
#     success_test <- "yes"
#   } else{
#     success_test <- "no"
#   }
# } else if (mean(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc) > mean(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc)){
#   model_type <- "from_successFreq"
#   formula.train <- paste("Response~", paste(final_feature$features_all_success,collapse = "+"))
#   mean_test_auc <- mean(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc)
#   final_features <- final_feature$features_all_success
#   if (all(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc>=0.6)) {
#     success_test <- "yes"
#   } else{
#     success_test <- "no"
#   }
# } else{
#   formula.train <- paste("Response~1")
#   mean_test_auc <- mean(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc)
#   final_features <- "NNa"
# }
# .data_for_logistic %>%
#   dplyr::filter(Biopsy_Time=="on-treatment" & Author=="Riaz")  %>% 
#   dplyr::mutate(data.ready = purrr::map2(GSVA,response,.f=function(.x,.y){
#     .x %>%
#       dplyr::inner_join(.y, by = "Run") %>%
#       dplyr::select(-Run) %>%
#       dplyr::mutate(Response = as.factor(Response))
#   })) %>%
#   dplyr::select(data.ready) %>%
#   tidyr::unnest() -> data.ready
# # data.ready <- na.omit(data.ready)
# data.ready[is.na(data.ready)] <- 0
# 
# model.train <- glm( as.formula(formula.train), data = data.ready, family = binomial)
# 
# # Make predictions
# 
# .data_for_logistic %>%
#   dplyr::mutate(n=purrr::map(response,.f=function(.x){
#     nrow(.x)
#   })) %>%
#   tidyr::unnest(n) %>%
#   dplyr::filter(n>5) %>%
#   dplyr::select(-n) %>%
#   dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model.train)) %>%
#   dplyr::select(-response,-GSVA) %>%
#   tidyr::unnest() -> validation_auc
