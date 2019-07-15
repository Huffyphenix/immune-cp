################################
# construct logistics model for GSVA score of all possible features
# GSVA was generate in each dataset
# for only pre-treatment samples prediction
# train data for pre-treatment data, is Riaz, pre-treatment

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

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
# res_path <- file.path(immune_res_path,"ICP_score/3.logistic-regression-clinical/by_gsva_score")
# score_path <- file.path(immune_res_path,"ICP_score/2.2.Clinical_validation-GSVA-ICPs_exp_site_5_feature")
score_path <- file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")
res_path <- file.path(score_path,"logistic_model_predict_Response/all_features")

# load data ---------------------------------------------------------------
gsva.score <- readr::read_rds(file.path(score_path,"ICP_GSVA_score_all-possible-features_all-specific.rds.gz"))

# sample data classfication --------------------------------------------------
sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::inner_join(readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","anti_ICB_study_info.txt")) %>%
                      dplyr::select(`data ID`,Author) %>%
                      unique(), by="data ID")

sample_info %>%
  dplyr::select(Run,Cancer.y,blockade,blockade,Biopsy_Time,Author,Response) %>%
  unique() -> Run_pubmed.id

# gsva.score %>%
#   dplyr::inner_join(Run_pubmed.id, by = "Run") %>%
#   # dplyr::filter(symbol %in% gene_list$symbol) %>%
#   dplyr::select(-blockade,-Cancer_type) %>%
#   tidyr::nest(-Cancer.y) -> exp_data.gather.genelist
# 
# exp_data.gather.genelist %>%
#   dplyr::rename("data_spread" = "data") -> exp_data.nest.genelist

Run_pubmed.id %>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no")) %>%
  tidyr::nest(Run,Response,.key = "response") %>%
  dplyr::inner_join(gsva.score, by = c("Cancer.y","blockade","Biopsy_Time","Author")) -> data_for_logistic

Response_statistic <- readr::read_tsv(file.path(score_path,"Response_statistic_each_dataset.tsv")) %>%
  dplyr::select(-Cancer.y)

############## fnctions to do logistic regression iteration ##################
fn_logistic <- function(Cancer_type,exp, response, sample_group, test_n, train_n, iteration, times = 15 ){
  # apply(exp[,-1],2,mad) %>% # mad() Compute the Median Absolute Deviation
  #   sort() %>%
  #   rev() %>%
  #   .[1:40] %>%
  #   names() -> top20_mad
  # print(Cancer_type)
  exp %>%
    dplyr::inner_join(response, by = "Run") %>%
    dplyr::inner_join(sample_group, by = "Run") %>%
    dplyr::filter(usage == "train") %>%
    dplyr::select(-Run, -usage) %>%
    dplyr::mutate(Response = as.factor(Response)) -> data.ready
  data.ready[is.na(data.ready)] <- 0
  colnames(data.ready) <- gsub(" ",".",colnames(data.ready))
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
    dplyr::select(features,binom.p) %>%
    tidyr::unnest() -> binomtest.p
  
  binomtest.p %>%
    dplyr::inner_join(.score,by="features")
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
      dplyr::select(Cancer.y, blockade,Biopsy_Time,Author, features, score_all,binom.p) %>% 
      dplyr::arrange(desc(binom.p),score_all) %>%
      dplyr::mutate(rank = 1:nrow(.))-> topn_tmp 
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
    dplyr::group_by(Cancer.y, blockade,Biopsy_Time,Author, features) %>%
    dplyr::mutate(rank_sum = sum( rank)) %>%
    dplyr::select(Cancer.y, blockade,Biopsy_Time,Author, features,rank_sum) %>%
    dplyr::ungroup() %>%
    unique() %>%
    dplyr::arrange(rank_sum) %>%
    # dplyr::filter(rank_sum < top_n*10) %>%
    # dplyr::select(Cancer.y, Cancer_type, features, rank_sum) %>%
    tidyr::nest(features,rank_sum) %>%
    dplyr::mutate(topn_count = purrr::map(data,function(.x){
      .x %>%
        top_n(top_n, rank_sum)
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

##################################### performing of the feature selected #####################################
# draw picture
#fn_logistc_draw <- function(Cancer.y,blockade,response,data_spread,sample_group,blockade,topn_count,dir,stepAIC=T){
# if(blockade.x != "anti–PD-1"){
#   sample_group %>%
#     dplyr::mutate(usage = "test") -> sample_group
# }
## data prepare
# data_spread %>%
#   dplyr::inner_join(response, by = "Run") %>%
#   dplyr::inner_join(sample_group, by = "Run") %>%
#   dplyr::filter(usage == "test") %>%
#   dplyr::select(-Run, -usage) %>%
#   dplyr::mutate(Response = as.factor(Response)) -> data.ready
# # data.ready <- na.omit(data.ready)
# data.ready[is.na(data.ready)] <- 0
# colnames(data.ready) <- gsub(" ",".",colnames(data.ready))

# fomula prepare
# for(i in 1:nrow(topn_count)){
#   if(i == 1){
#     formula <- paste("Response", topn_count$features[i], sep = "~")
#   } else{
#     formula <- paste(formula, topn_count$features[i], sep = "+")
#   }
# }
# for(i in 1:nrow(topn_count)){
#   if(i == 1){
#     formula <- paste("Response", topn_count$.[i], sep = "~")
#   } else{
#     formula <- paste(formula, topn_count$.[i], sep = "+")
#   }
# }
## do logistic regression
# model <- glm( as.formula(formula), data = data.ready, family = binomial)
# 
# ## stepwise regression model
# if(stepAIC){
#   step.model <- stepAIC(model, direction = "backward",trace = F)
#   model <- step.model
# } else{
#   model <-model
# }


# Make predictions
#   probabilities <- model %>% predict(type = "response")
#   predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
#   observed.classes <- data.ready$Response
#   # Model accuracy
#   accuracy <- mean(predicted.classes == data.ready$Response)
#   error <- mean(predicted.classes != data.ready$Response)
#   
#   # table(observed.classes, predicted.classes)
#   # confusionMatrix(as.factor(predicted.classes), observed.classes,
#   #                 positive = "yes")
#   ## ROC
#   res.roc <- roc(observed.classes, probabilities)
#   png(file.path(res_path,dir,paste(Cancer.y,blockade.x,"by-feature-from",Cancer_type,blockade,"ROC.png", sep="_")),width = 400, height = 300)
#   plot.roc(res.roc, print.auc = TRUE, print.thres = "best")
#   dev.off()
#   
#   formula <- model$formula
#   fn_formula_split(Reduce(paste, deparse(formula))) %>%
#     readr::write_tsv(file.path(res_path,dir,paste(Cancer.y,blockade.x,"by-feature-from",blockade,"features.tsv", sep="_")))
# }

# get AUC on validation data
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
########################################## logistic regression iteration result ------------------------------------------
data_for_logistic %>%
  dplyr::mutate(sample_group = purrr::map2(response,GSVA,fn_select_train_test,percent = 0.7)) -> data_for_logistic

data_for_logistic %>%
  dplyr::mutate(type = ifelse(Author=="Riaz", "train", "validation")) -> data_for_logistic
############## logistic regression
fn_Run <- function(stage,test_n, train_n){
  succss = 0 
  times = 0
  
  data_for_logistic %>%
    dplyr::filter(Biopsy_Time == stage) -> .data_for_logistic
  
  paste("auc.train","auc.test",
        paste(paste("auc_validation",.data_for_logistic$Author,sep="_"),collapse = " "),
        "Features",
        "Times","SuccssTime") %>%
    readr::write_lines(file.path(res_path,stage,"logistic_record.txt"),sep = "\n",append = TRUE)
  ##################### repeat 100 times get the best features
  while(succss < 10){
    times=times+1
    n=10
    logistic_feature_res_nrepeat <- list()
    for(i in 1:n){
      .data_for_logistic %>%
        dplyr::filter(type == "train") %>%
        dplyr::mutate(test_n = test_n, train_n = train_n) %>%
        # dplyr::mutate(test_n = 3, train_n = 4) %>%
        # dplyr::filter(! blockade == "anti–PD-1、CTLA-4") %>% # last two data sets of melanoma dont have enough response observation to train
        dplyr::mutate(itertion_res = purrr::pmap(list(Cancer.y,GSVA, response, sample_group, test_n, train_n),.f=fn_logistic,iteration = 500, times = 10)) %>%
        dplyr::select(-response,-GSVA)-> logistic_feature_res_nrepeat[[i]]
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
        dplyr::mutate(selection = purrr::map(itertion_res_datafrme,fn_feature_score,pos_auc = 0.8, neg_auc = 0.4)) %>%
        dplyr::select(-itertion_res_datafrme, -sample_group) -> logistic_feature_Selected_nrepeat[[i]]
    }
    ############## features selection
    fn_overlap_select(logistic_feature_Selected_nrepeat,top_n = 10, overlpa_n = 5) -> final_feature
    
    # final_feature %>%
    #   tidyr::unnest() %>%
    #   readr::write_tsv(file.path(res_path,"class_metastic_type-and-gsvascore_with_class_metastic_type","repeat_times_feature",paste(j,"final_feature.tsv",sep="_")))
    
    ############## draw picture, using stepAIC
    # get model from 70% train and use it on 30% Test and validation data
    .data_for_logistic %>%
      dplyr::filter(type == "train") %>%
      dplyr::mutate(data.ready = purrr::pmap(list(GSVA,response,sample_group),.f=function(.x,.y,.z){
        .x %>%
          dplyr::inner_join(.y, by = "Run") %>%
          dplyr::inner_join(.z, by = "Run") %>%
          dplyr::filter(usage == "train") %>%
          dplyr::select(-Run, -usage) %>%
          dplyr::mutate(Response = as.factor(Response))
      })) %>%
      dplyr::select(data.ready) %>%
      tidyr::unnest() -> data.ready
    # data.ready <- na.omit(data.ready)
    data.ready[is.na(data.ready)] <- 0
    colnames(data.ready) <- gsub(" ",".",colnames(data.ready))
    
    formula.train <- paste("Response~", paste(final_feature$topn_count[[1]]$features,collapse = "+"))
    model.train <- glm( as.formula(formula.train), data = data.ready, family = binomial)
    step.model.train <- stepAIC(model.train, direction = "backward",trace = F)
    
    # Make predictions
    probabilities <- step.model.train %>% predict(type = "response")
    predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
    observed.classes <- data.ready$Response
    auc.train <- roc(observed.classes, probabilities,quiet=TRUE)$auc[[1]]
    
    # 30% test data 
    .data_for_logistic %>%
      dplyr::filter(type == "train") %>%
      dplyr::mutate(data.ready = purrr::pmap(list(GSVA,response,sample_group),.f=function(.x,.y,.z){
        .x %>%
          dplyr::inner_join(.y, by = "Run") %>%
          dplyr::inner_join(.z, by = "Run") %>%
          dplyr::filter(usage == "test") %>%
          dplyr::select(-Run, -usage) %>%
          dplyr::mutate(Response = as.factor(Response))
      })) %>%
      dplyr::select(data.ready) %>%
      tidyr::unnest() -> data.ready.test
    data.ready.test[is.na(data.ready.test)] <- 0
    colnames(data.ready.test) <- gsub(" ",".",colnames(data.ready.test))
    
    # Make predictions
    probabilities <- step.model.train %>% predict(data.ready.test,type = "response")
    predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
    observed.classes <- data.ready.test$Response
    auc.test <- roc(observed.classes, probabilities,quiet=TRUE)$auc[[1]]
    
    .data_for_logistic %>%
      # dplyr::filter(type == "validation") %>%
      dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=step.model.train)) %>%
      dplyr::select(-response,-GSVA, -sample_group) %>%
      tidyr::unnest() -> validation_auc
    
    paste(signif(auc.train,2),signif(auc.test,2),
          paste(signif(validation_auc$auc,2),collapse = " "),
          paste(rownames(as.data.frame(summary(step.model.train)$coefficients))[-1],collapse = "/"),
          times,succss) %>%
      readr::write_lines(file.path(res_path,stage,"logistic_record.txt"),sep = "\n",append = TRUE)
    if(min(auc.train,auc.test)>=0.8){
      if(min(validation_auc$auc)>=0.7){
        succss = succss + 1
        step.model.train %>%
          readr::write_rds(file.path(res_path,stage,paste("Final_model",succss,".rds.gz",sep="_")),compress = "gz")
        summary(step.model.train)$coefficients %>%
          as.data.frame() %>%
          readr::write_tsv(file.path(res_path,stage,paste("Final_model_coef",succss,".rds.gz",sep="_")))
        
        # draw picture
        validation_auc %>%
          dplyr::inner_join(Response_statistic,by=c("Biopsy_Time","Author")) %>%
          dplyr::mutate(group = paste(Author,blockade.y,paste("(",Response,"*/",`non-Response`,"**);",sep=""),"AUC","=",signif(auc,2))) %>%
          dplyr::select(roc_data,group,type) %>%
          tidyr::unnest() -> plot_ready
        plot_ready %>%
          ggplot(aes(x=Specificity,y=Sensitivity)) +
          geom_path(aes(color=group,linetype = type)) + 
          scale_x_reverse()  +
          my_theme +
          theme(
            legend.position = c(0.75,0.25),
            legend.title = element_blank()
          )
        ggsave(file.path(res_path,stage,paste("ROC_plot",succss,"png",sep="_")),device = "png",width = 5, height = 5)
        ggsave(file.path(res_path,stage,paste("ROC_plot",succss,"pdf",sep="_")),device = "pdf",width = 5, height = 5)
        
      } else{
        succss = succss
      }
    }
  }
  succss
}

# cl <- 2
# cluster <- create_cluster(cores = cl)
# tibble::tibble(Biopsy_Time=c("pre-treatment","on-treatment")) %>%
#   partition(cluster = cluster) %>%
#   multidplyr::cluster_library("magrittr") %>%
#   multidplyr::cluster_library("tidyverse") %>%
#   multidplyr::cluster_library("stringr") %>%
#   multidplyr::cluster_library("caret") %>%
#   multidplyr::cluster_library("glmnet") %>%
#   multidplyr::cluster_library("broom") %>%
#   multidplyr::cluster_library("pROC") %>%
#   multidplyr::cluster_library("MASS") %>%
#   multidplyr::cluster_library("ggplot2") %>%
#   multidplyr::cluster_assign_value("fn_DE",fn_logistic) %>%
#   multidplyr::cluster_assign_value("fn_feature_filter",fn_feature_filter) %>%
#   multidplyr::cluster_assign_value("fn_auc",fn_auc) %>%
#   multidplyr::cluster_assign_value("fn_feature_score",fn_feature_score) %>%
#   multidplyr::cluster_assign_value("fn_overlap_select",fn_overlap_select) %>%
#   multidplyr::cluster_assign_value("fn_select_train_test",fn_select_train_test) %>%
#   multidplyr::cluster_assign_value("fn_auc_on_validation",fn_auc_on_validation) %>%
#   multidplyr::cluster_assign_value("my_theme",my_theme) %>%
#   multidplyr::cluster_assign_value("data_for_logistic",data_for_logistic) %>% 
#   multidplyr::cluster_assign_value("res_path",res_path) %>%
#   multidplyr::cluster_assign_value("fn_Run",fn_Run) %>%
#   dplyr::mutate(DE_res = purrr::map(Biopsy_Time,fn_Run)) %>%
#   collect() %>%
#   dplyr::as_tibble() %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-PARTITION_ID, -expr) -> gene_DE_between_high_low.all_features.cox_score.OS
# on.exit(parallel::stopCluster(cluster))
fn_Run(stage="pre-treatment",test_n=3, train_n=6)
# fn_Run(Biopsy_Time="on-treatment")

# save.image(file.path(res_path,"logistic_prediction.rda"))