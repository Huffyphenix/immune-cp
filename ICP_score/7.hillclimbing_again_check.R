
# re run the hill climbing to filtered features ---------------------------

################################
# construct logistics model for GSVA score of all possible features
# GSVA was generate in each dataset
# train data for all data, is Riaz
# seperate the pre- and on-treatment data
# filter 5 features from 78 features filtered from TCGA

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
basic_path <- file.path("S:/坚果云/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
exp_score_path <- file.path(immune_res_path,"ICP_score.new")
gsva_score_path <-  file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")
res_path <- file.path(exp_score_path,"logistic_model_predict_Response/hill_climbing_191223")

# load data ---------------------------------------------------------------
load( file.path(exp_score_path,"clinical_data_for_logistic.rda")) -> data_for_logistic

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
    dplyr::filter(! cv %in% CV) %>%
    # dplyr::select(-usage) %>%
    dplyr::mutate(Response = as.factor(Response)) -> data.ready
  data.ready[is.na(data.ready)] <- 0
  colnames(data.ready) <- gsub(" ",".",colnames(data.ready))
  iteration_n <- iteration
  res <- apply(array(1:iteration_n,dim = iteration_n),1,FUN=fn_feature_filter, data = data.ready, times = times, test_n = test_n, train_n = train_n, CV = CV)
  # print("res get done")
  
  res
  # print("res output done")
}

fn_feature_filter <- function(x, data, times, test_n = 3, train_n = 6, CV){
  # data prepare
  test_set <- sample_n(data %>% dplyr::filter(Response == "yes"),test_n) %>%
    rbind(sample_n(data %>% dplyr::filter(Response == "no"),test_n)) 
  
  data %>%
    dplyr::filter(!Run %in% test_set$Run) -> data_1
  
  train_set <- sample_n(data_1 %>% dplyr::filter(Response == "yes"),train_n) %>%
    rbind(sample_n(data_1 %>% dplyr::filter(Response == "no"),train_n)) 
  train_set %>%
    dplyr::mutate(class = "train") %>%
    rbind(test_set %>%
            dplyr::mutate(class = "test")) %>%
    dplyr::select(Run,class) %>%
    readr::write_tsv(file.path(res_path,"train_test_record",paste("I",x,"CV",CV,"train_test_set.tsv")))
  
  train_set <- train_set %>%
    dplyr::select(-Run)
  
  test_set <- test_set %>%
    dplyr::select(-Run)
  
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
    dplyr::group_by(features) %>%
    dplyr::mutate(score_all = sum(score_individual)) %>%
    dplyr::arrange(desc(score_all)) %>%
    dplyr::ungroup() %>%
    dplyr::select(features, score_all) %>%
    unique() -> .score
  
  for_score %>%
    dplyr::mutate(score_individual.class = ifelse(test_auc >= pos_auc, "success", "not_sure")) %>%
    dplyr::mutate(score_individual.class = ifelse(test_auc <= neg_auc, "failures", score_individual.class)) %>%
    dplyr::group_by(CV, features,score_individual.class) %>%
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
  set.seed(12356)
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
# if(blockade.x != "anti鈥揚D-1"){
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

# get AUC on validation data ---------------
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

########################################## logistic regression iteration result -----------------------
data_for_logistic %>%
  dplyr::filter(Author=="Riaz") %>%
  dplyr::mutate(sample_group = purrr::map2(response,GSVA,fn_select_train_test,percent = 0.7)) %>%
  dplyr::mutate(data = purrr::pmap(list(response,GSVA,sample_group),.f=function(.x,.y,.z){
    .x %>%
      dplyr::inner_join(.z,by="Run") %>%
      tidyr::nest(-usage,.key="response") -> response
    .y %>%
      dplyr::inner_join(.z,by="Run") %>%
      tidyr::nest(-usage,.key="GSVA") -> GSVA
    response %>%
      dplyr::inner_join(GSVA, by="usage")
  })) %>%
  dplyr::select(-response,-GSVA,-sample_group) %>%
  tidyr::unnest() -> Riaz.data 

data_for_logistic %>%
  dplyr::mutate(usage = c("validation","validation","validation","train+test","validation")) %>%
  rbind(Riaz.data) %>%
  dplyr::mutate(n = purrr::map2(Author,response,.f=function(.y,.x){
    print(.y)
    .x$Response %>%
      table()%>% 
      as.data.frame() %>%
      tidyr::spread(key=".",value="Freq") %>% 
      dplyr::as.tbl()})) %>%
  tidyr::unnest(n) -> data_for_logistic

############## logistic regression ---------------
# stage=c("pre-treatment")
dir=""
test_n=6
train_n=12
succss = 0
# fn_Run <- function(stage,dir,test_n, train_n){
for (times in 1:1000) {
  
  
  data_for_logistic  -> .data_for_logistic
  
  # SAMPLE INTO 5 FOLD
  .data_for_logistic %>%
    dplyr::filter(Biopsy_Time=="on-treatment" & Author=="Riaz") %>%
    dplyr::select(response) %>%
    tidyr::unnest() -> training_data
  fold_res <- tibble::tibble()
  training_data %>%
    dplyr::filter(Response == "yes") -> traing_yes
  training_data %>%
    dplyr::filter(Response == "no") -> traing_no
  n_per_cv_yes = round(nrow(traing_yes)/5)
  n_per_cv_no = round(nrow(traing_no)/5)
    
  for(cv in 1:5){
    traing_yes %>%
      dplyr::filter(! Run %in% fold_res$Run) -> traing_yes
    traing_no %>%
      dplyr::filter(! Run %in% fold_res$Run) -> traing_no
    if(cv < 5){
      fold_res_tmp <- tibble::tibble(Run = c(sample(traing_yes$Run,n_per_cv_yes), sample(traing_no$Run,n_per_cv_no)), fold = cv)
    } else{
      fold_res_tmp <- tibble::tibble(Run = c(traing_yes$Run, traing_no$Run), fold = cv)
    }
    fold_res <- rbind(fold_res, fold_res_tmp)
  }
  repeat_res_path <- file.path(res_path,paste("repeat",i))
  dir.create(repeat_res_path)
  fold_res %>% 
    readr::write_tsv(file.path(repeat_res_path,"cross_fold_info.tsv"))
  
  paste(paste(paste(.data_for_logistic$Author,.data_for_logistic$usage, sep="_"),collapse = "\t"),
        "Features",
        "Times","SuccssTime",sep="\t") %>%
    readr::write_lines(file.path(res_path,dir,"logistic_record.txt"),sep = "\n",append = TRUE)
  ##################### repeat 100 times get the best features
  while(succss <= 10){
    times=times+1
    n=5
    logistic_feature_res_nrepeat <- list()
    for(i in 1:n){
      .data_for_logistic %>%
        dplyr::filter(Biopsy_Time=="on-treatment" & Author=="Riaz") %>%
        dplyr::mutate(test_n = test_n, train_n = train_n) %>%
        dplyr::mutate(itertion_res = purrr::pmap(list(Cancer.y,GSVA, response, test_n, train_n),.f=fn_logistic,iteration = 500, times = 10, CV = i, fold_info)) %>%
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
        dplyr::select(-itertion_res_datafrme) -> logistic_feature_Selected_nrepeat[[i]]
    }
    ############## features selection
    fn_overlap_select(logistic_feature_Selected_nrepeat,top_n = 10, overlpa_n = 5) -> final_feature
    
    # final_feature %>%
    #   tidyr::unnest() %>%
    #   readr::write_tsv(file.path(res_path,"class_metastic_type-and-gsvascore_with_class_metastic_type","repeat_times_feature",paste(j,"final_feature.tsv",sep="_")))
    
    ############## draw picture, using stepAIC
    # get model from 70% train and use it on 30% Test and validation data
    .data_for_logistic %>%
      dplyr::filter(usage == "train") %>% 
      dplyr::mutate(data.ready = purrr::map2(GSVA,response,.f=function(.x,.y){
        .x %>%
          dplyr::inner_join(.y, by = "Run") %>%
          dplyr::select(-Run) %>%
          dplyr::mutate(Response = as.factor(Response))
      })) %>%
      dplyr::select(data.ready) %>%
      tidyr::unnest() -> data.ready
    # data.ready <- na.omit(data.ready)
    data.ready[is.na(data.ready)] <- 0
    colnames(data.ready) <- gsub(" ",".",colnames(data.ready))
    colnames(data.ready) <- gsub("-",".",colnames(data.ready))
    
    formula.train <- paste("Response~", paste(final_feature$topn_count[[1]]$features,collapse = "+"))
    model.train <- glm( as.formula(formula.train), data = data.ready, family = binomial)
    step.model.train <- stepAIC(model.train, direction = "backward",trace = F)
    final_features <- rownames(as.data.frame(summary(step.model.train)$coefficients))[-1]
    
    # Make predictions
    
    .data_for_logistic %>%
      dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=step.model.train)) %>%
      dplyr::select(-response,-GSVA) %>%
      tidyr::unnest() -> validation_auc
    
    paste(paste(signif(validation_auc$auc,2),collapse = "\t"),
          paste(final_features,collapse = "/"),
          times,succss,sep="\t") %>%
      readr::write_lines(file.path(res_path,dir,"logistic_record.txt"),sep = "\n",append = TRUE)
    if(min(validation_auc$auc)>=0.6){
      succss = succss + 1
      model.train %>%
        readr::write_rds(file.path(res_path,stage,paste("Final_model",succss,".rds.gz",sep="_")),compress = "gz")
      summary(model.train)$coefficients %>%
        as.data.frame() %>%
        readr::write_tsv(file.path(res_path,stage,paste("Final_model_coef",succss,".rds.gz",sep="_")))
      
      # draw picture
      validation_auc %>%
        dplyr::inner_join(Response_statistic,by=c("Author")) %>%
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
      ggsave(file.path(res_path,dir,paste("ROC_plot",Times,"png",sep="_")),device = "png",width = 5, height = 5)
      ggsave(file.path(res_path,dir,paste("ROC_plot",Times,"pdf",sep="_")),device = "pdf",width = 5, height = 5)
      
    } else{
      succss = succss
      
    }
    res <- readr::read_tsv(file.path(res_path,dir,"logistic_record.txt"))%>%
      dplyr::select(-Features,-SuccssTime)
    res %>%
      tidyr::gather(-Times,key="DataSets",value="AUC") %>%
      dplyr::mutate(AUC = as.numeric(AUC)) -> for_draw
    for_draw <- within(for_draw,DataSets<-factor(DataSets,levels = unique(for_draw$DataSets)))
    for_draw %>%
      dplyr::mutate(AUC_group=ifelse(AUC>0.7,"red","black")) %>%
      ggplot(aes(x=DataSets,y=AUC)) +
      geom_jitter(aes(color=AUC_group)) +
      scale_color_manual(values=c("black","red")) +
      theme(
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)
      )
    ggsave(file.path(res_path,"logistic_res.png"),device = "png",height = 5,width = 8)
  }
  succss
}