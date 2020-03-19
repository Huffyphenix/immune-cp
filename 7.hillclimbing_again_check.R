
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
# basic_path <- file.path("S:/坚果云/我的坚果云")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
exp_score_path <- file.path(immune_res_path,"ICP_score.new")
gsva_score_path <-  file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")
res_path <- file.path(exp_score_path,"logistic_model_predict_Response/hill_climbing_191226")

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
fn_logistic <- function(Cancer_type,exp, response, test_n, train_n, iteration, CV, times = 10, fold_info, out_dir){
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
  res <- apply(array(1:iteration_n,dim = iteration_n),1,FUN=fn_feature_filter, data = data.ready, times = times, test_n = test_n, train_n = train_n, CV = CV, out_dir=out_dir)
  # print("res get done")
  
  res
  # print("res output done")
}

fn_feature_filter <- function(x, data, times, test_n = 3, train_n = 6, CV, out_dir){
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
    readr::write_tsv(file.path(out_dir,paste("I",x,"CV",CV,"train_test_set.tsv",sep=".")))
  
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
fn_feature_select_accross_CV <- function(.data, top_n = 10, overlpa_n = 5, out_dir){
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
    readr::write_tsv(file.path(out_dir, "binomP_accross_CV.tsv"))
  success %>%
    dplyr::filter(`1`>=0.6 & `2`>=0.6 & `3` >=0.6 & `4` >=0.6 &`5`>=0.6) -> features_all_success_0.6
  success %>%
    readr::write_tsv(file.path(out_dir, "successFreq_accross_CV.tsv"))
  list(features_all_binomP_sig = features_all_binomP_sig$features,
       features_all_success = features_all_success_0.6$features)
}
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

############## logistic regression ---------------
# stage=c("pre-treatment")
dir=""
test_n=4
train_n=8
succss = 0
# fn_Run <- function(stage,dir,test_n, train_n){
fn_repeat <- function (times) {
  
  readr::write_lines(paste("repreat",times), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)
  start_point <- date()
  
  readr::write_lines(paste(times,start_point), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)
  
  # print(start_point)
  data_for_logistic %>%
    dplyr::mutate(n=purrr::map(response,.f=function(.x){
      nrow(.x)
    })) %>%
    tidyr::unnest(n) %>%
    dplyr::filter(n>5) %>%
    dplyr::select(-n) -> .data_for_logistic
  
  # SAMPLE INTO 5 FOLD
  readr::write_lines(paste(times,"# SAMPLE INTO 5 FOLD"), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)
  
  # print("# SAMPLE INTO 5 FOLD")
  .data_for_logistic %>%
    dplyr::filter(Author=="Riaz") %>%
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
  repeat_res_path <- file.path(res_path,paste("repeat",times,sep="_"))
  dir.create(repeat_res_path, recursive = T)
  fold_res %>% 
    readr::write_tsv(file.path(repeat_res_path,"cross_fold_info.tsv"))
  
  # if(times==1){
  #   paste(paste(paste(.data_for_logistic$Author,.data_for_logistic$Biopsy_Time, sep="_"),collapse = "\t"),
  #         "Features",
  #         "Repeats","Model_type","Success_test","Success_valid",sep="\t") %>%
  #     readr::write_lines(file.path(res_path,"logistic_record.txt"),sep = "\n",append = TRUE)
  # }
  
  ##################### repeat 100 times get the best features
  .data_for_logistic %>%
    dplyr::filter(Author=="Riaz") %>%
    dplyr::select(-GSVA,-Biopsy_Time)  %>%
    tidyr::unnest() %>%
    tidyr::nest(- Cancer.y,-blockade,-Author) %>%
    dplyr::rename("response"="data") -> train_response
  
  .data_for_logistic %>%
    dplyr::filter(Author=="Riaz") %>%
    dplyr::select(-response,-Biopsy_Time)  %>%
    tidyr::unnest() %>%
    tidyr::nest(- Cancer.y,-blockade,-Author) %>%
    dplyr::rename("GSVA"="data") -> train_GSVA
  
  readr::write_lines(paste(times,"# 500 iteration among 5 CV"), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)
  readr::write_lines(paste(times,"# hill climbing"), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)
  # print("# 500 iteration among 5 CV")
  # print("# hill climbing")
    n=5
    logistic_feature_res_nrepeat <- list()
    for(i in 1:n){
      train_response %>%
        dplyr::inner_join(train_GSVA) %>%
        dplyr::mutate(test_n = test_n, train_n = train_n) %>%
        dplyr::mutate(itertion_res = purrr::pmap(list(Cancer.y,GSVA, response, test_n, train_n),.f=fn_logistic,iteration = 500, CV = i, times = 20, fold_info=fold_res, out_dir = repeat_res_path)) %>%
        dplyr::select(-response,-GSVA)-> logistic_feature_res_nrepeat[[i]]
    }
    ############## score for regression iteration features
    readr::write_lines(paste(times,"# get hill climbing score"), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)
    
    # print("# get hill climbing score")
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
    readr::write_lines(paste(times,"# features selection"), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)
    
    # print("# features selection")
    fn_feature_select_accross_CV(logistic_feature_Selected_nrepeat,top_n = 10, overlpa_n = 5, out_dir = repeat_res_path) -> final_feature
    
    final_feature %>%
      readr::write_rds(file.path(repeat_res_path,"final_feature.rds.gz"),compress = "gz")
    
    ############## draw picture, using stepAIC
    # get model from 70% train and use it on 30% Test and validation data
    readr::write_lines(paste(times,"# perform of features selected among 5 CV test data"), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)
    
    # print("# perform of features selected among 5 CV test data")
    TEST_AUC <- tibble::tibble()
    for (i in 1:n) {
      train_response %>%
        dplyr::inner_join(train_GSVA) %>% 
        dplyr::mutate(data.ready = purrr::map2(GSVA,response,.f=function(.x,.y){
          .x %>%
            dplyr::inner_join(.y, by = "Run") %>%
            dplyr::inner_join(fold_res, by= "Run") %>%
            dplyr::filter(! fold %in% i) %>%
            dplyr::select(-fold) %>%
            dplyr::select(-Run) %>%
            dplyr::mutate(Response = as.factor(Response))
        })) %>%
        dplyr::select(data.ready) %>%
        tidyr::unnest() -> data.ready 
      data.ready[is.na(data.ready)] <- 0
      
      # Make predictions
      train_response %>%
        dplyr::inner_join(train_GSVA) %>%
        dplyr::mutate(GSVA = purrr::map(GSVA,.f=function(.x){
          .x %>%
            dplyr::inner_join(fold_res, by= "Run") %>%
            dplyr::filter(fold %in% i) %>%
            dplyr::select(-fold) 
        })) %>%
        dplyr::mutate(response = purrr::map(response,.f=function(.x){
          .x %>%
            dplyr::inner_join(fold_res, by= "Run") %>%
            dplyr::filter(fold %in% i) %>%
            dplyr::select(-fold) 
        }))  -> data.ready.test
      
      # select p all significant
      if(length(final_feature$features_all_binomP_sig)>0){
        formula.train_binom <- paste("Response~", paste(final_feature$features_all_binomP_sig,collapse = "+"))
        model.train_binom <- glm( as.formula(formula.train_binom), data = data.ready, family = binomial)
        step.model.train_binom <- stepAIC(model.train_binom, direction = "backward",trace = F)
        final_feature_binom <- rownames(as.data.frame(summary(step.model.train_binom)$coefficients))[-1]
        
        data.ready.test %>%
          dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model.train_binom)) %>%
          dplyr::select(-response,-GSVA) %>%
          tidyr::unnest() %>%
          dplyr::mutate(CV = i, model = "from_binom") %>%
          dplyr::select(-roc_data) -> TEST_auc_binom
      }else{
        TEST_auc_binom <- tibble::tibble(Cancer.y="melanoma",blockade="anti–PD-1",Author="Riaz",auc=0,CV= i,model= "from_binom")
      }
      
      # select frequency all > 0.6
      if(length(final_feature$features_all_success)>0){
        formula.train_success <- paste("Response~", paste(final_feature$features_all_success,collapse = "+"))
        model.train_success <- glm( as.formula(formula.train_success), data = data.ready, family = binomial)
        step.model.train_success <- stepAIC(model.train_success, direction = "backward",trace = F)
        final_features_success <- rownames(as.data.frame(summary(step.model.train_success)$coefficients))[-1]
        
        data.ready.test %>%
          dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model.train_success)) %>%
          dplyr::select(-response,-GSVA) %>%
          tidyr::unnest() %>%
          dplyr::mutate(CV = i, model = "from_successFreq")  %>%
          dplyr::select(-roc_data) -> TEST_auc_success
      }else{
        TEST_auc_success <- tibble::tibble(Cancer.y="melanoma",blockade="anti–PD-1",Author="Riaz",auc=0,CV= i,model= "from_successFreq")
      }
      
      TEST_AUC <- rbind(TEST_AUC, TEST_auc_binom, TEST_auc_success)
    }
    TEST_AUC %>% 
      readr::write_tsv(file.path(repeat_res_path, "TEST_AUC.tsv"))
    if(mean(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc) > mean(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc)){
      model_type <- "from_binom"
      formula.train <- paste("Response~", paste(final_feature$features_all_binomP_sig,collapse = "+"))
      mean_test_auc <- mean(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc)
      final_features <- final_feature$features_all_binomP_sig
      if (all(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc>=0.6)) {
         success_test <- "yes"
      } else{
        success_test <- "no"
      }
    } else if (mean(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc) > mean(TEST_AUC %>% dplyr::filter(model=="from_binom") %>% .$auc)){
      model_type <- "from_successFreq"
      formula.train <- paste("Response~", paste(final_feature$features_all_success,collapse = "+"))
      mean_test_auc <- mean(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc)
      final_features <- final_feature$features_all_success
      if (all(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc>=0.6)) {
        success_test <- "yes"
      } else{
        success_test <- "no"
      }
    } else{
      formula.train <- paste("Response~1")
      mean_test_auc <- mean(TEST_AUC %>% dplyr::filter(model=="from_successFreq") %>% .$auc)
      final_features <- "NNa"
      model_type <- "NNa"
      success_test <- "no"
    }
    train_response %>%
      dplyr::inner_join(train_GSVA) %>%
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
    
    model.train <- glm( as.formula(formula.train), data = data.ready, family = binomial)
    model.train %>%
      readr::write_rds(file.path(repeat_res_path,"logic.model.rds.gz"))
    
    # Make predictions
    
    .data_for_logistic %>%
      dplyr::mutate(n=purrr::map(response,.f=function(.x){
        nrow(.x)
      })) %>%
      tidyr::unnest(n) %>%
      dplyr::filter(n>5) %>%
      dplyr::select(-n) %>%
      dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model.train)) %>%
      dplyr::select(-response,-GSVA) %>%
      tidyr::unnest() -> validation_auc
    if(all(validation_auc$auc >= 0.6)){
      success_valid<-"yes"
    }else{
      success_valid<-"no"
    }
    # paste(paste(signif(validation_auc$auc,2),collapse = "\t"),
    #       paste(final_features,collapse = "/"),
    #       times,model_type,success_test,success_valid,sep="\t") %>%
    #   readr::write_lines(file.path(res_path,dir,"logistic_record.txt"),sep = "\n",append = TRUE)
   
      summary(model.train)$coefficients %>%
        as.data.frame() %>%
        readr::write_tsv(file.path(repeat_res_path,paste("Final_model_coef",".rds.gz",sep="_")))
      
      # draw picture
      validation_auc %>%
        dplyr::inner_join(Response_statistic,by=c("Author","Biopsy_Time")) %>%
        dplyr::mutate(Biopsy_Time.short=strsplit(Biopsy_Time,split = "-")[[1]][1]) %>%
        dplyr::mutate(group = paste(Author,Biopsy_Time.short, blockade.y,paste("(",Response,"*/",`non-Response`,"**);",sep=""),"AUC","=",signif(auc,2))) %>%
        dplyr::select(roc_data,group) %>%
        tidyr::unnest() -> plot_ready
      plot_ready %>%
        ggplot(aes(x=Specificity,y=Sensitivity)) +
        geom_path(aes(color=group)) + 
        scale_x_reverse()  +
        my_theme +
        theme(
          legend.position = c(0.75,0.25),
          legend.title = element_blank()
        )
      ggsave(file.path(repeat_res_path,paste("ROC_plot_repeat",times,".png",sep="_")),device = "png",width = 5, height = 5)
      ggsave(file.path(repeat_res_path,paste("ROC_plot_repeat",times,".pdf",sep="_")),device = "pdf",width = 5, height = 5)
      
      end_point <- date()
      if(times==1){
        tibble::tibble(times=times,
                       AUC=signif(validation_auc$auc,2),
                       features=paste(final_features,collapse = "/"),
                       model_type=model_type,success_test=success_test,success_valid=success_valid,
                       dataset=paste(.data_for_logistic$Author,.data_for_logistic$Biopsy_Time, sep="_")) %>%
          tidyr::spread(key = "dataset", value="AUC") %>%
          as.data.frame() %>%
          colnames() %>%
          as.character() %>%
          paste0(collapse = "\t") %>%
          readr::write_lines(file.path(res_path,"logistic.record"),sep = "\n",append = TRUE)
        
        tibble::tibble(times=times,
                       AUC=signif(validation_auc$auc,2),
                       features=paste(final_features,collapse = "/"),
                       model_type=model_type,success_test=success_test,success_valid=success_valid,
                       dataset=paste(.data_for_logistic$Author,.data_for_logistic$Biopsy_Time, sep="_")) %>%
          tidyr::spread(key = "dataset", value="AUC") %>%
          as.data.frame() %>%
          as.character() %>%
          paste0(collapse = "\t") %>%
          readr::write_lines(file.path(res_path,"logistic.record"),sep = "\n",append = TRUE)
      }else{
        tibble::tibble(times=times,
                       AUC=signif(validation_auc$auc,2),
                       features=paste(final_features,collapse = "/"),
                       model_type=model_type,success_test=success_test,success_valid=success_valid,
                       dataset=paste(.data_for_logistic$Author,.data_for_logistic$Biopsy_Time, sep="_")) %>%
          tidyr::spread(key = "dataset", value="AUC") %>%
          as.data.frame() %>%
          as.character() %>%
          paste0(collapse = "\t") %>%
          readr::write_lines(file.path(res_path,"logistic.record"),sep = "\n",append = TRUE)
      }
      
      tibble::tibble(times=times,
                     AUC=signif(validation_auc$auc,2),
                     features=paste(final_features,collapse = "/"),
                     model_type=model_type,success_test=success_test,success_valid=success_valid,
                     dataset=paste(.data_for_logistic$Author,.data_for_logistic$Biopsy_Time, sep="_")) 
        
}

library(multidplyr)
cl <- parallel::detectCores()
cluster <- multidplyr::new_cluster(10)

multidplyr::cluster_assign(cluster,"fn_logistic"=fn_logistic)  
multidplyr::cluster_assign(cluster,"fn_feature_filter"=fn_feature_filter) 
multidplyr::cluster_assign(cluster,"fn_auc"=fn_auc) 
multidplyr::cluster_assign(cluster,"fn_feature_score"=fn_feature_score) 
multidplyr::cluster_assign(cluster,"fn_formula_split"=fn_formula_split) 
multidplyr::cluster_assign(cluster,"fn_feature_selection"=fn_feature_selection) 
multidplyr::cluster_assign(cluster,"fn_feature_select_accross_CV"=fn_feature_select_accross_CV) 
multidplyr::cluster_assign(cluster,"fn_auc_on_validation"=fn_auc_on_validation) 
multidplyr::cluster_assign(cluster,"fn_repeat"=fn_repeat) 


multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_library(cluster,"pROC") 
multidplyr::cluster_library(cluster,"tidyverse") 
multidplyr::cluster_library(cluster,"caret") 
multidplyr::cluster_library(cluster,"glmnet") 
multidplyr::cluster_library(cluster,"broom") 
multidplyr::cluster_library(cluster,"ggplot2") 
multidplyr::cluster_library(cluster,"pROC") 
multidplyr::cluster_library(cluster,"MASS")

multidplyr::cluster_assign(cluster,"res_path"=res_path)
multidplyr::cluster_assign(cluster,"data_for_logistic"=data_for_logistic)
multidplyr::cluster_assign(cluster,"test_n"=test_n) 
multidplyr::cluster_assign(cluster,"train_n"=train_n) 
multidplyr::cluster_assign(cluster,"Response_statistic"=Response_statistic)
multidplyr::cluster_assign(cluster,"my_theme"=my_theme)


tibble::tibble(x=c(1:10000)) %>%
  dplyr::group_by(x) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(
    final_features = purrr::map(x,fn_repeat)
  ) %>%
  dplyr::collect() -> final_feature_res

final_feature_res %>%
  tidyr::unnest() %>%
  tidyr::spread(key="dataset",value="AUC") -> final_feature_res

readr::write_lines(paste("# result output"), path = file.path(res_path,"nohup_out.record"),sep = "\n",append = TRUE)

print("result output")
final_feature_res %>%
  readr::write_rds(file.path(res_path,"final_feature_res_5_CV_AUC.rds.gz"), compress = "gz")

parallel::stopCluster(cluster)
