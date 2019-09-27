################################
# construct logistics model for GSVA score of all possible features
# GSVA was generate in each dataset
# train data for all data, is Riaz
# seperate the pre- and on-treatment data
# features are permutation and combination of 44 features filtered from TCGA

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
exp_score_path <- file.path(immune_res_path,"ICP_score.new")
gsva_score_path <-  file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")
res_path <- file.path(exp_score_path,"logistic_model_predict_Response/use_filtered_signatures_permutation_and_combination")

# load data ---------------------------------------------------------------
gsva.score <- readr::read_rds(file.path(gsva_score_path,"ICP_GSVA_score_all-possible-features_all-togather.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(-tissue)
exp_ratio <- readr::read_rds(file.path(exp_score_path, "clinical_mean_fold_ratio_features_value.rds.gz"))
exp <- readr::read_tsv(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp/ICPs_FPKM_expression_2.txt")) %>%
  dplyr::select(-gene_id) %>%
  tidyr::gather(-symbol,key="Run",value="exp") %>%
  dplyr::mutate(exp = ifelse(is.na(exp),0,exp)) %>%
  tidyr::spread(key="symbol",value="exp")

gsva.score %>%
  tidyr::unnest()  %>%
  dplyr::inner_join(exp_ratio %>%
                      dplyr::rename("Run"="barcode"),by="Run") %>%
  dplyr::inner_join(exp, by="Run") -> combine_GSVA.exp_ratio

filtered_features <- readr::read_tsv(file.path(immune_res_path,"TCGA_GSVAScore/GSVA_add_exp_ratio/overall_feature_filter.tsv")) %>%
  dplyr::filter(Counts_of_significant_analysis_in_all_cancers>=3) %>%
  .$Features

colnames(combine_GSVA.exp_ratio) <- gsub(colnames(combine_GSVA.exp_ratio),pattern = " ",replacement = "_")
colnames(combine_GSVA.exp_ratio) <- gsub(colnames(combine_GSVA.exp_ratio),pattern = "-",replacement = ".")
combine_GSVA.exp_ratio <- combine_GSVA.exp_ratio[,c("Run",filtered_features)]

# sample data classfication --------------------------------------------------
sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::inner_join(readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","anti_ICB_study_info.txt")) %>%
                      dplyr::select(`data ID`,Author) %>%
                      unique(), by="data ID")

sample_info %>%
  dplyr::select(Run,Cancer.y,blockade,Biopsy_Time,Author,Response) %>%
  unique() -> Run_pubmed.id

combine_GSVA.exp_ratio %>%
  dplyr::inner_join(Run_pubmed.id, by = c("Run")) %>%
  dplyr::select(-Response) %>%
  tidyr::nest(-Cancer.y,-blockade, -Author,.key = "GSVA") -> gsva.score


Run_pubmed.id %>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no"))  %>%
  tidyr::nest(Run,Response,.key = "response") %>%
  dplyr::inner_join(gsva.score, by = c("Cancer.y","blockade","Author")) -> data_for_logistic

Response_statistic <- readr::read_tsv(file.path(gsva_score_path,"Response_statistic_each_dataset.tsv")) %>%
  dplyr::select(-Cancer.y) %>%
  dplyr::mutate(Response = ifelse(is.na(Response),0,Response)) %>%
  dplyr::group_by(Author) %>%
  dplyr::mutate(`non-Response`=sum(`non-Response`),Response=sum(Response)) %>%
  dplyr::select(-Biopsy_Time,-`Response_percentage(%)`) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::mutate(`Response_percentage(%)` = 100*Response/(Response+`non-Response`))
############## fnctions to do logistic regression iteration ##################

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
fn_prediction <- function(group,final_feature){
  
  data_for_logistic  -> .data_for_logistic
  
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
  
  formula.train <- paste("Response~", paste(final_feature$features,collapse = "+"))
  model.train <- glm(as.formula(formula.train), data = data.ready, family = binomial)
  
  # Make predictions
  
  .data_for_logistic %>%
    dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model.train)) %>%
    dplyr::select(-response,-GSVA) %>%
    tidyr::unnest() -> validation_auc
  
  if(min(validation_auc$auc)>=0.6){
    summary(model.train)$coefficients %>%
      as.data.frame() -> coef
    
    # draw picture
    validation_auc %>%
      dplyr::inner_join(Response_statistic,by=c("Author")) %>%
      dplyr::mutate(group = paste(Author,blockade.y,paste("(",Response,"*/",`non-Response`,"**);",sep=""),"AUC","=",signif(auc,2))) %>%
      dplyr::select(roc_data,group,usage) %>%
      tidyr::unnest() -> plot_ready
    plot_ready %>%
      ggplot(aes(x=Specificity,y=Sensitivity)) +
      geom_path(aes(color=group,linetype = usage)) + 
      scale_x_reverse()  +
      my_theme +
      theme(
        legend.position = c(0.75,0.25),
        legend.title = element_blank()
      )
    ggsave(file.path(res_path,paste("ROC_plot",group,"png",sep=".")),device = "png",width = 5, height = 5)
    ggsave(file.path(res_path,paste("ROC_plot",group,"pdf",sep=".")),device = "pdf",width = 5, height = 5)
    sucess="yes"
  } else{
    sucess="no"
  }
  validation_auc %>%
    dplyr::mutate(sucess=sucess) %>%
    dplyr::mutate(AUC_mean = mean(auc)) %>%
    tidyr::nest(-sucess,-AUC_mean,.key="AUC")
}


library(caret)
fn_select_train_test <- function(response,data_spread,percent = 0.7){
  
  set.seed(12356)
  trainIndex <- createDataPartition(response$Response, p = percent, 
                                    list = FALSE, 
                                    times = 1)
  dataTrain <- data_spread[ trainIndex,]
  data.responseTrain <- response[ trainIndex,]
  dataTest  <- data_spread[-trainIndex,]
  data.responseTest <- response[-trainIndex,]
  
  tibble::tibble(Run = c(dataTrain$Run,dataTest$Run), 
                 usage = c(rep("train",length(dataTrain$Run)),rep("test",length(dataTest$Run))))
}

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
  colnames(data.ready) <- gsub("-",".",colnames(data.ready))
  
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
########################################## logistic regression iteration result -----------------------
data_for_logistic %>%
  dplyr::filter(Author=="Riaz" & Biopsy_Time == "on-treatment") %>%
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
  dplyr::mutate(n=purrr::map(response,.f=function(.x){
    nrow(.x)
  })) %>%
  tidyr::unnest(n) %>%
  dplyr::filter(n>5) %>%
  dplyr::select(-n) %>%
  dplyr::mutate(usage = c("validation","validation","validation","validation","validation","train+test","validation")) %>%
  rbind(Riaz.data) %>%
  dplyr::mutate(n = purrr::map2(Author,response,.f=function(.y,.x){
    print(.y)
    .x$Response %>%
      table()%>% 
      as.data.frame() %>%
      tidyr::spread(key=".",value="Freq") %>% 
      dplyr::as.tbl()})) %>%
  tidyr::unnest(n) -> data_for_logistic


############## logistic regression
# stage=c("pre-treatment")

fn_choose_feature <- function(n,f){
  combn(filtered_features,2) %>% 
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    tidyr::gather(key="group",value="features") %>%
    dplyr::mutate(group = paste(group,n,sep="_"))
}

tibble::tibble(n=c(2:8)) %>%
  dplyr::mutate(feature_comn = purrr::map(n,.f=fn_choose_feature,f=filtered_features)) %>%
  tidyr::unnest() -> feature_group
feature_group %>%
  dplyr::select(-n) %>%
  tidyr::nest(-group,.key="final_feature") -> feature_group
library(multidplyr)
cl <- parallel::detectCores()
cluster <- multidplyr::new_cluster(32)

multidplyr::cluster_assign(cluster,"data_for_logistic"=data_for_logistic) 
multidplyr::cluster_assign(cluster,"fn_prediction"=fn_prediction)  
multidplyr::cluster_assign(cluster,"fn_auc_on_validation"=fn_auc_on_validation) 
multidplyr::cluster_library(cluster,"magrittr")
multidplyr::cluster_library(cluster,"pROC") 
multidplyr::cluster_library(cluster,"ggplot2") 
multidplyr::cluster_assign(cluster,"Response_statistic"=Response_statistic)
multidplyr::cluster_assign(cluster,"my_theme"=my_theme)
multidplyr::cluster_assign(cluster,"res_path"=res_path)



feature_group %>%
  dplyr::group_by(group) %>%
  multidplyr::partition(cluster = cluster) %>%
  dplyr::mutate(auc = purrr::map2(group,final_feature, fn_prediction)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  tidyr::unnest(auc) -> feature_group_AUC



feature_group_AUC %>%
  readr::write_rds(file.path(res_path,"AUC_res_all.rds.gz"),compress = "gz")

feature_group_AUC %>%
  dplyr::select(group,AUC_mean,sucess) %>%
  dplyr::filter(sucess == "yes") %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"AUC_res_yes.tsv"))

parallel::stopCluster(cluster)
