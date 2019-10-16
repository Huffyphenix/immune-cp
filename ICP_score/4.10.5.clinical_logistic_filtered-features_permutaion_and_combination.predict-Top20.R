################################
# construct logistics model for GSVA score of all possible features
# GSVA was generate in each dataset
# train data for all data, is Riaz
# seperate the pre- and on-treatment data
# features are permutation and combination of 44 features filtered from TCGA

################################
start <- Sys.time()
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
# res_path <- file.path(exp_score_path,"logistic_model_predict_Response/use_filtered_signatures_permutation_and_combination")
res_path <- file.path(exp_score_path,"logistic_model_predict_Response/use_filtered_signatures_permutation_and_combination-from_GSVA_add_exp_ratio_cancerSpecific-evenTop20")
# res_path <- file.path(exp_score_path,"logistic_model_predict_Response/test")

# load data ---------------------------------------------------------------
gsva.score <- readr::read_rds(file.path(gsva_score_path,"new-ICP_GSVA_score_all-possible-features_all-togather.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::select(-tissue)
exp_ratio <- readr::read_rds(file.path(exp_score_path, "new-clinical_mean_fold_ratio_features_value.rds.gz"))
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

# filtered_features <- readr::read_tsv(file.path(immune_res_path,"TCGA_GSVAScore/GSVA_add_exp_ratio/overall_feature_filter.tsv")) %>%
#   dplyr::filter(Counts_of_significant_analysis_in_all_cancers>=3) %>%
#   .$Features
filtered_features <- readr::read_tsv(file.path("/home/huff/project/immune_checkpoint/result_20171025/TCGA_GSVAScore/GSVA_add_exp_ratio_cancerSpecific-new/top_50_Feature_filtered_by_cor_FC-noMean.tsv")) %>%
  .$Features
colnames(combine_GSVA.exp_ratio) <- gsub(colnames(combine_GSVA.exp_ratio),pattern = " ",replacement = "_")
colnames(combine_GSVA.exp_ratio) <- gsub(colnames(combine_GSVA.exp_ratio),pattern = "-",replacement = ".")
filtered_features <- intersect(colnames(combine_GSVA.exp_ratio),filtered_features)
combine_GSVA.exp_ratio <- combine_GSVA.exp_ratio[,c("Run",filtered_features)]

# sample data classfication --------------------------------------------------
sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::inner_join(readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","anti_ICB_study_info.txt")) %>%
                      dplyr::select(`data ID`,Author) %>%
                      unique(), by="data ID")

sample_info %>%
  dplyr::select(Run,Cancer.y,blockade,Biopsy_Time,Author,Response, Response_standard,Second_Response_standard) %>%
  unique() -> Run_pubmed.id

combine_GSVA.exp_ratio %>%
  dplyr::inner_join(Run_pubmed.id, by = c("Run")) %>%
  dplyr::select(-Response) %>%
  tidyr::nest(-Cancer.y,-blockade, -Author,-Biopsy_Time,.key = "GSVA") -> gsva.score


Run_pubmed.id %>%
  dplyr::filter(! Response %in% c("NE")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no"))  %>%
  dplyr::mutate(Response = ifelse(blockade == "anti-CTLA-4" & Second_Response_standard %in% c("R"), "yes", Response)) %>%
  dplyr::mutate(Response = ifelse(blockade == "anti-CTLA-4" & Second_Response_standard %in% c("long-survival","NR"), "no", Response)) %>%
  dplyr::mutate(x = ifelse(is.na(Second_Response_standard) & blockade == "anti-CTLA-4", "1","2")) %>%
  dplyr::filter(x=="2") %>%
  dplyr::select(-Response_standard,-Second_Response_standard,-x) %>%
  tidyr::nest(Run,Response,.key = "response") %>%
  dplyr::inner_join(gsva.score, by = c("Cancer.y","blockade","Author","Biopsy_Time")) -> data_for_logistic

# Response_statistic <- readr::read_tsv(file.path(gsva_score_path,"Response_statistic_each_dataset.tsv")) %>%
#   dplyr::select(-Cancer.y) %>%
#   dplyr::mutate(Response = ifelse(is.na(Response),0,Response)) %>%
#   dplyr::group_by(Author,Biopsy_Time) %>%
#   dplyr::mutate(`non-Response`=sum(`non-Response`),Response=sum(Response)) %>%
#   dplyr::select(-`Response_percentage(%)`) %>%
#   dplyr::ungroup() %>%
#   unique() %>%
#   dplyr::mutate(`Response_percentage(%)` = 100*Response/(Response+`non-Response`))
############## fnctions to do logistic regression iteration ##################
print("######################### fnctions to do logistic regression")
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
        # dplyr::select(-Run) %>%
        dplyr::mutate(Response = as.factor(Response))
    })) %>%
    dplyr::select(data.ready) %>%
    tidyr::unnest() -> data.ready
  # data.ready <- na.omit(data.ready)
  data.ready %>%
    tidyr::gather(-Run,-Response,key="f",value="score") %>%
    dplyr::filter(is.na(score)) %>%
    dplyr::select(Run) %>%
    unique() -> NA.train
  data.ready<-data.ready %>% dplyr::filter(! Run %in% NA.train$Run) %>% dplyr::select(-Run)
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
      dplyr::inner_join(Response_statistic,by=c("Author","Biopsy_Time")) %>%
      dplyr::mutate(group = paste(Author,blockade.y,Biopsy_Time,paste("(",Response,"*/",`non-Response`,"**);",sep=""),"AUC","=",signif(auc,2),usage)) %>%
      dplyr::select(roc_data,group,usage) %>%
      tidyr::unnest() -> plot_ready
    plot_ready %>%
      ggplot(aes(x=Specificity,y=Sensitivity)) +
      geom_path(aes(color=group)) + 
      scale_x_reverse()  +
      my_theme +
      theme(
        legend.position = c(0.72,0.17),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.key.height = unit(0.15,"inches"),
        legend.key.width = unit(0.15,"inches"),
        legend.key = element_rect(colour="white")
      )
    ggsave(file.path(res_path,paste("ROC_plot",group,"png",sep=".")),device = "png",width = 8, height = 5)
    ggsave(file.path(res_path,paste("ROC_plot",group,"pdf",sep=".")),device = "pdf",width = 8, height = 5)
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
  # tibble::tibble(roc = res.roc)
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

############## choose form top50 features----------------
# get model from 70% train and use it on 30% Test and validation data
print("########################### get top20 features")
data_for_logistic %>%
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

formula.train <- paste("Response~", paste(filtered_features,collapse = "+"))
model.train <- glm(as.formula(formula.train), data = data.ready, family = binomial)

  model.train$coefficients %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    dplyr::mutate(names = names(model.train$coefficients)) %>%
    dplyr::filter(names!="(Intercept)") %>%
    dplyr::arrange(desc(abs(.))) %>%
    head(20) %>%
    .$names -> top_20features

  model.train$coefficients %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    dplyr::mutate(names = names(model.train$coefficients)) %>%
    dplyr::filter(names!="(Intercept)") %>%
    dplyr::arrange(desc(abs(.))) %>%
    head(20) %>%
  readr::write_tsv(file.path(res_path,"AUC_top_20features.tsv"))
############## logistic regression---------------------
# stage=c("pre-treatment")
print("####################### get top20 features combination")
  
fn_choose_feature <- function(n,f){
  combn(f,n) %>% 
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    tidyr::gather(key="group",value="features") %>%
    dplyr::mutate(group = paste(group,n,sep="_"))
}

tibble::tibble(n=c(2:10)) %>%
  dplyr::mutate(feature_comn = purrr::map(n,.f=fn_choose_feature,f=top_20features)) %>%
  tidyr::unnest() -> feature_group
feature_group %>%
  dplyr::select(-n) %>%
  tidyr::nest(-group,.key="final_feature") -> feature_group

print("########################## parallel anasis")

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


print("###################### parallel analysis finished")

print("###################### output results")

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
  dplyr::filter(sucess == "yes") %>%
  dplyr::select(group,AUC_mean,sucess) %>%
  tidyr::unnest() %>%
  readr::write_tsv(file.path(res_path,"AUC_res_yes.tsv"))

feature_group_AUC %>%
  dplyr::filter(sucess == "yes") %>%
  dplyr::select(group,AUC_mean,sucess,AUC) %>%
  tidyr::unnest() %>%
  dplyr::group_by(group) %>%
  dplyr::filter(!usage %in% c("train+test","train")) %>%
  dplyr::mutate(max_auc = max(auc), min_auc = min(auc), mid_auc = quantile(auc,0.5), mean_auc = mean(auc)) %>%
  dplyr::select(group,max_auc, min_auc, mid_auc, mean_auc) %>%
  dplyr::arrange(desc(mean_auc, mid_auc,max_auc, min_auc)) %>%
  unique() %>% 
  readr::write_tsv(file.path(res_path,"AUC_res_yes_rank.tsv"))

print("###################### all analysis finished")

############## validate the results
# AUC_res_yes_rank
# # The first ine : V483182_5
# feature_group_AUC %>%
#   dplyr::filter(group %in% c("V483182_5")) -> final_res
# 
# final_res$final_feature[[1]]$features -> final_features
final_features <- c("Immune_and_tumor_cell_almost","Fold.All_gene.PDCD1","Fold.Pair_14.BTN2A1","Fold.Pair_3.ICOSLG","Fold.TwoSide.PVR")

data_for_logistic %>%
  dplyr::mutate(filter_score = purrr::map(GSVA, .f =function(.x){
    .x[,c("Run",final_features)]
  })) %>%
  dplyr::select(-GSVA) %>%
  readr::write_rds(file.path(res_path,"select_best_and_compare/clinical_out_score.rds.gz"))
# 
# # get model from 70% train and use it on 30% Test and validation data
data_for_logistic %>%
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
# 
formula.train <- paste("Response~", paste(final_features,collapse = "+"))
model.train <- glm(as.formula(formula.train), data = data.ready, family = binomial)
model.train %>%
  readr::write_rds(file.path(file.path(res_path,"select_best_and_compare/our_ICP_lm_model.rds.gz")))
# # Make predictions
# 
data_for_logistic %>%
  dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model.train)) %>%
  dplyr::select(-response,-GSVA) %>%
  tidyr::unnest() -> validation_auc

# draw picture

validation_auc %>%
  dplyr::inner_join(Response_statistic,by=c("Author","Biopsy_Time")) %>%
  dplyr::mutate(group = paste(Author,blockade.y,Biopsy_Time,paste("(",yes,"*/",`no`,"**);",sep=""),"AUC","=",signif(auc,2),usage)) %>%
  dplyr::select(roc_data,group,usage) %>%
  tidyr::unnest() -> plot_ready
plot_ready %>%
  ggplot(aes(x=Specificity,y=Sensitivity)) +
  geom_path(aes(color=group)) + 
  scale_x_reverse()  +
  my_theme +
  theme(
    legend.position = c(0.72,0.17),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches"),
    legend.key = element_rect(colour="white")
  )
ggsave(file.path(res_path,paste("ROC_plot","BEST","png",sep=".")),device = "png",width = 8, height = 5)
ggsave(file.path(res_path,paste("ROC_plot","BEST","pdf",sep=".")),device = "pdf",width = 8, height = 5)

start
Sys.time()
parallel::stopCluster(cluster)
