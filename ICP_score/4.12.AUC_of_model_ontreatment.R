################################
# get the AUC result of the model we get from ICP_score/4.10.5.clinical_logistic_use_all-features.predict.Res-nonRes-ontreatment.R
# dealing the result of logistics model for GSVA score of all possible features, only on-treatment
# GSVA was generate in each dataset
# for only on-treatment samples prediction results
# train data for on-treatment data, is Riaz, on-treatment
# validation data is the on-treat data
# there is no significance in other pre-treat data
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
res_path <- file.path(exp_score_path,"logistic_model_predict_Response/on_treatment")

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

# sample data classfication --------------------------------------------------
sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::inner_join(readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","anti_ICB_study_info.txt")) %>%
                      dplyr::select(`data ID`,Author) %>%
                      unique(), by="data ID") %>%
  dplyr::mutate(blockade=ifelse(blockade=="anti–PD-1、CTLA-4","anti–PD-1/CTLA-4",blockade))

sample_info %>%
  dplyr::select(Run,Cancer.y,blockade,blockade,Biopsy_Time,Author,Response) %>%
  unique() -> Run_pubmed.id

combine_GSVA.exp_ratio %>%
  dplyr::inner_join(Run_pubmed.id, by = c("Run")) %>%
  dplyr::select(-Response) %>%
  tidyr::nest(-Cancer.y,-blockade, -Biopsy_Time, -Author,.key = "GSVA") -> gsva.score


Run_pubmed.id %>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "yes", "no"))  %>%
  tidyr::nest(Run,Response,.key = "response") %>%
  dplyr::inner_join(gsva.score, by = c("Cancer.y","blockade","Biopsy_Time","Author")) -> data_for_logistic


# load model  -------------------------------------------------------------

model <- readr::read_rds(file.path(res_path,"Final_model_3_.rds.gz"))


# function ----------------------------------------------------------------

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
# data process ------------------------------------------------------------

data_for_logistic[-3,] %>%
  dplyr::filter(Biopsy_Time == "on-treatment") %>%
  dplyr::mutate(sample_group = purrr::map2(response,GSVA,fn_select_train_test,percent = 0.7)) -> data_for_logistic


data_for_logistic %>%
  dplyr::filter(Author == "Riaz") %>%
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
  dplyr::select(-sample_group) %>%
  dplyr::mutate(usage = c("validation","train")) %>%
  rbind(Riaz.data) %>%
  dplyr::mutate(n = purrr::map2(Author,response,.f=function(.y,.x){
    print(.y)
    .x$Response %>%
      table()%>% 
      as.data.frame() %>%
      tidyr::spread(key=".",value="Freq") %>% 
      dplyr::as.tbl()})) %>%
  tidyr::unnest(n) -> data_for_logistic


# model prediction --------------------------------------------------------

data_for_logistic %>%
  dplyr::filter(usage != "train+test") %>%
  dplyr::mutate(auc = purrr::map2(response,GSVA,fn_auc_on_validation,model=model)) %>%
  dplyr::select(-response,-GSVA) %>%
  tidyr::unnest() -> validation_auc


# AUC plot ----------------------------------------------------------------

validation_auc %>%
  dplyr::mutate(group = paste(Author,blockade,paste("(",yes,"*/",no,"**); ",sep=""),"AUC"," = ",signif(auc,2),"; ",usage,sep="")) %>%
  dplyr::select(roc_data,group,usage) %>%
  tidyr::unnest() -> plot_ready.auc

tibble::tibble(group=c("Base_line","Base_line"),usage=c("Base_line","Base_line"),
               Sensitivity=c(0,1),Specificity=c(1,0),thresholds=c(0,1), best=c("no","no")) -> base_line
plot_ready.auc %>%
  rbind(base_line) -> plot_ready
plot_ready <- within(plot_ready,group<- factor(group,levels = c("Auslanderanti–PD-1/CTLA-4(1*/22**); AUC = 0.91; validation",
                                                  "Riazanti–PD-1(12*/29**); AUC = 0.83; train",
                                                  "Riazanti–PD-1(6*/10**); AUC = 0.9; test",
                                                  "Base_line")))
plot_ready %>%
  ggplot(aes(x=Specificity,y=Sensitivity)) +
  geom_path(aes(color=group)) + 
  scale_x_reverse()  +
  scale_color_manual(values = c(c("#CD5555", "#1C86EE", "#EE7621", "#030303"))) +
  my_theme +
  theme(
    legend.position = c(0.64,0.15),
    legend.title = element_blank(),
    legend.text = element_text(size=8),
    legend.key.height = unit(0.15,"inches"),
    legend.key.width = unit(0.15,"inches")
  )
ggsave(file.path(res_path,"Best_model.AUC.png"),device = "png",height = 4, width = 6)
ggsave(file.path(res_path,"Best_model.AUC.pdf"),device = "pdf",height = 4, width = 6)
