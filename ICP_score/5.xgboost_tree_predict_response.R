################################
# construct logistics model for GSVA score of all possible features
# GSVA was generate in all samples, all datasets as one
# for both pre- and on-treatment samples prediction
# train data is Riaz, on-treatment
################################
library(magrittr)
library(tidyverse)
library(xgboost)
library(caret)

# data path ---------------------------------------------------------------
# server 1 basic path
basic_path <- file.path("/home/huff/project")

immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path("/data/TCGA/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
exp_score_path <- file.path(immune_res_path,"ICP_score.new")
gsva_score_path <-  file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")

res_path <- file.path(exp_score_path,"xgboost")

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
  dplyr::inner_join(exp, by="Run")-> combine_GSVA.exp_ratio

# sample data classfication --------------------------------------------------
sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::inner_join(readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","anti_ICB_study_info.txt")) %>%
                      dplyr::select(`data ID`,Author) %>%
                      unique(), by="data ID")

sample_info %>%
  dplyr::select(Run,Cancer.y,blockade,blockade,Biopsy_Time,Author,Response) %>%
  unique() -> Run_pubmed.id

combine_GSVA.exp_ratio %>%
  dplyr::inner_join(Run_pubmed.id, by = c("Run")) %>%
  dplyr::select(-Response) %>%
  tidyr::nest(-Cancer.y,-blockade, -Biopsy_Time, -Author,.key = "GSVA") -> gsva.score


Run_pubmed.id %>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), 1, 0))  %>%
  tidyr::nest(Run,Response,.key = "response") %>%
  dplyr::inner_join(gsva.score, by = c("Cancer.y","blockade","Biopsy_Time","Author")) -> data_for_logistic


# test xgboost in one of our dataset --------------------------------------------------

data_for_logistic$response[[1]] %>%
  dplyr::select(-Run)-> data.response
data_for_logistic$GSVA[[1]] %>%
  dplyr::select(-Run)-> data

set.seed(3456)
trainIndex <- createDataPartition(data.response$Response, p = .8, 
                                  list = FALSE, 
                                  times = 1)
head(trainIndex)
dataTrain <- data[ trainIndex,]
data.responseTrain <- data.response[ trainIndex,]
dataTest  <- data[-trainIndex,]
data.responseTest <- data.response[-trainIndex,]

Train_set <- xgb.DMatrix(data = as.matrix(dataTrain), label = data.responseTrain$Response)
Test_set <- xgb.DMatrix(data = as.matrix(dataTest), label = data.responseTest$Response)

# use xgboost
xgb_model <- xgboost(data = Train_set,max.depth = 2, eta = 1, nthread = 2, nrounds = 2, 
        objective = "binary:logistic", verbose = 2)
pred <- predict(xgb_model, as.matrix(dataTest))

prediction <- as.numeric(pred > 0.5)
print(head(prediction))
err <- mean(as.numeric(pred > 0.5) != data.responseTest$Response)
print(paste("test-error=", err))
res.roc <- roc(data.responseTest$Response, pred,quiet=TRUE)
auc <- res.roc$auc[[1]]

# use xgb.train
watchlist <- list(train=Train_set, test=Test_set)
bst <- xgb.train(data=Train_set, max.depth=10, eta=1, nthread = 2, nrounds=2, watchlist=watchlist, objective = "binary:logistic")
bst.pred <- predict(bst,as.matrix(dataTest))
bst.res.roc <- roc(data.responseTest$Response,bst.pred,quiet = T)
bst.res.roc$auc


# use function to serve each dataset as train data ------------------------

fn_xgboost <- function(i){
  data_for_logistic$response[[i]] -> response
  data_for_logistic$GSVA[[i]] %>%
    dplyr::inner_join(response %>% dplyr::select(Run),by="Run") %>%
    dplyr::select(-Run)-> data
  
  set.seed(3456)
  trainIndex <- createDataPartition(data.response$Response, p = .8, 
                                    list = FALSE, 
                                    times = 1)
  
  # data prepare
  trainIndex <- createDataPartition(response$Response, p = .8, 
                                    list = FALSE, 
                                    times = 1)
  head(trainIndex)
  dataTrain <- data[trainIndex,]
  data.responseTrain <- response[trainIndex,]
  dataTest  <- data[-trainIndex,]
  data.responseTest <- response[-trainIndex,]
  
  Train_set <- xgb.DMatrix(data = as.matrix(dataTrain), label = data.responseTrain$Response)
  Test_set <- xgb.DMatrix(data = as.matrix(dataTest), label = data.responseTest$Response)
  
  # use xgboost
  xgb_model <- xgboost(data = Train_set,max.depth = 2, eta = 1, nthread = 2, nrounds = 2, 
                       objective = "binary:logistic", verbose = 2)
  pred <- predict(xgb_model, as.matrix(dataTest))
  
  prediction <- as.numeric(pred > 0.5)
  print(head(prediction))
  err <- mean(as.numeric(pred > 0.5) != data.responseTest$Response)
  print(paste("test-error=", err))
  res.roc <- roc(data.responseTest$Response, pred,quiet=TRUE)
  auc <- res.roc$auc[[1]]
  
  # get AUC on validation set
  data_for_logistic[-c(3,i),] %>%
    dplyr::mutate(AUC = purrr::map2(GSVA,response,fn_AUC,model=xgb_model)) %>%
    dplyr::select(-response,-GSVA) %>%
    tidyr::unnest() %>%
    rbind(tibble::tibble(Cancer.y=data_for_logistic$Cancer.y[[i]],
                  blockade="train_test",
                  Biopsy_Time="train_test",
                  Author=data_for_logistic$Author[[i]],
                  AUC=auc))
}

fn_xgb.train <- function(i){
  data_for_logistic$response[[i]] -> response
  data_for_logistic$GSVA[[i]] %>%
    dplyr::inner_join(response %>% dplyr::select(Run),by="Run") %>%
    dplyr::select(-Run)-> data
  
  set.seed(3456)
  trainIndex <- createDataPartition(data.response$Response, p = .8, 
                                    list = FALSE, 
                                    times = 1)
  
  # data prepare
  trainIndex <- createDataPartition(response$Response, p = .8, 
                                    list = FALSE, 
                                    times = 1)
  head(trainIndex)
  dataTrain <- data[trainIndex,]
  data.responseTrain <- response[trainIndex,]
  dataTest  <- data[-trainIndex,]
  data.responseTest <- response[-trainIndex,]
  
  Train_set <- xgb.DMatrix(data = as.matrix(dataTrain), label = data.responseTrain$Response)
  Test_set <- xgb.DMatrix(data = as.matrix(dataTest), label = data.responseTest$Response)
  watchlist <- list(train=Train_set, test=Test_set)
  
 
  # use xgb.train
  xgb_model <- xgb.train(data=Train_set, max.depth=2, eta=1, nthread = 2, nrounds=2, watchlist=watchlist, objective = "binary:logistic")
  pred <- predict(xgb_model, as.matrix(dataTest))
  
  prediction <- as.numeric(pred > 0.5)
  print(head(prediction))
  err <- mean(as.numeric(pred > 0.5) != data.responseTest$Response)
  print(paste("test-error=", err))
  res.roc <- roc(data.responseTest$Response, pred,quiet=TRUE)
  auc <- res.roc$auc[[1]]
  
  # get AUC on validation set
  data_for_logistic[-c(3,i),] %>%
    dplyr::mutate(AUC = purrr::map2(GSVA,response,fn_AUC,model=xgb_model)) %>%
    dplyr::select(-response,-GSVA) %>%
    tidyr::unnest() %>%
    rbind(tibble::tibble(Cancer.y=data_for_logistic$Cancer.y[[i]],
                         blockade="train_test",
                         Biopsy_Time="train_test",
                         Author=data_for_logistic$Author[[i]],
                         AUC=auc))
}
fn_AUC <- function(data,response,model){
  data %>%
    dplyr::inner_join(response %>% dplyr::select(Run),by="Run") %>%
    dplyr::select(-Run) -> data
  pred <- predict(model, as.matrix(data))
  
  prediction <- as.numeric(pred > 0.5)
  print(head(prediction))
  err <- mean(as.numeric(pred > 0.5) != response$Response)
  print(paste("test-error=", err))
  res.roc <- roc(response$Response, pred,quiet=TRUE)
  auc <- res.roc$auc[[1]]
  
  auc
}

# calculation -----------------------
# xgboost
tibble::tibble(i=c(1:8))[-c(3,4,5),] %>%
  dplyr::mutate(res = purrr::map(i,fn_xgboost)) -> res # non significance
res.1 %>%
  tidyr::unnest() %>%
  view()

# xgb.train
tibble::tibble(i=c(1:8))[-c(3,4,5),] %>%
  dplyr::mutate(res = purrr::map(i,fn_xgb.train)) -> res.1# non significance

res.1 %>%
  tidyr::unnest() %>%
  view()
