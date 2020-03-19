# compare our model perform and other published expression factors
library(magrittr)
library(glmnet)
library(caret)
library(pROC)

# path config -------------------------------------------------------------
basic_path <- file.path("/home/huff/project")
immune_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
res_path <- file.path(immune_path, "ICP_score.new/logistic_model_predict_Response/use_filtered_signatures_permutation_and_combination-from_GSVA_add_exp_ratio_cancerSpecific-evenTop20/select_best_and_compare")
res_path_final <- file.path(immune_path, "ICP_score.new/logistic_model_predict_Response/hill_climbing_191223/repeat_858_best/select_best_and_compare")


# load data ---------------------------------------------------------------

publish_score <- readr::read_rds(file.path("/home/huff/project/immune_checkpoint/result_20171025/ICP_score.new/logistic_model_predict_Response/use_filtered_signatures_permutation_and_combination-from_GSVA_add_exp_ratio_cancerSpecific/select_best_and_compare","CYT_IFN_IMPRESS_score.rds.gz"))
our_score <- readr::read_rds(file.path(res_path_final,"clinical_out_score.rds.gz"))
# our_model <- readr::read_rds(file.path(res_path_final,"logic.model.rds.gz"))

our_score %>%
  dplyr::inner_join(publish_score %>% dplyr::select(-response), by=c("Cancer.y","blockade", "Author")) %>%
  dplyr::mutate(IFN_score = purrr::map2(IFN_score,filter_score, .f=function(.x,.y){
    .x %>%
      dplyr::filter(Run %in% .y$Run)
  })) %>%
  dplyr::mutate(expanded_immune_score = purrr::map2(expanded_immune_score,filter_score, .f=function(.x,.y){
    .x %>%
      dplyr::filter(Run %in% .y$Run)
  })) %>%
  dplyr::mutate(CYT_score = purrr::map2(CYT_score,filter_score, .f=function(.x,.y){
    .x %>%
      dplyr::filter(Run %in% .y$Run)
  })) %>%
  dplyr::mutate(n=purrr::map(response,.f=function(.x){
  nrow(.x)})) %>%
  tidyr::unnest(n) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(-n)  -> combine_data

# functions ---------------------------------------------------------------

fn_select_train_test <- function(response,percent = 0.6){
  response %>%
    dplyr::filter(Response=="yes") -> yes_response
  response %>%
    dplyr::filter(Response=="no") -> no_response
  yes_trainIndex <- createDataPartition(yes_response$Response, p = percent, 
                                        list = FALSE, 
                                        times = 1)
  no_trainIndex <- createDataPartition(no_response$Response, p = percent, 
                                       list = FALSE, 
                                       times = 1)
  
  dataTrain <- response[c(yes_trainIndex,no_trainIndex),]
  
  dataTest  <- response[-c(yes_trainIndex,no_trainIndex),]
  
  tibble::tibble(Run = c(dataTrain$Run,dataTest$Run), 
                 usage = c(rep("train",length(dataTrain$Run)),rep("test",length(dataTest$Run))))
}

fn_auc <- function(data,response,type){
  combine_data %>% 
    dplyr::filter(Author == "Riaz" & Biopsy_Time=="on-treatment" ) %>%
    dplyr::select(response,type) %>%
    dplyr::rename("score"=type) %>%
    dplyr::mutate(combine = purrr::map2(response,score,.f=function(.x,.y){
      .x %>%
        dplyr::inner_join(.y, by="Run")
    })) %>%
    dplyr::select(combine) %>%
    tidyr::unnest() %>%
    dplyr::mutate(Response = as.factor(Response))-> train_set
  train_set %>%
    tidyr::gather(-Run,-Response,key="f",value="score") %>%
    dplyr::filter(is.na(score)) %>%
    dplyr::select(Run) %>%
    unique() -> NA.train
  train_set<-train_set %>% dplyr::filter(! Run %in% NA.train$Run) %>% dplyr::select(-Run)
  colnames(train_set) <- gsub(" ",".",colnames(train_set))
  colnames(train_set) <- gsub("-",".",colnames(train_set))
  data %>%
    dplyr::inner_join(response, by="Run") %>%
    dplyr::mutate(Response = as.factor(Response)) -> test_set
  test_set %>%
    tidyr::gather(-Run,-Response,key="f",value="score") %>%
    dplyr::filter(is.na(score)) %>%
    dplyr::select(Run) %>%
    unique() -> NA.test
  test_set<-test_set %>% dplyr::filter(! Run %in% NA.test$Run) %>% dplyr::select(-Run)
  colnames(test_set) <- gsub(" ",".",colnames(test_set))
  colnames(test_set) <- gsub("-",".",colnames(test_set))
  
  formula <-  paste("Response~", paste(colnames(data)[-1],collapse = "+"))
  # Fit the model
  model <- glm( formula, data = train_set, family = binomial)
  # Make predictions
  probabilities <- model %>% predict(test_set, type = "response")
  predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
  observed.classes <- test_set$Response
  # Model accuracy
  accuracy <- mean(predicted.classes == test_set$Response)
  error <- mean(predicted.classes != test_set$Response)
  
  res.roc <- roc(observed.classes, probabilities,quiet=TRUE)
  auc <- res.roc$auc[[1]]
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


fn_cross_validation <- function(Author, response, filter_score, IFN_score, expanded_immune_score, CYT_score, IMPRESS){
  print(Author)
  # tmp <- tibble::tibble()
  
  # for (i in 1:n) {
    
    # usage <- fn_select_train_test(response)
    
    # auc_ICP <- fn_auc_model(filter_score,  usage, response, our_model)
    
    auc_IFN <- fn_auc(IFN_score, response, type="IFN_score")
    
    auc_EIS <- fn_auc(expanded_immune_score, response, type="expanded_immune_score")
    
    auc_CYT <- fn_auc(CYT_score,  response, type="CYT_score")
    
    # auc_IMPRESS <- fn_auc_IMPRESS(IMPRESS,  usage, response)
    rbind(auc_IFN,auc_EIS) %>%
      rbind(auc_CYT) %>%
      dplyr::mutate(type = c("IFN","EIS","CYT"))
  # }
  # tmp
  # tmp %>%
  #   tidyr::gather(key="feature_group",value="auc") %>%
  #   dplyr::group_by(feature_group) %>%
  #   dplyr::mutate(mean_auc = mean(auc)) 
}


# calculation -------------------------------------------------------------

combine_data %>%
  # dplyr::filter(Author!="Riaz") %>%
  dplyr::mutate(AUC = purrr::pmap(list(Author,response, filter_score, IFN_score, expanded_immune_score, CYT_score,IMPRESS),fn_cross_validation)) %>%
  dplyr::select(-response, -filter_score, -IFN_score, -expanded_immune_score, -CYT_score, -IMPRESS) %>%
  tidyr::unnest()  -> other_AUC


# plot --------------------------------------------------------------------
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
Response_statistic <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/result_20171025","ICP_score/5.GSVA-ICPs_exp_site-all_possible","Response_statistic_each_dataset.tsv")) %>%
  dplyr::select(-Cancer.y) %>%
  dplyr::mutate(Response = ifelse(is.na(Response),0,Response)) %>%
  dplyr::group_by(Author,Biopsy_Time) %>%
  dplyr::mutate(`non-Response`=sum(`non-Response`),Response=sum(Response)) %>%
  dplyr::select(-`Response_percentage(%)`) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::mutate(`Response_percentage(%)` = 100*Response/(Response+`non-Response`))

other_AUC %>%
  tidyr::unnest() %>%
  dplyr::inner_join(Response_statistic,by=c("Author","Biopsy_Time")) %>%
  dplyr::mutate(Biopsy_Time.short=strsplit(Biopsy_Time,split = "-")[[1]][1]) %>%
  dplyr::mutate(group = paste(Author,Biopsy_Time.short, blockade.y,paste("(",Response,"*/",`non-Response`,"**);",sep=""))) %>%
  dplyr::mutate(`1-Specificity`=1-Specificity) %>%
  dplyr::filter(best=="best") %>%
  dplyr::select(auc,type,Sensitivity, `1-Specificity`,group) -> plot_text

other_AUC %>%
  tidyr::unnest() %>%
  dplyr::inner_join(Response_statistic,by=c("Author","Biopsy_Time")) %>%
  dplyr::mutate(Biopsy_Time.short=strsplit(Biopsy_Time,split = "-")[[1]][1]) %>%
  dplyr::mutate(group = paste(Author,Biopsy_Time.short, blockade.y,paste("(",Response,"*/",`non-Response`,"**);",sep=""))) %>%
  dplyr::mutate(`1-Specificity`=1-Specificity) %>%
  ggplot(aes(x=`1-Specificity`,y=Sensitivity)) +
  geom_path(aes(color=group)) + 
  # scale_x_reverse()  +
  facet_grid(.~type) +
  geom_text(aes(x=`1-Specificity`,y =Sensitivity+0.02, label=signif(auc,2), color=group), data=plot_text) +
  my_theme +
  theme(
    # legend.position = c(0.75,0.25),
    legend.title = element_blank()
  )
ggsave(file.path(immune_path, "ICP_score.new/publish_model_roc","publish_model_roc.pdf"),width = 12,height = 6)
