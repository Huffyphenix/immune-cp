# compare our model perform and other published expression factors
library(magrittr)
library(glmnet)
library(caret)
library(pROC)

# path config -------------------------------------------------------------
basic_path <- file.path("/home/huff/project")
immune_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
res_path <- file.path(immune_path, "ICP_score.new/logistic_model_predict_Response/use_filtered_signatures_permutation_and_combination-from_GSVA_add_exp_ratio_cancerSpecific/select_best_and_compare")


# load data ---------------------------------------------------------------

publish_score <- readr::read_rds(file.path(res_path,"CYT_IFN_score.rds.gz"))
our_score <- readr::read_rds(file.path(res_path,"clinical_out_score.rds.gz"))
our_model <- readr::read_rds(file.path(res_path,"our_ICP_lm_model.rds.gz"))

our_score %>%
  dplyr::filter(!usage %in% c("train","test")) %>%
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
  dplyr::filter(Author != "Auslander") %>%
  dplyr::select(-usage)-> combine_data


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

fn_auc <- function(data, usage,response){
  data %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(usage, by="Run") %>%
    dplyr::inner_join(response, by="Run") %>%
    dplyr::mutate(Response = as.factor(Response)) %>%
    dplyr::filter(usage == "train" ) %>%
    dplyr::select(-usage) -> train_set
  train_set %>%
    tidyr::gather(-Run,-Response,key="f",value="score") %>%
    dplyr::filter(is.na(score)) %>%
    dplyr::select(Run) %>%
    unique() -> NA.train
  train_set<-train_set %>% dplyr::filter(! Run %in% NA.train$Run) %>% dplyr::select(-Run)
  colnames(train_set) <- gsub(" ",".",colnames(train_set))
  colnames(train_set) <- gsub("-",".",colnames(train_set))
  data %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(usage, by="Run") %>%
    dplyr::inner_join(response, by="Run") %>%
    dplyr::mutate(Response = as.factor(Response)) %>%
    dplyr::filter(usage == "test" ) %>%
    dplyr::select(-usage) -> test_set
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
  auc
}

fn_auc_model <- function(data, usage,response, model){
  data %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(usage, by="Run") %>%
    dplyr::inner_join(response, by="Run") %>%
    dplyr::mutate(Response = as.factor(Response)) %>%
    dplyr::filter(usage == "train" ) %>%
    dplyr::select(-usage) -> train_set
  train_set %>%
    tidyr::gather(-Run,-Response,key="f",value="score") %>%
    dplyr::filter(is.na(score)) %>%
    dplyr::select(Run) %>%
    unique() -> NA.train
  train_set<-train_set %>% dplyr::filter(! Run %in% NA.train$Run) %>% dplyr::select(-Run)
  colnames(train_set) <- gsub(" ",".",colnames(train_set))
  colnames(train_set) <- gsub("-",".",colnames(train_set))
  data %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(usage, by="Run") %>%
    dplyr::inner_join(response, by="Run") %>%
    dplyr::mutate(Response = as.factor(Response)) %>%
    dplyr::filter(usage == "test" ) %>%
    dplyr::select(-usage) -> test_set
  test_set %>%
    tidyr::gather(-Run,-Response,key="f",value="score") %>%
    dplyr::filter(is.na(score)) %>%
    dplyr::select(Run) %>%
    unique() -> NA.test
  test_set<-test_set %>% dplyr::filter(! Run %in% NA.test$Run) %>% dplyr::select(-Run)
  colnames(test_set) <- gsub(" ",".",colnames(test_set))
  colnames(test_set) <- gsub("-",".",colnames(test_set))
  
  # Make predictions
  probabilities <- model %>% predict(test_set, type = "response")
  predicted.classes <- ifelse(probabilities > 0.5, "yes", "no")
  observed.classes <- test_set$Response
  # Model accuracy
  accuracy <- mean(predicted.classes == test_set$Response)
  error <- mean(predicted.classes != test_set$Response)
  
  res.roc <- roc(observed.classes, probabilities,quiet=TRUE)
  auc <- res.roc$auc[[1]]
  auc
}
fn_cross_validation <- function(n,Author, response, filter_score, IFN_score, expanded_immune_score, CYT_score){
  print(Author)
  tmp <- tibble::tibble()
  
  for (i in 1:n) {
    
    usage <- fn_select_train_test(response)
    
    auc_ICP <- fn_auc_model(filter_score,  usage, response, our_model)
    
    auc_IFN <- fn_auc(IFN_score,  usage, response)
    
    auc_EIS <- fn_auc(expanded_immune_score,  usage, response)
    
    auc_CYT <- fn_auc(CYT_score,  usage, response)
    
    tmp<- rbind(tmp,tibble::tibble(auc_ICP=auc_ICP,auc_IFN=auc_IFN,auc_EIS=auc_EIS,auc_CYT=auc_CYT))
  }
  tmp %>%
    tidyr::gather(key="feature_group",value="auc") %>%
    dplyr::group_by(feature_group) %>%
    dplyr::mutate(mean_auc = mean(auc)) 
}

fn_cross_validation_VA <- function(n,Author, response, filter_score, IFN_score, expanded_immune_score, CYT_score){
  print(Author)
  tmp <- tibble::tibble()
  
  for (i in 1:n) {
    
    usage <- fn_select_train_test(response)
    
    auc_ICP <- fn_auc(filter_score,  usage, response)
    
    auc_IFN <- fn_auc(IFN_score,  usage, response)
    
    auc_EIS <- fn_auc(expanded_immune_score,  usage, response)
    
    auc_CYT <- fn_auc(CYT_score,  usage, response)
    
    tmp<- rbind(tmp,tibble::tibble(auc_ICP=auc_ICP,auc_IFN=auc_IFN,auc_EIS=auc_EIS,auc_CYT=auc_CYT))
  }
  tmp %>%
    tidyr::gather(key="feature_group",value="auc") %>%
    dplyr::group_by(feature_group) %>%
    dplyr::mutate(mean_auc = mean(auc)) 
}
# Run ---------------------------------------------------------------------
res_VanAllen <- tibble::tibble()
for (j in 1:100) {
  combine_data %>%
    dplyr::filter(Author == "Van Allen") %>%
    dplyr::mutate(AUC = purrr::pmap(list(Author,response, filter_score, IFN_score, expanded_immune_score, CYT_score),fn_cross_validation_VA,n=5)) %>%
    dplyr::select(-response, -filter_score, -IFN_score, -expanded_immune_score, -CYT_score) %>%
    tidyr::unnest() %>%
    dplyr::select(-auc) %>%
    unique() %>%
    dplyr::mutate(repeats = j) -> .tmp
  rbind(res_VanAllen,.tmp) -> res_VanAllen
}
res_VanAllen %>%
  readr::write_tsv(file.path(res_path,"VanAllen_dataset_compare_res.tsv"))

res2 <- tibble::tibble()
for (j in 1:100) {
  combine_data %>%
    dplyr::mutate(AUC = purrr::pmap(list(Author,response, filter_score, IFN_score, expanded_immune_score, CYT_score),fn_cross_validation,n=5)) %>%
    dplyr::select(-response, -filter_score, -IFN_score, -expanded_immune_score, -CYT_score) %>%
    tidyr::unnest() %>%
    dplyr::select(-auc) %>%
    unique() %>%
    dplyr::mutate(repeats = j) -> .tmp
  rbind(res2,.tmp) -> res2
}
res2 %>%
  readr::write_tsv(file.path(res_path,"All_dataset_compare_res.tsv"))
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
tibble::tibble(feature_group = unique(res2$feature_group),
               feature_group2 = c("Immune regulators","IFNy","Expanded immune","CYT")) -> feature_correspond
res2 %>%
  dplyr::inner_join(feature_correspond,by="feature_group") %>%
  dplyr::mutate(group = paste(Author,Biopsy_Time,blockade,sep=",")) %>%
  ggplot(aes(x=feature_group2,y=mean_auc)) +
  geom_boxplot(aes(fill = group)) +
  labs(y="Mean AUC of 5 fold cross validation\n(100 repetitions)") +
  my_theme +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    legend.key = element_rect(colour="white")
  )

ggsave(file.path(res_path,"All_compare_plot.png"),device = "png",height = 4, width = 8)
ggsave(file.path(res_path,"All_compare_plot.pdf"),device = "pdf",height = 4, width = 8)

comp_list <- list(c("CYT","Immune regulators"),
                  c("Expanded immune","Immune regulators"),
                  c("IFNy","Immune regulators"))
res_VanAllen %>%
  dplyr::inner_join(feature_correspond,by="feature_group") %>%
  dplyr::mutate(group = paste(Author,Biopsy_Time,blockade,sep=",")) %>%
  ggplot(aes(x=feature_group2,y=mean_auc)) +
  geom_boxplot(aes(fill = group)) +
  labs(y="Mean AUC of 5 fold cross validation\n(100 repetition)") +
  my_theme+
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "top",
    legend.key = element_rect(colour="white"),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  ) +
  ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")

ggsave(file.path(res_path,"VanAllen_compare_plot.png"),device = "png",height = 4, width = 4)
ggsave(file.path(res_path,"VanAllen_compare_plot.pdf"),device = "pdf",height = 4, width = 4)

# res2 %>% dplyr::filter(feature_group=="auc_ICP") %>% dplyr::mutate(res="2") %>%
#   rbind(res %>% dplyr::filter(feature_group=="auc_ICP") %>% dplyr::mutate(res="1")) %>%
#   dplyr::mutate(group = paste(Author,Biopsy_Time,blockade,sep=",")) %>%
#   ggplot(aes(x=feature_group,y=mean_auc)) +
#   geom_boxplot(aes(fill = res)) 
