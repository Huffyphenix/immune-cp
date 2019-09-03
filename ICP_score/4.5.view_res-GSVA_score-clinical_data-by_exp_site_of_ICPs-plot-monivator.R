library(ggplot2)
library(magrittr)

# monivator of the logistic results
basic_path <- file.path("/home/huff/project/immune_checkpoint/result_20171025")
res_path <- file.path(basic_path,"ICP_score.new/logistic_model_predict_Response/all_sample_togather_GSVA")

time_1 <- gsub(" ",":",gsub("\\/",":",format(Sys.time(), "%d/%m/%y %H:%M:%OS3"))) %>% 
  strsplit(":") %>% .[[1]] %>% as.numeric() %>% .[-3]


  res <- read_csv(file.path(res_path,"logistic_record.txt"),quote = " ") %>%
    tidyr::separate(col =`auc.train auc.test auc_validation_Kim auc_validation_Hugo auc_validation_Auslander auc_validation_Auslander auc_validation_Riaz auc_validation_Riaz auc_validation_Van_Allen Features Times SuccssTime`, into = c("auc.train_on", "auc.test_on", "auc_validation_Kim_pre", "auc_validation_Hugo_pre", "auc_validation_Auslander_pre", "auc_validation_Auslander_on", "auc_validation_Riaz_pre", "auc_validation_Riaz_on", "auc_validation_Van_Allen_pre", "Features", "Times", "SuccssTime"),sep=" ") %>%
    dplyr::filter(auc.train_on!="auc.train") %>%
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


