########## get cox score of clinical samples
########## sum of GSVA scre * coef of COX[from TCGA data analysis:4.10.1.TCGA_cancers.GSVAscore.relationship-TIL-MB.R]
library(magrittr)
library(tidyverse)

# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
score_path <- file.path(immune_res_path,"ICP_score/5.GSVA-ICPs_exp_site-all_possible")
res_path <- file.path(score_path,"cox_score")


# load data ---------------------------------------------------------------
load(file.path(res_path,"cox_score.difference_between_response_non-response.rda"))

# coef of Coxph model from tcga
uni_cox.PFS <- readr::read_tsv(file.path(immune_res_path,"TCGA_GSVAScore/all_togather/3.survival_with_GSVA_score.new","GSVA.score.univarite.surv.PFS.tsv"))
multi_cox_PFS <- readr::read_tsv(file.path(immune_res_path,"TCGA_GSVAScore/all_togather/3.survival_with_GSVA_score.new","GSVA.score.multi-varite.surv.PFS.tsv"))

uni_cox.OS <- readr::read_tsv(file.path(immune_res_path,"TCGA_GSVAScore/all_togather/3.survival_with_GSVA_score.new","2.GSVA.score.univarite.surv.OS.tsv"))
multi_cox_OS <- readr::read_tsv(file.path(immune_res_path,"TCGA_GSVAScore/all_togather/3.survival_with_GSVA_score.new","2.GSVA.score.multi-varite.surv.OS.tsv"))

# clinical GSVA score
gsva_clinical <- readr::read_rds(file.path(score_path,"ICP_GSVA_score_all-possible-features_all-togather.rds.gz")) %>%
  dplyr::mutate(cancer_types = "all_cancer")
# gsva_clinical <- readr::read_rds(file.path(score_path,"ICP_GSVA_score_all-possible-features_all-specific.rds.gz")) %>%
#   dplyr::mutate(cancer_types = ifelse(tissue == "skin", "SKCM" ,"STAD"))
# clinical response data
sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::filter(! Response %in% c("NE", "X")) %>%
  dplyr::mutate(Response = ifelse(Response %in% c("CR", "PR", "PRCR", "R"), "Response", "Non-Response"))%>%
  dplyr::select(Run, Biopsy_Time, Response,`data ID`)  %>%
  dplyr::inner_join(readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","anti_ICB_study_info.txt"))%>% dplyr::filter(`data type`=="RNA-seq"),by=c("data ID")) %>%
  dplyr::mutate(group = paste(blockade,Biopsy_Time,Author,sep=","))

# function--------
fn_cox_score <- function(GSVA,survival){
  GSVA %>%
    tidyr::gather(-Run,key="Features",value="GSVAscore") %>%
    dplyr::mutate(Features = gsub(" ","_",Features)) %>%
    dplyr::filter(Features %in% survival$Features) %>%
    dplyr::inner_join(survival,by="Features") %>%
    dplyr::mutate(GSVA_muliple_coef = GSVAscore*coef) %>%
    dplyr::group_by(Run) %>%
    dplyr::mutate(cox_score = sum(GSVA_muliple_coef)) %>%
    dplyr::ungroup() %>%
    dplyr::select(Run,cox_score) %>%
    unique()
}

# 1. multi-cox model coef ------
# 1.1 PFS --------
# 1.1.1 get score ----
multi_cox_PFS %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(gsva_clinical,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> clinical.PFS.multi.coxscore

# clinical.PFS.multi.coxscore %>%
#   readr::write_tsv(file.path(res_path,"1.1.clinical.PFS.multi.coxscore(GSVA*coef_of_multi_cox).tsv"))

# 1.1.2 score difference between response and non-response samples ----
clinical.PFS.multi.coxscore %>%
  dplyr::inner_join(sample_info,by="Run") %>%
  ggplot(aes(x=Response,y=cox_score)) +
  # geom_violin() +
  geom_boxplot(width = 0.1) +
  facet_wrap(".~group",scales = "free") +
  ggpubr::stat_compare_means()
ggsave(file.path(res_path,"clinical.PFS.multi.coxscore_diff-Res-nonRes.png"),device = "png",height = 4,width = 5)

# 1.2 OS --------
# 1.2.1 get score ----
multi_cox_OS %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(gsva_clinical,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> clinical.OS.multi.coxscore

# clinical.PFS.multi.coxscore %>%
#   readr::write_tsv(file.path(res_path,"1.1.clinical.PFS.multi.coxscore(GSVA*coef_of_multi_cox).tsv"))

# 1.2.2 score difference between response and non-response samples ----
clinical.OS.multi.coxscore %>%
  dplyr::inner_join(sample_info,by="Run") %>%
  ggplot(aes(x=Response,y=cox_score)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  facet_wrap(".~group",scales = "free") +
  ggpubr::stat_compare_means()
ggsave(file.path(res_path,"clinical.OS.multi.coxscore_diff-Res-nonRes.png"),device = "png",height = 4,width = 5)

# 2. uni-cox model coef ------
# 2.1 PFS --------
# 2.1.1 get score ----
uni_cox.PFS %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(gsva_clinical,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> clinical.PFS.uni.coxscore

# clinical.PFS.uni.coxscore %>%
#   readr::write_tsv(file.path(res_path,"1.1.clinical.PFS.multi.coxscore(GSVA*coef_of_multi_cox).tsv"))

# 2.1.2 score difference between response and non-response samples ----
clinical.PFS.uni.coxscore %>%
  dplyr::inner_join(sample_info,by="Run") %>%
  ggplot(aes(x=Response,y=cox_score)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  facet_wrap(".~group",scales = "free") +
  ggpubr::stat_compare_means()
ggsave(file.path(res_path,"clinical.PFS.uni.coxscore_diff-Res-nonRes.png"),device = "png",height = 4,width = 5)

# 2.2 OS --------
# 2.2.1 get score ----
uni_cox.OS %>%
  tidyr::nest(-cancer_types,.key="PFS") %>%
  dplyr::inner_join(gsva_clinical,by="cancer_types") %>%
  dplyr::mutate(cox_score = purrr::map2(GSVA,PFS,fn_cox_score)) %>%
  dplyr::select(cancer_types,cox_score) %>%
  tidyr::unnest() -> clinical.OS.uni.coxscore

# clinical.PFS.uni.coxscore %>%
#   readr::write_tsv(file.path(res_path,"1.1.clinical.PFS.multi.coxscore(GSVA*coef_of_multi_cox).tsv"))

# 2.1.2 score difference between response and non-response samples ----
clinical.OS.uni.coxscore %>%
  dplyr::inner_join(sample_info,by="Run") %>%
  ggplot(aes(x=Response,y=cox_score)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  facet_wrap(".~group",scales = "free") +
  ggpubr::stat_compare_means()
ggsave(file.path(res_path,"clinical.OS.uni.coxscore_diff-Res-nonRes.png"),device = "png",height = 4,width = 5)


# 3. GSVAscore with response and non-response -----------------------------
gsva_clinical %>%
  tidyr::unnest() %>%
  dplyr::select(-tissue,-cancer_types) %>%
  tidyr::gather(-Run,key="features",value="gsva") %>%
  tidyr::nest(-features) %>%
  dplyr::mutate(DE_res = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::inner_join(sample_info%>%
                          dplyr::select(Run,group,Response),by="Run") %>%
      dplyr::group_by(group,Response) %>%
      dplyr::mutate(n=n()) %>%
      dplyr::filter(n>5) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(n=length(table(Response))) %>%
      dplyr::filter(n==2) %>%
      dplyr::select(-n) %>%
      tidyr::nest(-group) %>%
      dplyr::mutate(DE = purrr::map(data,.f=function(.x){
        broom::tidy(wilcox.test(gsva ~ Response, data = .x, alternative = "two.sided"))
      })) %>%
      dplyr::select(-data) %>%
      tidyr::unnest()
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> GSVA_score_diff_in_R_NR
  
save.image(file.path(res_path,"cox_score.difference_between_response_non-response.rda"))
