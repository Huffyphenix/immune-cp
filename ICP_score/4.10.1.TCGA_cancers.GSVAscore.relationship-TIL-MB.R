########################## GSVA score of features are: expression site of ICPs, gene family
########################## get the correlation between GSVA score and TIL, mutation burden
library(magrittr)
library(tidyverse)

# data path ---------------------------------------------------------------
# server 1
basic_path <- file.path("/home/huff/project")
immune_res_path <- file.path(basic_path,"immune_checkpoint/result_20171025")
TCGA_path <- file.path(basic_path,"immune_checkpoint/data/TCGA_data")
gene_list_path <- file.path(basic_path,"immune_checkpoint/checkpoint/20171021_checkpoint")
res_path <- file.path(immune_res_path,"TCGA_GSVAScore")


# load GSVA score ---------------------------------------------------------

GSVA.score.onlytumor <- 
  readr::read_rds(file.path(res_path,"TCGA_cancer_specific.onlyTumor_GSVA.score_ICPs_features.rds.gz"))

# 1. GSVA score correlation with TIL -----------------------------------------
# 1.1. load TIL ----
immunity_path_2 <- "/project/huff/huff/immune_checkpoint/data/immunity"
TIMER_immunity <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic)

# 1.2.function --------
fn_correlation_and_DE <- function(GSVA,TIL){
  GSVA %>%
    dplyr::mutate(barcode = substr(barcode,1,15)) %>%
    tidyr::gather(-barcode,key="Features",value="GSVA_score")%>%
    dplyr::inner_join(TIL,by="barcode") %>%
    dplyr::group_by(Features) %>%
    dplyr::mutate(group = ifelse(GSVA_score >= quantile(GSVA_score,0.5),"GSVA_high","GSVA_low")) %>%
    tidyr::nest(-Features) -> .tmp
  
  # correlation
  .tmp %>%
    dplyr::mutate(cor = purrr::map(data,.f=function(.x){
      cor.test(.x$TIL,.x$GSVA_score,method = "spearman") %>%
        broom::tidy()
    })) %>%
    dplyr::select(-data) -> cor.res
  
  # DE
  .tmp %>%
    dplyr::mutate(DE = purrr::map2(data,Features,.f=function(.x,.y){
      print(.y)
      .x %>%
        dplyr::filter(group == "GSVA_high") %>%
        .$GSVA_score %>%
        mean() -> mean_high
      .x %>%
        dplyr::filter(group == "GSVA_low") %>%
        .$GSVA_score %>%
        mean() -> mean_low
      broom::tidy(wilcox.test(TIL  ~ group, data = .x, alternative = "two.sided")) %>%
        dplyr::mutate(mean_high = mean_high, mean_low = mean_low) %>%
        dplyr::mutate(`log2FC(High/Low)` = log2((mean_high+1)/(mean_low+1))) 
    })) %>%
    dplyr::select(-data)  -> DE.res
  
  cor.res %>%
    dplyr::inner_join(DE.res,by="Features")
}

# 1.3.calculation ------
GSVA.score.onlytumor %>% 
  head(1) %>%
  dplyr::mutate(res = purrr::map(GSVA,.f=fn_correlation_and_DE,
                                 TIL = TIMER_immunity %>% 
                                   dplyr::select(barcode,TIL))) %>%
  dplyr::select(-GSVA) %>%
  tidyr::unnest() -> GSVA.TIL.res

# 2. GSVA score with mutation burden ----------------------------
# 2.1. load mutation data ------
mutation_burden_class <- readr::read_rds(file.path("/project/huff/huff/data/TCGA","classfication_of_26_cancers_by_mutation_burden192.rds.gz"))
snv <- readr::read_rds(file.path(tcga_path, "pancan33_snv_from_syn7824274_spread.rds.gz"))
# calculation --------
