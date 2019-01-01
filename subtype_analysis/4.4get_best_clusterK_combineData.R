#load library
library(methods)
library(magrittr)
library(CancerSubtypes)
library(ConsensusClusterPlus)
library(SNFtool)
library(survival)
library(survminer)
library(ggplot2)

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))

results <- readr::read_rds(file.path("/project/huff/huff/immune_checkpoint/genelist_data",".rds_PanCan28_combine-expr-cnv-methy_snf_20.rds.gz"))

# get best K for cluster results  -----------------------------------------
source("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")

result_path <- file.path("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/combine/Get_best_clutser_20")
W <- results$distanceMatrix
# pdf(file.path(result_path,"Get_best_clutser_20.pdf"))
for (i in 2:20) {
  C <- i
  group = spectralClustering(W,C)
  # Figure 1
  # displayClusters(W, group)
  
  # Figure 2
  group_statistic <- group %>% table()
  less_than_10<- names(group_statistic[group_statistic<10])
  all_clusters <- names(group_statistic)
  more_than_10 <- setdiff(all_clusters,less_than_10)
  data.frame(sample = names(group),group=group,time = time.combine,status = status.combine) %>%
    dplyr::as.tbl() %>%
    dplyr::filter(group %in% more_than_10) %>%
    dplyr::rename("barcode"="sample") %>%
    dplyr::left_join(combined.cancer_info,by="barcode") -> group_survival_data
  color_list = rainbow(length(more_than_10))
  title <- paste("Combined Survival for",C,"Clusters",sep=" ")
  sur_name <-  paste("Combined_Survival_for",C,"Clusters.png",sep="_")
  if (length(more_than_10)>=2) {
    fn_survival(group_survival_data,title,color_list,sur_name)
  }else{
    print("Too small groups")
  }
  
  
  # Figure 3
  group_survival_data %>%
    dplyr::left_join(mutation_burden_class,by="barcode") %>%
    dplyr::mutate(sm_count = ifelse(is.na(sm_count),0,sm_count)) %>%
    dplyr::rename("cancer_types"="cancer_types.x")-> group_cluster_mutation
  comp_mat <- combn(more_than_10,2)
  comp_list <- list()
  m_name <-  paste("Combined_mutation_for",C,"Clusters.png",sep="_")
  for (col in 1:ncol(comp_mat)) {
    comp_list[[col]] <- comp_mat[,col]
  }
  fn_mutation_burden(group_cluster_mutation,color_list,comp_list,m_name)
  
  # Figure 4
  m_a_name <-  paste("Combined_all_mutation_for",C,"Clusters.png",sep="_")
  fn_mutation_burden_all(group_cluster_mutation,color_list,comp_list,m_a_name)
}

