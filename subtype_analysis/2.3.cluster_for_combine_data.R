###### do subtype analysis for each data type
#load library
library(methods)
library(magrittr)
library(CancerSubtypes)
library(ConsensusClusterPlus)
library(SNFtool)

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))


# all data -----------------------------------------------------------------
# data.checkDistribution(PanCan26_gene_list_cnv_matrix)
index=which(is.na(gene_list_expr_with_mutation_load.matrix.combine)) # no NA value

index=which(is.na(cnv_merge_snv_data.matrix.combine)) # no NA value


index=which(is.na(PanCan26_gene_list_methy_matrix.combine)) # have NA value
PanCan26_gene_list_methy_matrix.combine=data.imputation(PanCan26_gene_list_methy_matrix.combine,fun="median")


## Concencus clustering
setwd(data_result_path)
# expr_survival <- FSbyCox(gene_list_expr_with_mutation_load.matrix.combine,time.combine,status.combine,cutoff=0.05)
# cnv_survival <- FSbyCox(cnv_merge_snv_data.matrix.combine,time.combine,status.combine,cutoff=0.05)
# methy_survival <- FSbyCox(PanCan26_gene_list_methy_matrix.combine,time.combine,status.combine,cutoff=0.05)

combine_data =list(GeneExp=gene_list_expr_with_mutation_load.matrix.combine,methy=PanCan26_gene_list_methy_matrix.combine,cnv=cnv_merge_snv_data.matrix.combine)
results <- ExecuteSNF(combine_data,clusterNum=20,K=10)
results %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_combine-expr-cnv-methy_snf_20.rds.gz"), compress = 'gz')
W <- results$distanceMatrix

# get best K for cluster results  -----------------------------------------
# load("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")
# 
# result_path <- file.path("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/combined/")
# pdf(file.path(result_path,"Get_best_clutser_20.pdf"))
# par(mfrow=c(2,2))
# for (i in 2:20) {
#   C <- i
#   group = spectralClustering(W,C)
#   # Figure 1
#   displayClusters(W, group)
#   
#   # Figure 2
#   
#   group_statistic <- group %>% table()
#   less_than_10<- names(group_statistic[group_statistic<10])
#   all_clusters <- names(group_statistic)
#   more_than_10 <- setdiff(all_clusters,less_than_10)
#   data.frame(sample = colnames(W),group=group,time = expr_time,status = expr_status) %>%
#     dplyr::as.tbl() %>%
#     dplyr::filter(group %in% more_than_10)-> group_survival_data
#   color_list = rainbow(length(more_than_10))
#   if (length(more_than_10)>=2) {
#     fn_survival(group_survival_data,paste("Combined PFS for",C,"Clusters",sep=""),color_list)
#   }else{
#     print("Too small groups")
#   }
#   
#   
#   # Figure 3
#   group_survival_data %>%
#     dplyr::rename("barcode" = "sample") %>%
#     dplyr::left_join(mutation_burden_class,by="barcode") %>%
#     dplyr::mutate(sm_count = ifelse(is.na(sm_count),0,sm_count)) -> group_cluster_mutation
#   comp_mat <- combn(more_than_10,2)
#   comp_list <- list()
#   for (col in 1:ncol(comp_mat)) {
#     comp_list[[col]] <- comp_mat[,col]
#   }
#   fn_mutation_burden(group_cluster_mutation,color_list,comp_list)
#   
#   # Figure 4
#   fn_mutation_burden_all(group_cluster_mutation,color_list,comp_list)
#   
# }
# dev.off()



