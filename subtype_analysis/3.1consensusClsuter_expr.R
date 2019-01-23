#load library
library(methods)
library(magrittr)
library(CancerSubtypes)
library(ConsensusClusterPlus)
library(SNFtool)

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
# load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))
PanCan26_gene_list_expr_matrix <- readr::read_rds(file.path(data_result_path,"PanCan26_gene_list_expr_matrix.scaled.rds.gz"))

# reduce the dataset to the top 5,000 most variable genes, measured by median absolute deviation(mad). 绝对中位差
d <- PanCan26_gene_list_expr_matrix
# mads=apply(d,1,mad)
# d=d[rev(order(sort(mads)))[1:5000],] # 排序，获得rank值，倒置，取前5000

dc = sweep(d,1, apply(d,1,median,na.rm=T)) ## median center genes

# title=tempdir()
# 
# results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
#                                title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
# results[[2]][['consensusClass']] #  sample distribution in two clusters

# same as above but with pre-computed distance matrix, useful for large datasets (>1,000's of items)
# setwd("/project/huff/huff/github/immune-cp/subtype_analysis/ConsensusClusterResult")
dt = as.dist(1-cor(dc,method="pearson"))
results = ConsensusClusterPlus(dt,maxK=20,reps=100,pItem=0.8,pFeature=1,title="expr_CC",distance="pearson",clusterAlg="hc",seed=1262118388.71279)

results %>%
  readr::write_rds(file.path("/project/huff/huff/immune_checkpoint/genelist_data","genelist_expr_scaled_CC_20.rds.gz"),compress = "gz")


# get best K for cluster results  -----------------------------------------
# load("/project/huff/huff/github/immune-cp/subtype_analysis/funtions_to_draw_pic.R")
# 
# result_path <- file.path("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/expr/")
# for (i in 2:20) {
#   C <- i
#   W <- results[[C]][['consensusMatrix']]
#   group <- results[[i]][['consensusClass']]
#   # Figure 1
#   displayClusters(W, group)
#   
#   # Figure 2
#   group_statistic <- group %>% table()
#   less_than_10<- names(group_statistic[group_statistic<10])
#   all_clusters <- names(group_statistic)
#   more_than_10 <- setdiff(all_clusters,less_than_10)
#   data.frame(sample = names(group),group=group,time = expr_time,status = expr_status) %>%
#     dplyr::as.tbl() %>%
#     dplyr::filter(group %in% more_than_10)-> group_survival_data
#   color_list = rainbow(length(more_than_10))
#   if (length(more_than_10)>=2) {
#     fn_survival(group_survival_data,paste("Expr Survival for",C,"Clusters",sep=""),color_list)
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
