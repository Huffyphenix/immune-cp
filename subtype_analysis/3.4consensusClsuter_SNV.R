#load library
library(methods)
library(magrittr)
library(CancerSubtypes)
library(ConsensusClusterPlus)

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
PanCan26_gene_list_snv_matrix <- readr::read_rds(file.path(data_result_path,"PanCan26_gene_list_snv_matrix.rds.gz"))

# reduce the dataset to the top 5,000 most variable genes, measured by median absolute deviation(mad). 绝对中位差
d <- PanCan26_gene_list_snv_matrix

index=which(is.na(d)) # no NA value
d[index]=0
# mads=apply(d,1,mad)
# d=d[rev(order(sort(mads)))[1:5000],] # 排序，获得rank值，倒置，取前5000

dc = sweep(d,1, apply(d,1,median,na.rm=T)) ## median center genes

# title=tempdir()
# 
# results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,
#                                title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
# results[[2]][['consensusClass']] #  sample distribution in two clusters

# same as above but with pre-computed distance matrix, useful for large datasets (>1,000's of items)
setwd("/project/huff/huff/github/immune-cp/subtype_analysis/ConsensusClusterResult")
dt = as.dist(1-cor(dc,method="pearson"))
results = ConsensusClusterPlus(dt,maxK=10,reps=100,pItem=0.8,pFeature=1,title="SNV_CC",distance="binary",clusterAlg="hc",seed=1262118388.71279)

results %>%
  readr::write_rds(file.path("/project/huff/huff/immune_checkpoint/genelist_data","genelist_SNV_CC_10.rds.gz"),compress = "gz")