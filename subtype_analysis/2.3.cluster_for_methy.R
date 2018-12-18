###### do subtype analysis for each data type
#load library
library(methods)
library(magrittr)
library(CancerSubtypes)

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))


# methy -----------------------------------------------------------------
# data.checkDistribution(PanCan26_gene_list_methy_matrix)

index=which(is.na(PanCan26_gene_list_methy_matrix)) # no NA value
PanCan26_gene_list_methy_matrix=data.imputation(PanCan26_gene_list_methy_matrix,fun="median")

## Concencus clustering
methy_survival <- FSbyCox(PanCan26_gene_list_methy_matrix,methy_time,methy_status,cutoff=0.05)
methy_cc <- ExecuteCC(clusterNum=3,d=methy_survival,maxK=10,clusterAlg="hc",distance="pearson",title = "methy_sur_cc")
group_cc=methy_cc$group
distanceMatrix_cc=methy_cc$distanceMatrix
p_value=survAnalysis(mainTitle="methy_sur_result",methy_time,methy_status,group_cc,
                     distanceMatrix_cc,similarity=TRUE)

methy_cc %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_methy_cc.rds.gz"), compress = 'gz')


## SNF
# methy_snf <- ExecuteSNF(PanCan26_gene_list_methy_matrix, clusterNum=3, K=20, alpha=0.5, t=20,plot = F)
# methy_snf %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_methy_snf.rds.gz"), compress = 'gz')

## SNF and CC
# methy_snfcc=ExecuteSNF.CC(PanCan26_gene_list_methy_matrix, clusterNum=3, K=20, alpha=0.5, t=20,maxK = 10, pItem = 0.8,reps=500, 
#                         plot = F, finalLinkage ="average") # what the meaning of result pic 11?12?13?
# methy_snfcc %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_methy_snfcc.rds.gz"), compress = 'gz')
