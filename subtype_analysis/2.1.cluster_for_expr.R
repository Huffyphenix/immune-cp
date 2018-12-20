###### do subtype analysis for each data type
#load library
library(methods)
library(magrittr)
library(CancerSubtypes)

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))


# expr -----------------------------------------------------------------
# data.checkDistribution(PanCan26_gene_list_expr_matrix)

# index=which(is.na(PanCan26_gene_list_expr_matrix)) # no NA value

## Concencus clustering
expr_survival <- FSbyCox(PanCan26_gene_list_expr_matrix,expr_time,expr_status,cutoff=0.05)
expr_cc <- ExecuteCC(clusterNum=10,d=expr_survival,maxK=10,clusterAlg="hc",distance="pearson",title = "expr_sur_cc_10")
group_cc=expr_cc$group
distanceMatrix_cc=expr_cc$distanceMatrix
expr_cc %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_expr_cc_10.rds.gz"), compress = 'gz')
pdf("expr_sur_cc_10.pdf")
p_value=survAnalysis(mainTitle="expr_sur_result",expr_time,expr_status,group_cc,
                     distanceMatrix_cc,similarity=TRUE)
dev.off()


## SNF
# expr_snf <- ExecuteSNF(PanCan26_gene_list_expr_matrix, clusterNum=3, K=20, alpha=0.5, t=20,plot = F)
# expr_snf %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_expr_snf.rds.gz"), compress = 'gz')

## SNF and CC
# expr_snfcc=ExecuteSNF.CC(PanCan26_gene_list_expr_matrix, clusterNum=3, K=20, alpha=0.5, t=20,maxK = 10, pItem = 0.8,reps=500, 
#                      plot = F, finalLinkage ="average") # what the meaning of result pic 11?12?13?
# expr_snfcc %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_expr_snfcc.rds.gz"), compress = 'gz')
