###### do subtype analysis for each data type
#load library
library(methods)
library(magrittr)
library(CancerSubtypes)

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))


# cnv -----------------------------------------------------------------
# data.checkDistribution(PanCan26_gene_list_cnv_matrix)

index=which(is.na(PanCan26_gene_list_cnv_matrix)) # no NA value
PanCan26_gene_list_cnv_spread %>%
  tidyr::gather(-symbol,key = "sample", value = "cnv") %>%
  .$cnv %>% table()

## Concencus clustering
setwd(data_result_path)
cnv_survival <- FSbyCox(PanCan26_gene_list_cnv_matrix,cnv_time,cnv_status,cutoff=0.05)
cnv_cc <- ExecuteCC(clusterNum=10,d=cnv_survival,maxK=10,clusterAlg="hc",distance="pearson",title = "cnv_sur_cc_10")
cnv_cc %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_cnv_cc_10.rds.gz"), compress = 'gz')

group_cc=cnv_cc$group
distanceMatrix_cc=cnv_cc$distanceMatrix
pdf("cnv_sur_cc_10.pdf")
p_value=survAnalysis(mainTitle="cnv_sur_result",cnv_time,cnv_status,group_cc,
                     distanceMatrix_cc,similarity=TRUE)
dev.off()



## SNF
# cnv_snf <- ExecuteSNF(PanCan26_gene_list_cnv_matrix, clusterNum=3, K=20, alpha=0.5, t=20,title = "CNV_SNF")
# cnv_snf %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_cnv_snf.rds.gz"), compress = 'gz')

## SNF and CC
# cnv_snfcc=ExecuteSNF.CC(PanCan26_gene_list_cnv_matrix, clusterNum=3, K=20, alpha=0.5, t=20,maxK = 10, pItem = 0.8,reps=500, 
#                         title = "CNV_SNF.CC", finalLinkage ="average") # what the meaning of result pic 11?12?13?
# cnv_snfcc %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_cnv_snfcc.rds.gz"), compress = 'gz')
