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
index=which(is.na(PanCan26_gene_list_expr_matrix)) # no NA value

index=which(is.na(PanCan26_gene_list_cnv_matrix)) # no NA value


index=which(is.na(PanCan26_gene_list_methy_matrix)) # have NA value
PanCan26_gene_list_methy_matrix=data.imputation(PanCan26_gene_list_methy_matrix,fun="median")


## Concencus clustering
setwd(data_result_path)
# expr_survival <- FSbyCox(gene_list_expr_with_mutation_load.matrix.combine,time.combine,status.combine,cutoff=0.05)
# cnv_survival <- FSbyCox(cnv_merge_snv_data.matrix.combine,time.combine,status.combine,cutoff=0.05)
# methy_survival <- FSbyCox(PanCan26_gene_list_methy_matrix.combine,time.combine,status.combine,cutoff=0.05)

combine_data =list(GeneExp=gene_list_expr_with_mutation_load.matrix.combine,methy=PanCan26_gene_list_methy_matrix.combine,cnv=cnv_merge_snv_data.matrix.combine)
combine_snf <- ExecuteSNF(combine_data,clusterNum=5,K=10)
combine_snf %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_combine_snf_5.rds.gz"), compress = 'gz')

group_snf=combine_snf$group
group_snf %>% table()
# distanceMatrix_cc=combine_cc$distanceMatrix
# p_value=survAnalysis(mainTitle="combine_sur_result",time.combine,status.combine,group_cc,
#                      distanceMatrix_cc,similarity=TRUE)



## SNF
# combine_snf <- ExecuteSNF(combine_data, clusterNum=3, K=20, alpha=0.5, t=20,title = "combine_SNF")
# combine_snf %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_combine_snf.rds.gz"), compress = 'gz')

# ## SNF and CC
# combine_snfcc=ExecuteSNF.CC(combine_data, clusterNum=3, K=20, alpha=0.5, t=20,maxK = 10, pItem = 0.8,reps=500, 
#                         title = "combine_SNF.CC", finalLinkage ="average") # what the meaning of result pic 11?12?13?
# combine_snfcc %>% readr::write_rds(file.path(data_result_path, ".rds_PanCan28_combine_snfcc.rds.gz"), compress = 'gz')
