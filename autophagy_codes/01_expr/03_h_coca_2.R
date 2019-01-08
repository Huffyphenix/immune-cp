#load library
library(methods)
library(magrittr)
library(CancerSubtypes)


# Integrative clustering analysis
# Use mRNA, methylation and copy number
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")


load(file = file.path(expr_path_a, ".rda_03_h_coca.rda"))
cnv_cc <- ExecuteCC(clusterNum = 3, d = cnv_matrix, maxK = 10, clusterAlg = "hc", distance = "pearson", title = "cnv_matrix", plot = F)
cnv_cc %>% readr::write_rds(file.path(expr_path_a, ".rds_03_h_coca_cnv_cc.rds.gz"), compress = 'gz')
