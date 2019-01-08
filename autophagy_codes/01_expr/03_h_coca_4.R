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
coca_snf <- ExecuteSNF(coca, clusterNum=3, K=20, alpha=0.5, t=20, plot = F)
coca_snf %>% readr::write_rds(path = file.path(expr_path_a, ".rds_03_h_coca_coca_snf.rds.gz"), compress = "gz")
