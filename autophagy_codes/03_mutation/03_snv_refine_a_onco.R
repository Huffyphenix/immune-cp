library(magrittr)
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")
snv_path <- "/home/cliu18/liucj/projects/6.autophagy/04_snv"
suppressPackageStartupMessages(require(maftools))

gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz"))

atg <- gene_list %>% dplyr::filter(status == "a")
lys <- gene_list %>% dplyr::filter(status == "l")

readr::read_rds(path = file.path(tcga_path, "syn_mutation_syn7824274_mc3_public.maf.rds.gz")) -> mc3_maf



atg_maf <- subsetMaf(maf = mc3_maf, mafObj = T, genes = atg$symbol)




save.image(file = file.path(snv_path, ".rda_03_snv_refine_a_onco.rda"))
load(file = file.path(snv_path, ".rda_03_snv_refine_a_onco.rda"))
