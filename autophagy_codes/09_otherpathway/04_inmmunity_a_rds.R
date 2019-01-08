library(magrittr)
maker_gene_path <- "/extraspace/yye1/share_data/markerGenes"
tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
expr_path_a <- file.path(expr_path, "03_a_gene_expr")


expr <- readr::read_rds(file.path(tcga_path, "pancan_expr_20160513.rds.gz"))
gene_list <- readr::read_csv(file.path(maker_gene_path, "immnuo.marker.csv")) 

readr::write_rds(gene_list, file.path(tcga_path, "immunity_raw.rds.gz"), compress = "gz")
