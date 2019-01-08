library(magrittr)
url_gff <- "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3"
mir_gff <- readr::read_tsv(url(url_gff), comment = "#", col_names = F) 
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"


mir_gff %>% 
  dplyr::select(gf = X9) %>% 
  dplyr::filter(stringr::str_detect(gf, "MIMAT")) %>% 
  tidyr::separate(col = gf, into = c("ID", "Alias", "Name", "Derives"), sep = ";") %>% 
  dplyr::mutate_all(.funs = dplyr::funs(stringr::str_split(.,"=", simplify = T)[,2])) ->
  mir_acc_name

readr::write_rds(mir_acc_name, file.path(tcga_path, "mirna_acc_name.rds.gz"), compress = "gz")
  
