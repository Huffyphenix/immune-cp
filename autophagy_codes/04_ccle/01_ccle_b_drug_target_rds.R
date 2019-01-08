library(magrittr)
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
drug_path <- "/extraspace/yye1/share_data/DrugData"

t_gdsc <- readr::read_rds(file.path(drug_path, "GDSC_compounds_target_GenePathway.rds")) %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(
    drug_id = DRUG.ID,
    drug_name = DRUG.NAME,
    synonyms = SYNONYMS,
    target_gene = TARGET,
    target_pathway = TARGET.PATHWAY
  ) %>% 
  tidyr::separate_rows(target_gene) %>% 
  tidyr::nest(target_gene, target_pathway, .key = target)

readr::write_rds(t_gdsc, path = file.path(tcga_path, "drug_target_gdsc.rds.gz"), compress = "gz")

# t_ctrp <- read.csv(file.path(drug_path, "CTRP_Drug.target.genes.pathway.csv"))
t_ctrp <- read.csv(file.path(drug_path, "CTRP_Drug.target.genes.csv"), stringsAsFactors = F) %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate(synonyms = NA_character_, target_pathway = NA_character_) %>% 
  dplyr::rename(
    drug_id = master_cpd_id,
    drug_name = cpd_name,
    target_gene = gene_symbol_of_protein_target
  ) %>% 
  dplyr::filter(target_gene != "") %>% 
  tidyr::separate_rows(target_gene) %>% 
  tidyr::nest(target_gene, target_pathway, .key = target)
readr::write_rds(t_ctrp, path = file.path(tcga_path, "drug_target_ctrp.rds.gz"), compress = "gz")
