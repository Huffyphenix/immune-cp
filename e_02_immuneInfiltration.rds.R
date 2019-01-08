##############################
# get the data combine gene expression and immune infiltration
##############################

# load data
pancan_immune_infiltration_by_TCAP <- readr::read_rds(file.path("/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction","pancan33_immune_infiltration_by_TCAP.rds.gz")) %>%
  dplyr::select(-names) %>%
  dplyr::ungroup() %>%
  tidyr::unnest()

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))

# data combination
gene_list_expr %>%
  dplyr::select(-cancer_types) %>%
  dplyr::inner_join(pancan_immune_infiltration_by_TCAP,by="barcode") %>%
  dplyr::select(-PFS,-PFS.time) %>%
  tidyr::nest(-cancer_types,-symbol,-entrez_id) -> genelist_exp_inmmune_infiltration

genelist_exp_inmmune_infiltration %>%
  readr::write_rds(file.path("/project/huff/huff/immune_checkpoint/genelist_data","genelist_exp_immuneInfiltration.rds.gz"),compress = "gz")

