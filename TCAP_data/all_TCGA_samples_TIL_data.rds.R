
load("/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction_for_all_samples/All_TCGA_sample_TIL.Rdata")

all_cancer_result_matrix %>%
  as.data.frame() %>%
  dplyr::mutate(barcode = rownames(all_cancer_result_matrix)) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::mutate(barcode = gsub("\\.","-",barcode)) -> all_cancer_TIL

all_cancer_TIL %>%
  readr::write_rds("/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction_for_all_samples/All_TCGA_sample_TIL.rds.gz",compress = "gz")
