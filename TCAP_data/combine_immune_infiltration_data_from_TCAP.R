file.name <- list.files(path = "/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction2", pattern = "_predicted_result_tumor.Rdata$")
cancer_types <- file.name %>% stringr::str_split(pattern = "\\_", simplify = T) %>% .[,1]

process_raw_data <- function(.x){
  print(.x)
  path <- "/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction2"
  load(file.path(path,.x))
  rownames <- gsub("\\.","-",substr(rownames(predict_result),1,12))
  cancer_type <- strsplit(x = .x,split = "_")[[1]][1]
  predict_result %>%
    as.data.frame() %>%
    dplyr::as.tbl() %>%
    dplyr::mutate(barcode = rownames) %>%
    dplyr::select(barcode, InfiltrationScore) -> .tmp
  return(.tmp)
}

cancers_names <- tibble::tibble(names = file.name, cancer_types = cancer_types)

cancers_names %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(Infiltration=purrr::map(names,process_raw_data)) -> pancan_immune_infiltration_by_TCAP

pancan_immune_infiltration_by_TCAP %>%
  readr::write_rds(file.path("/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction","pancan33_immune_infiltration_by_TCAP.rds.gz"),compress = "gz")
