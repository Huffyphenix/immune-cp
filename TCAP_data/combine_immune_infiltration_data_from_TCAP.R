file.name <- list.files(path = "/home/huff/project/data/TCGA/immune_infiltration/miao_TCAP_prediction", pattern = "_immune_infiltration_ref.Rdata$")
cancer_types <- file.name %>% stringr::str_split(pattern = "\\_", simplify = T) %>% .[,1]

process_raw_data <- function(.x){
  print(.x)
  path <- "/home/huff/project/data/TCGA/immune_infiltration/miao_TCAP_prediction"
  load(file.path(path,.x))
  rownames <- gsub("\\.","-",substr(rownames(fra),1,12))
  cancer_type <- strsplit(x = .x,split = "_")[[1]][1]
  fra %>%
    dplyr::as.tbl() %>%
    dplyr::mutate(barcode = rownames) %>%
    tidyr::gather(-barcode,key="TIL",value="value") %>%
    dplyr::mutate(value = as.numeric(value)) %>%
    tidyr::spread(key="TIL",value="value")-> .tmp
  return(.tmp)
}

cancers_names <- tibble::tibble(names = file.name, cancer_types = cancer_types)

cancers_names %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(Infiltration=purrr::map(names,process_raw_data)) -> pancan_immune_infiltration_by_TCAP

pancan_immune_infiltration_by_TCAP %>%
  dplyr::select(-names) %>%
  dplyr::ungroup() %>%
  readr::write_rds(file.path("/home/huff/project/data/TCGA/immune_infiltration/miao_TCAP_prediction","pancan33_immune_infiltration_by_TCAP.rds.gz"),compress = "gz")
