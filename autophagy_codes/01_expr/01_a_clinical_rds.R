`%>%` <- magrittr::`%>%`
path <- "/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp"
files.names <- list.files(path="/extraspace/TCGA/TCGA_clinical/",pattern = "_clinical_clean.txt$")
tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"



cancer_types <- files.names %>% stringr::str_split(pattern = "\\_", simplify = T) %>% .[,1]


process_raw_data <- function(.x){
  print(.x)
  path <- "/extraspace/TCGA/TCGA_clinical/"
  d <- readr::read_tsv(file = file.path(path, .x), progress = F)

}

# read all clinical data
cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)
cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("process_raw_data", process_raw_data) -> cancers_names_shards
cancers_names_shards %>%
  dplyr::mutate(clinical = purrr::map(.x = names, process_raw_data)) %>%
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(names, PARTITION_ID)) -> pancan_clinical
on.exit(parallel::stopCluster(cluster))
pancan_clinical %>% readr::write_rds(path = file.path(tcga_path, 'pancan34_clinical.rds.gz'), compress = "gz")


# stage data
stage_fun <- function(.x){
  stage_type <- c("Stage I","Stage II","Stage III","Stage IV")
  if(!tibble::has_name(.x, "pathologic_stage")){
    return(tibble::tibble())
  } else{
    .x %>% 
      dplyr::select(
        barcode,
        gender,
        race, 
        stage = pathologic_stage,
        time = os_days,
        status = os_status
      ) %>% 
      dplyr::filter(stage %in% stage_type) %>% 
      dplyr::mutate(status = plyr::revalue(status, c("Alive" = 0, "Dead" = 1)))
  }
}

pancan_clinical %>% 
  dplyr::mutate(stage = purrr::map(.x = clinical, stage_fun)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(n = nrow(stage)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n > 0) %>% 
  dplyr::select(cancer_types, stage, n) -> pancan_clinical_stage

pancan_clinical_stage %>% 
  readr::write_rds(path = file.path(tcga_path, 'pancan34_clinical_stage.rds.gz'), compress = "gz")

# age data
age_fun <- function(.x){
  if(!tibble::has_name(.x, "age_at_initial_pathologic_diagnosis")){
    return(tibble::tibble())
  } else{
    .x %>% 
      dplyr::select(
        barcode,
        gender,
        race, 
        age = age_at_initial_pathologic_diagnosis ,
        time = os_days,
        status = os_status
      ) %>% 
      dplyr::filter(!is.na(age)) %>% 
      dplyr::mutate(status = plyr::revalue(status, c("Alive" = 0, "Dead" = 1)))
  }
}

pancan_clinical %>% 
  dplyr::mutate(age = purrr::map(.x = clinical, age_fun)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(n = nrow(age)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n > 0) %>% 
  dplyr::select(cancer_types, age, n) -> pancan_clinical_age

pancan_clinical_age %>% 
  readr::write_rds(path = file.path(tcga_path, 'pancan34_clinical_age.rds.gz'), compress = "gz")

#subtype
subtype <- readr::read_tsv(file = file.path(tcga_path, "cancer_subtype.txt"))
subtype_fun <- function(.y, .x, subtype = subtype){
  if(! .y %in% subtype$cancer){
    return(tibble::tibble())
  } else{
    sel <- dplyr::filter(subtype, cancer == .y)
    select_catgory <- sel %>% .$category 
    sb <- .x %>% dplyr::select_("barcode", select_catgory)
    names(sb) <- c("barcode", "subtype")
    
    .x %>% 
      dplyr::select(
        barcode,
        gender,
        race, 
        time = os_days,
        status = os_status
      ) %>% 
      dplyr::mutate(status = plyr::revalue(status, c("Alive" = 0, "Dead" = 1))) %>% 
      dplyr::left_join(sb, by = "barcode") %>% 
      dplyr::filter(!is.na(subtype))
  }
}

pancan_clinical %>% 
  dplyr::mutate(subtype = purrr::map2(cancer_types, clinical, subtype_fun, subtype = subtype)) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(n = nrow(subtype)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(n > 0) %>% 
  dplyr::select(cancer_types, subtype, n) -> pancan_clinical_subtype

pancan_clinical_subtype %>% 
  readr::write_rds(path = file.path(tcga_path, 'pancan34_clinical_subtype.rds.gz'), compress = "gz")








