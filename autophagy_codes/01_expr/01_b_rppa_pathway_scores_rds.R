path <- "/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp"
files.names <- list.files("/extraspace/TCGA/TCGA_protein/", pattern="*re_20160627")

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
library(magrittr)

tcga_path = c("/project/huff/huff/immune_checkpoint/data/TCGA_data")
result_path<-c("/project/huff/huff/immune_checkpoint/result_20171025")
expr_path<-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}
fun_median_cluster <- function(cancer_types, filter_expr){
  cancer <- cancer_types
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type == "01") %>%
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(-type) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::mutate(cluster = ifelse(expr >= median(expr), 2, 1)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(cancer_types = cancer) %>% 
    tidyr::nest(-cancer_types, .key = median_cluster)
}

gene_list_expr %>% 
  purrr::pmap(fun_median_cluster) %>% 
  dplyr::bind_rows() -> gene_list_expr_median_cluster

gene_list_expr_median_cluster %>% 
  readr::write_rds(file.path(expr_path, ".rds_01_b_rppa_median_cluster_expr.rds.gz"), compress = "gz")


>>>>>>>>>>>Use ready made data by liucj.
cancer_types <- files.names %>% stringr::str_split(pattern = "\\_", simplify = T) %>% .[,1]

fun_z_score <- function(rppa){
  med <- median(rppa)
  med_sd <- sd(rppa)
  z_score <- (rppa - med) / med_sd
}
fun_cor <- function(pat, name, .d){
  # print(pat)
  .d %>% 
    dplyr::select(-rppa) %>% 
    dplyr::filter(stringr::str_detect(string = protein, pattern = pat)) %>% 
      tidyr::spread(key = protein, value = z_score) %>% 
      dplyr::select(-barcode) %>% 
      cor() -> .d_cor
  
  switch (name,
          EBP1.65.37 = .d_cor[1, 2],
          EBP1.65.70 = .d_cor[1, 3],
          EBP1.37.70 = .d_cor[2, 3],
          .d_cor[1,2])
}
fun_check <- function(pat, name, type, pathway, .d){
  dup <- c("Akt_", "GSK3.alpha.beta_|GSK3_pS9", "EGFR_", "S6_")
  if(pat %in% dup){
    .d %>% 
      dplyr::filter(stringr::str_detect(string = protein, pattern = pat)) %>% 
      tidyr::spread(key = protein, value = rppa) %>% 
      dplyr::select(-barcode) %>% 
      cor() -> .d_cor
    
    switch (name,
            EBP1.65.37 = .d_cor[1, 2],
            EBP1.65.70 = .d_cor[1, 3],
            EBP1.37.70 = .d_cor[2, 3],
            .d_cor[1,2]) -> .d_cor_val
    
    .d %>% 
      dplyr::filter(stringr::str_detect(string = protein, pattern = pat)) %>% 
      dplyr::group_by(barcode) %>% 
      dplyr::summarise(rppa = ifelse(.d_cor_val > 0.85, mean(rppa), sum(rppa))) %>% 
      dplyr::mutate(name = name, type = type, pathway = pathway)
    
  } else{
    .d %>% 
      dplyr::filter(stringr::str_detect(string = protein, pattern = pat)) %>% 
      dplyr::group_by(barcode) %>% 
      dplyr::summarise(rppa = sum(rppa)) %>% 
      dplyr::mutate(name = name, type = type, pathway = pathway)
  }
}

process_raw_data <- function(.x){
  print(.x)
  path <- "/extraspace/TCGA/TCGA_protein/"
  d <- 
    readr::read_tsv(file = file.path(path, .x), progress = F) %>% 
    dplyr::select(- dplyr::starts_with("X"))

  d %>% 
    tidyr::gather(key = barcode, value = rppa, -protein) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type %in% c("01", "06")) %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::distinct(protein, barcode, .keep_all = T) %>% 
    tidyr::replace_na(replace = list(rppa = 0)) %>% 
    dplyr::group_by(protein) %>% 
    dplyr::mutate(z_score = (rppa - median(rppa)) / sd(rppa)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(protein, barcode, rppa, z_score) -> d_z
  
  score <- 
    tibble::tibble(pat = c("Akt_", "GSK3-alpha-beta_|GSK3_pS9", "p27_", "EGFR_", "S6_", "4E-BP1_p", "4E-BP1_p", "4E-BP1_p", "Src_", "ER-alpha"), 
                          name = c("AKT", "GSK", "p27", "pEGFR", "pS6", "EBP1.65.37", "EBP1.65.70", "EBP1.37.70", "pSrc", "ERalpha"))
  
  score %>% 
    dplyr::mutate(p_cor = purrr::map2_dbl(.x = pat, .y = name, .f = fun_cor, .d = d_z)) %>% 
    dplyr::select(-pat) -> score_cor
#----------------------------------------------------------------------------------------------------------  
  d_z %>% 
    dplyr::select(-z_score) %>% 
    dplyr::mutate(protein = gsub("-R-C|-R-V|-M-C|-M-V|-R-E|-G-C|-R-E|-M-E","",protein)) %>% 
    dplyr::group_by(protein, barcode) %>% 
    dplyr::mutate(rppa = mean(rppa)) %>% 
    dplyr::ungroup() %>% 
    dplyr::distinct() -> d_m
  
  # AKT <- purrr::pmap(.x = "Akt_", .y = "AKT", type = "p", pathway = "PI3KAKT", .f = fun_check, .d = d_m)
  # AKT <- fun_check(pat = "Akt_", name = "AKT", type = "p", pathway = "PI3KAKT")
  # GSK <- fun_check(pat = "GSK3.alpha.beta|GSK3_pS9", name = "GSK", type = "p", pathway = "PI3KAKT")
  # pEGFR <- fun_check(pat = "EGFR_", name = "pEGFR", type = "p", pathway = "PI3KAKT")
  # pS6 <- fun_check(pat = "S6_", name = "pS6", type = "p", pathway = "PI3KAKT")
  
  #pathways -------------------------------------------------------------------------------------
  # pi3kakt
  pi3kakt_pat <- c("Akt_", "GSK3.alpha.beta|GSK3_pS9", "p27_", "PRAS40_", "Tuberin_pT", "INPP4B", "PTEN")
  pi3kakt_name <- c("AKT", "GSK", "p27", "PRAS40", "Tuberin_pT", "INPP4B", "PTEN")
  pi3kakt_type <- c("p", "p", "p", "p", "p", "n", "n")
  pi3kakt <- rep("PI3KAKT", length(pi3kakt_pat))
  
  # RASMAPK
  rasmapk_pat <- c("A.Raf_|c\\.Jun_pS73|C.Raf_|JNK_|MAPK_|MEK1_|p38_p|p90RSK_p|Shc_p|YB-1_p")
  rasmapk_name <- c("cj")
  rasmapk_type <- c("p")
  rasmapk <- c("RASMAPK")
  
  # RTK
  rtk_pat <- c("HER2_|HER3_|Ret_|Shc_p|Src_p", "EGFR_")
  rtk_name <- c("cj", "EGFR_")
  rtk_type <- c("p", "p")
  rtk <- c("RTK", "RTK")
  
  # TSCmTOR
  tscmtor_pat <- c("4E.BP1_p|mTOR_|p70S6K_|Rictor_", "S6_")
  tscmtor_name <- c("cj", "pS6")
  tscmtor_type <- c("p", "p")
  tscmtor <- c("TSCmTOR", "TSCmTOR")
  
  #Apoptosis 
  apoptosis_pat <- c("Bak|Bid|Bim|Caspase.|Bax", "Bcl.2|Bcl.xL|Bad_|cIAP")
  apoptosis_name <- c("cj", "n")
  apoptosis_type <- c("p", "n")
  apoptosis <- c("Apoptosis", "Apoptosis")
  
  # Cell cycle
  cellcycle_pat <- c("CDK1|Cyclin|p27_|PCNA")
  cellcycle_name <- c("cj")
  cellcycle_type <- c("p")
  cellcycle <- c("CellCycle")
  
  #DNA damage
  dnad_pat <- c("53BP1|ATM|BRCA2|Chk1_|Chk2_|Ku80|Mre11|PARP|Rad50|Rad51|XRCC1|p53")
  dnad_name <- c("cj")
  dnad_type <- c("p")
  dnad <- c("DNADamage")
  
  # EMT 
  emt_pat <- c("Collagen|Fibronectin|N.Cadherin", "Claudin.7|E.Cadherin")
  emt_name <- c("cj", "cjl")
  emt_type <- c("p", "n")
  emt <- c("EMT", "EMT")
  
  # ER
  er_pat <- c("ER-|PR")
  er_name <- c("cj")
  er_type <- c("p")
  er <- c("Hormone ER")
  
  # AR
  ar_pat <- c("AR|INPP4B|GATA3|Bcl.2")
  ar_name <- c("cj")
  ar_type <- c("p")
  ar <- c("Hormone AR")
  
  
  path_way <-
    tibble::tibble(
      pat = c(pi3kakt_pat, rasmapk_pat, rtk_pat, tscmtor_pat, apoptosis_pat, cellcycle_pat, dnad_pat, emt_pat, er_pat, ar_pat), 
      name = c(pi3kakt_name, rasmapk_name, rtk_name, tscmtor_name, apoptosis_name, cellcycle_name, dnad_name, emt_name, er_name, ar_name),
      type = c(pi3kakt_type, rasmapk_type, rtk_type, tscmtor_type, apoptosis_type, cellcycle_type, dnad_type, emt_type, er_type, ar_type),
      pathway = c(pi3kakt, rasmapk, rtk, tscmtor, apoptosis, cellcycle, dnad, emt, er, ar),
      .d = list(d_m)
    )
  
  path_way %>% 
    purrr::pmap(fun_check) %>% 
    dplyr::bind_rows() %>% 
    tidyr::spread(key = type, value = rppa) %>% 
    # dplyr::mutate(n = ifelse(tibble::has_name(., "n"), n, 0)) %>% 
    tidyr::replace_na(replace = list(n = 0, p = 0)) %>% 
    dplyr::group_by(barcode, pathway) %>% 
    dplyr::summarise(score = sum(p) - sum(n)) %>% 
    dplyr::ungroup()
  
  }

cancers_names <- tibble::tibble(names = files.names, cancer_types = cancer_types)
cancers_names %>%
  head(1) %>%
  dplyr::mutate(rppa = purrr::map(names, process_raw_data)) -> te

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
cancers_names %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_z_score", fun_z_score)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type)  %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("process_raw_data", process_raw_data)  %>%
  multidplyr::cluster_assign_value("fun_cor", fun_cor) %>% 
  multidplyr::cluster_assign_value("fun_check", fun_check) %>% 
  dplyr::mutate(rppa = purrr::map(names, process_raw_data)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID, -names) -> pancan_rppa
on.exit(parallel::stopCluster(cluster))


pancan_rppa %>% readr::write_rds(path = file.path(tcga_path, 'pancan_rppa_score.rds.gz'), compress = "gz")


