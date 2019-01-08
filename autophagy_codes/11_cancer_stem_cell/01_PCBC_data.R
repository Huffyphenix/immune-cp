# PCBC data was downloaded from synapse syn2246521
# data from stem cell report paper.

#1. use PCBC data to calculate a stem cell index (sci) based on mRNA expressision
#2. using one class logic regression (OCLR) (Sokolov et al., 2016)

# Process
# 1. map Ensembl ID to HUGO, drop gene that had no mapping.
# 2. resulting trainning with 12945 gene expression measure across PCBC samples.
# 3. mean-center data, apply OCLR to the samples labeled SC (which included both embryonic stem cells, as well as induced pluripotent stem cells)
library(magrittr)
pcbc_dir <- "/home/cliu18/liucj/projects/6.autophagy/synapse/PCBC"
tcga_dir <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
all_filename <- list.files(path = pcbc_dir, pattern = "^SC.*_genes.fpkm_tracking")

tcga_gene_symbol <- readr::read_rds(path = file.path(tcga_dir, "tcga_gene_symbol.rds.gz"))

tibble::tibble(filename = all_filename) %>% 
  dplyr::mutate(sc = stringr::str_split(string = filename, pattern = "-", simplify = T)[,1]) %>% 
  dplyr::mutate(sample = stringr::str_replace(string = filename, pattern = "SC\\d*-", replacement = "")  %>% 
                  stringr::str_replace(pattern = "_genes.fpkm_tracking", replacement = "")) -> samples_info

readr::read_tsv(file.path(pcbc_dir, "SC11-002A.133.1.2_genes.fpkm_tracking")) %>% 
  dplyr::pull(gene_id) %>% 
  stringr::str_replace(pattern = "\\.\\d*", replacement = "") -> pcbc_ens_id

# map ensembl to hgnc
GENES = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
GENES.ATTRIBUTES <- biomaRt::listAttributes(GENES)
GENES.FILTERS <- biomaRt::listFilters(GENES)
biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"), filters = c("ensembl_gene_id"), values = list(ensembl_gene_id = pcbc_ens_id), mart = GENES) -> ens2hgnc

ens2hgnc %>% 
  tibble::as_tibble() %>% 
  dplyr::rename(gene_id = ensembl_gene_id) %>% 
  dplyr::distinct(gene_id, .keep_all = T) -> ens2hgnc_uniqe

fn_load_gene_expr <- function(.x, .dir = pcbc_dir, .map = ens2hgnc_uniqe){
 # .x <- "SC11-002A.133.1.2_genes.fpkm_tracking"
  readr::read_tsv(file = file.path(.dir, .x)) -> .d
  
  .d %>% 
    dplyr::select(gene_id, FPKM, FPKM_status) %>% 
    dplyr::mutate(gene_id = stringr::str_replace(string = gene_id, pattern = "\\.\\d*", replacement = "")) %>% 
    dplyr::inner_join(.map, by = "gene_id") 
  
}

samples_info %>% 
  dplyr::mutate(expr = purrr::map(.x = filename, .f = fn_load_gene_expr, .dir = pcbc_dir, .map = ens2hgnc_uniqe)) -> 
  samples_expr

readr::write_rds(x = samples_expr, path = file.path(pcbc_dir ,"01_pcbc_data_samples_expr.rds.gz"), compress = "gz")
samples_expr <- readr::read_rds(path = file.path(pcbc_dir ,"01_pcbc_data_samples_expr.rds.gz"))



syn2701943 <- readr::read_rds(path = file.path(pcbc_dir, "rnaseq_syn2701943_normlized.rds.gz"))
syn2701943 %>% dplyr::pull(tracking_id) %>% stringr::str_replace(pattern = "\\.\\d*", replacement = "") -> .tmp_name
syn2701943 %>% dplyr::select(-tracking_id) %>% as.matrix() -> syn2701943_mat
rownames(syn2701943_mat) <- .tmp_name

ens2hgnc %>% 
  dplyr::filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>% 
  dplyr::filter(hgnc_symbol %in% tcga_gene_symbol$symbol) %>% 
  dplyr::distinct(ensembl_gene_id, .keep_all = T) %>% 
  dplyr::select(ens_id = ensembl_gene_id, symbol = hgnc_symbol) %>% 
  dplyr::filter(ens_id %in% .tmp_name) %>% 
  tibble::deframe() -> .syn2701943_colnames

syn2701943_mat[,-1] -> syn2701943_mat

.syn2701943_rownames <- 
  colnames(syn2701943_mat)[colnames(syn2701943_mat) %>% stringr::str_detect(pattern = "SC\\d\\d-")]

syn2701943_mat %>% t() -> syn2701943_mat_t

syn2701943_mat_t[.syn2701943_rownames, names(.syn2701943_colnames)] -> syn2701943_matrix
colnames(syn2701943_matrix) <- .syn2701943_colnames
readr::write_rds(x = syn2701943_matrix, path = file.path(pcbc_dir, "syn2701943_matrix_filter.rds.gz"))

fn_filter <- function(.expr, .tcga_gene = tcga_gene_symbol){
  # .expr <- .te$expr[[1]]
  
  .expr %>% 
    dplyr::filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>% 
    dplyr::filter(hgnc_symbol %in% .tcga_gene$symbol) %>% 
    dplyr::select(symbol = hgnc_symbol, ens_id = gene_id, expr = FPKM, status = FPKM_status) %>% 
    dplyr::distinct(symbol, .keep_all = T) %>% 
    dplyr::filter(expr > 0)
    
}

samples_expr %>% 
  dplyr::mutate(expr_mapped = purrr::map(.x = expr, .f = fn_filter)) ->
  samples_expr_mapped



samples_expr_mapped$expr_mapped %>% 
  purrr::map(.f = function(.x){dplyr::filter(.x, status == "OK") %>% dplyr::pull(symbol)}) %>% 
  purrr::reduce(.f = intersect) -> symbol_exclude_status
readr::write_rds(symbol_exclude_status, path = file.path(pcbc_dir, "02_OCLR_gene_order.rds.gz"), compress = "gz")

samples_expr_mapped %>% 
  dplyr::mutate(expr_mapped = purrr::map(.x = expr_mapped, .f = function(.d, .s){dplyr::filter(.d, symbol %in% .s )}, .s = symbol_exclude_status)) %>% 
  dplyr::mutate(sample = plyr::revalue(filename, dplyr::select(samples_info, filename, sample) %>% tibble::deframe())) -> samples_expr_filter

readr::write_rds(x = samples_expr_filter, path = file.path(pcbc_dir ,"01_pcbc_data_samples_expr_filter.rds.gz"), compress = "gz")



