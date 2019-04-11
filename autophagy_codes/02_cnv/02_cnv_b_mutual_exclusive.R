library(methods)
library(magrittr)
library(utils)

out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
cnv_path <- file.path(out_path, "m_1_cnv")
tcga_path <- "/home/huff/project/immune_checkpoint/data/TCGA_data"
expr_path <-c("/home/huff/project/immune_checkpoint/result_20171025/expr_rds")
gene_list_cnv <- readr::read_rds(file.path(cnv_path, '.rds_02_cnv_a_gene_list.rds.gz'))

fn_ex <- function(V1, V2, .data, cancer_types){
  # .data <- filter_cnv
  # V1 <- 'CD276'
  # V2 <- 'VTCN1'
  .data %>% 
    dplyr::filter(symbol %in% c(V1, V2)) %>% 
    tidyr::gather(key = barcode, value = gistic, -symbol) -> plot_ready
  plot_ready %>% 
    tidyr::spread(key = symbol, value = gistic) %>% 
    dplyr::select(-barcode) -> .d
  .g_name <- colnames(.d)
  # colnames(.d) <- c("A", "B")
  name <- paste(c(cancer_types, .g_name), collapse = "_")
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(. == 0)) %>% 
    nrow() -> nn
  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(. != 0)) %>% 
    nrow()-> aa
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. != 0)) %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == 0)) -> .d_an
  
  sum(.d_an %>% dplyr::pull(1) != 0) -> an
  sum(.d_an %>% dplyr::pull(2) != 0) -> na
  c(nn = nn, an = an, na = na, aa = aa) %>% 
    cometExactTest::comet_exact_test(mutmatplot = F) -> p_val
  
  tibble::tibble(te = name, nn = nn, an = an ,na = na, aa = aa, p_val = p_val)
}
fn_mutal_exclusive <- function(cancer_types, filter_cnv, cluster){
  # cancer_types <- te$cancer_types
  # filter_cnv <- te$filter_cnv[[1]]
  filter_cnv %>% 
    dplyr::pull(symbol) %>% 
    combn(m = 2) %>% 
    t() %>% 
    dplyr::as_data_frame() -> .gene_pairs
  
  .gene_pairs %>% 
    # multidplyr::partition(cluster = cluster) %>%
    # multidplyr::cluster_library("magrittr") %>%
    # multidplyr::cluster_assign_value("fn_ex", fn_ex) %>%
    # multidplyr::cluster_assign_value("filter_cnv", filter_cnv) %>%
    # multidplyr::cluster_assign_value("cancer_types", cancer_types) %>%
    dplyr::mutate(rs = purrr::map2(V1, V2, .f = fn_ex, .data = filter_cnv, cancer_types = cancer_types)) %>% 
    dplyr::collect() %>%
    dplyr::as_tibble() %>%
    dplyr::ungroup() %>%
    #dplyr::select(-PARTITION_ID) %>% 
    dplyr::select(rs) %>% 
    tidyr::unnest() %>% 
    tidyr::separate(col = te, into = c('cancer_types', 'g1', 'g2')) -> .gene_pairs_pval
  
  .gene_pairs_pval %>% 
    dplyr::mutate(fdr = p.adjust(p_val, method = 'fdr'))
}

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_cnv %>% 
  # dplyr::filter(cancer_types == "OV") %>% 
  purrr::pmap(.f = fn_mutal_exclusive, cluster = cluster) %>% 
  dplyr::bind_rows() -> mutual_exclusive
parallel::stopCluster(cluster)

mutual_exclusive %>% 
  readr::write_rds(path = file.path(cnv_path, ".rds_02_cnv_b_mutual_exclusive.rds.gz"), compress = 'gz')

save.image(file = file.path(cnv_path, ".rds_02_cnv_b_mutual_exclusive.rds.gz"))
load(file = file.path(cnv_path, ".rd_02_cnv_b_mutual_exclusive.rda"))
