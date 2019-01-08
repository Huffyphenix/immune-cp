# 3. mean-center data, apply OCLR to the samples labeled SC (which included both embryonic stem cells, as well as induced pluripotent stem cells)
library(magrittr)
pcbc_dir <- "/home/cliu18/liucj/projects/6.autophagy/synapse/PCBC"
tcga_dir <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
all_filename <- list.files(path = pcbc_dir, pattern = "^SC.*_genes.fpkm_tracking")

samples_expr_filter <- readr::read_rds(path = file.path(pcbc_dir ,"01_pcbc_data_samples_expr_filter.rds.gz"))
# samples_expr_mapped <- readr::read_rds(path = file.path(pcbc_dir, "01_pcbc_data_samples_expr_mapped.rds.gz"))

gene_order <- readr::read_rds(path = file.path(pcbc_dir, "02_OCLR_gene_order.rds.gz"))

setNames(as.list(rep(0, length(gene_order))), gene_order) -> replace_na
  
samples_expr_filter %>% 
  dplyr::select(sample, expr_mapped) %>% 
  tidyr::unnest() %>% 
  dplyr::select(-ens_id, -status) %>% 
  tidyr::spread(key = symbol, value = expr) %>% 
  dplyr::select(-sample) %>% 
  tidyr::replace_na(replace = replace_na) ->
  samples_data

# use gelnet to train one class model
# scale = False -> Mean-center data
log2(as.matrix(samples_data) + 0.1)  -> feature_matrix
feature_matrix[is.na(feature_matrix)] <- 0
readr::write_rds(feature_matrix, path = file.path(pcbc_dir, "02_OCLR_feature_matrix.rds.gz"), compress = "gz")
gene_order <- colnames(feature_matrix)
readr::write_rds(gene_order, path = file.path(pcbc_dir, "02_OCLR_gene_order.rds.gz"), compress = "gz")

# use syn2701943 normalized data
readr::read_rds(path = file.path(pcbc_dir, "syn2701943_matrix_filter.rds.gz")) -> feature_matrix
gene_order <- colnames(feature_matrix)
readr::write_rds(gene_order, path = file.path(pcbc_dir, "02_OCLR_gene_order.rds.gz"), compress = "gz")

#----------------------------------------------
# create interaction matrix
# ppi from BIOGRID
ppi <- readr::read_rds(path = "/home/cliu18/liucj/reference/PPI/source/all_PPI_data.rds") %>% 
  dplyr::filter(P1_species == 9606, P2_species == 9606) %>% 
  dplyr::select(p1 = P1_symbol, p2 = P2_symbol)

fn_interaction <- function(.d, .ppi = ppi, .gene = gene_order){
  .ppi %>% 
    dplyr::filter(p1 == .d | p2 == .d) %>% 
    tidyr::gather(key = name, value = p) %>% 
    dplyr::filter(p != .d) %>% 
    dplyr::pull(p) %>% 
    unique() -> .inter
  ifelse(.gene %in% .inter, 1, 0) %>% 
    tibble::enframe() %>% 
    tidyr::spread(key = name, value = value)
}

tibble::tibble(gene = gene_order) -> gene_order_tibble

# gene_order_tibble %>% 
#   dplyr::mutate(out = purrr::map(.x = gene, .f = fn_interaction, .ppi = ppi, .gene = gene_order)) %>% 
#   dplyr::select(-gene) %>% 
#   tidyr::unnest() -> inter_matrix

cl <- 40
cluster <- multidplyr::create_cluster(ifelse(nrow(gene_order_tibble) > cl, cl, nrow(gene_order_tibble)))
gene_order_tibble %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_interaction", fn_interaction)  %>%
  multidplyr::cluster_assign_value("gene_order", gene_order)  %>%
  multidplyr::cluster_assign_value("ppi", ppi)  %>%
  dplyr::mutate(out = purrr::map(.x = gene, .f = fn_interaction, .ppi = ppi, .gene = gene_order)) %>%
  dplyr::collect() %>%
  dplyr::ungroup() %>%
  dplyr::select(-gene, -PARTITION_ID) -> inter_matrix_tibble
parallel::stopCluster(cluster)

inter_matrix_tibble %>% 
  readr::write_rds(path = file.path(pcbc_dir, "02_OCLR_inter_matrix_tibble.rds.gz"), compress = "gz")

inter_matrix_tibble %>% tidyr::unnest() -> inter_matrix

colnames(inter_matrix) <- gene_order

readr::write_rds(x = inter_matrix, path = file.path(pcbc_dir, "02_OCLR_inter_matrix.rds.gz"), compress = "gz")

inter_matrix %>% as.data.frame() %>% as.matrix() -> A

gelnet::adj2nlapl(A & t(A)) -> gene_interaction_matrix

readr::write_rds(x = gene_interaction_matrix, path = file.path(pcbc_dir, "02_OCLR_gene_interaction_matrix.rds.gz"), compress = "gz")

# gelnet( X = feature_matrix, y = NULL, P = gene_interaction_matrix )
# stem cell score is dot product
# scores= X %*% model2bal$w + model2bal$b
# gene_interaction_matrix <- readr::read_rds(path = file.path(pcbc_dir, "02_OCLR_gene_interaction_matrix.rds.gz"))
model <- gelnet::gelnet(X = feature_matrix, l1 = 0.1, l2 = 1, y = NULL)
.s1 <- feature_matrix %*% model$w
hist(.s1)
.s2 <- apply(feature_matrix, 1, function(.m){cor(.m, model$w, method = "spearman")})
hist(.s2)

res_model <- list(gene_order = gene_order, model = model)
readr::write_rds(res_model, path = file.path(pcbc_dir, "02_OCLR_res_model.rds.gz"), compress = "gz")

save.image(file = file.path(pcbc_dir, "02_OCLR.rda"))
load(file.path(pcbc_dir, "02_OCLR.rda"))
