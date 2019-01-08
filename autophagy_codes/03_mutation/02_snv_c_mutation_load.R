library(magrittr)
library(ggplot2)
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
# two main issue needs to be addressed.
# 1. mutation load for every cancer types.
# 2. gene list hypermutated | hypomutated | non significant

# mutation load
# load mutation matrix

mut_mat <- readr::read_rds(file.path(tcga_path, 'syn_mutation_syn7824274_mc3_public.pass.mat.rds.gz'))
sample_info <- readr::read_rds(file.path(tcga_path, "sample_info.rds.gz"))
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
} 

mut_mat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column() %>% 
  tibble::as_tibble() %>% 
  tidyr::gather(key = barcode, value = count, -rowname) %>% 
  dplyr::mutate(sample = stringr::str_replace_all(barcode, pattern = "\\.", replacement = "-") %>% 
                  stringr::str_sub(start = 1, end = 12)) %>% 
  dplyr::select(rowname, count, sample) %>% 
  tidyr::spread(key = rowname, value = count) -> 
  mut_mat_tibble
mut_mat_tibble %>% 
  readr::write_rds(path = file.path(snv_path, ".rds_02_snv_c_mut_mat_tibble.rds.gz"))

mut_mat_tibble %>% 
  dplyr::filter(rowSums(.[-1]) < 1000) -> 
  mut_mat_tibble_1000

sample_info %>% 
  dplyr::select(sample, cancer_types) %>%
  dplyr::arrange(cancer_types) %>% 
  dplyr::distinct() %>% 
  tidyr::nest(sample)  -> sample_info_lite

fn_merge_mut <- function(data, mut_mat_tibble_1000 = mut_mat_tibble_1000){
  # data <- te$data[[1]]
  data %>% 
    dplyr::inner_join(mut_mat_tibble_1000, by = "sample")
}
fn_mutation_load <- function(data,  gene = atg_gene){
  # data <- te$data[[1]]
  n_sample <- nrow(data)
   # gene <- atg_gene
  data %>% 
    tidyr::gather(key = symbol, value = count, -sample) %>%
    dplyr::filter(symbol %in% gene) %>% 
    dplyr::summarise(load = sum(count) / n_sample / length(unique(gene))) %>% 
    dplyr::pull(load)
}


sample_info_lite %>%
  head(2) %>% 
  dplyr::mutate(data = purrr::map(.x = data, .f = fn_merge_mut, mut_mat_tibble_1000)) %>% 
  dplyr::mutate(mutation_load = purrr::map_dbl(.x = data, .f = fn_mutation_load, gene = atg_gene))

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
sample_info_lite %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fn_merge_mut", fn_merge_mut)  %>%
  multidplyr::cluster_assign_value("fn_mutation_load", fn_mutation_load) %>% 
  multidplyr::cluster_assign_value("mut_mat_tibble_1000", mut_mat_tibble_1000) %>% 
  multidplyr::cluster_assign_value("atg_gene", atg_gene) %>% 
  dplyr::mutate(data = purrr::map(.x = data, .f = fn_merge_mut, mut_mat_tibble_1000)) %>% 
  dplyr::mutate(mutation_load = purrr::map_dbl(.x = data, .f = fn_mutation_load, gene = atg_gene)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> sample_info_mutation_load_atg

sample_info_lite %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fn_merge_mut", fn_merge_mut)  %>%
  multidplyr::cluster_assign_value("fn_mutation_load", fn_mutation_load) %>% 
  multidplyr::cluster_assign_value("mut_mat_tibble_1000", mut_mat_tibble_1000) %>% 
  multidplyr::cluster_assign_value("lys_gene", lys_gene) %>% 
  dplyr::mutate(data = purrr::map(.x = data, .f = fn_merge_mut, mut_mat_tibble_1000)) %>% 
  dplyr::mutate(mutation_load = purrr::map_dbl(.x = data, .f = fn_mutation_load, gene = lys_gene)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) -> sample_info_mutation_load_lys

parallel::stopCluster(cluster)

# sample_info_mutation_load_atg -> atg_genes
# sample_info_mutation_load_lys -> lys_genes

sample_info_mutation_load_atg %>% 
  ggplot(aes(x = reorder(cancer_types, - mutation_load), y = mutation_load)) +
  geom_bar(stat = 'identity')

sample_info_mutation_load_lys %>% 
  ggplot(aes(x = reorder(cancer_types, - mutation_load), y = mutation_load)) +
  geom_bar(stat = 'identity')
  

# hyper mutated
# obtain gene length
all_protein_coding_genes <- row.names(mut_mat) 
gencodev19 <- readr::read_tsv('/home/cliu18/liucj/reference/GENCODE/GENODEv25/gencode.v19.gene_base.annotation.gtf', col_names = F) 

fn_extract_symbol <- function(.s){
  # .s <- te$X9[1]
  .s %>% 
    stringr::str_split( ";", simplify = T)  %>% 
    .[, 5] %>% 
    stringr::str_replace('gene_name', "") %>% 
    stringr::str_replace_all('\\"', "") %>% 
    stringr::str_trim()
}

gencodev19 %>% 
  dplyr::select(chrom = X1, start = X4, end = X5, X9) %>% 
  dplyr::filter(stringr::str_detect(X9, "protein_coding")) %>% 
  dplyr::mutate(symbol = purrr::map_chr(.x = X9, .f = fn_extract_symbol)) %>% 
  dplyr::mutate(gene_length = end - start) %>% 
  dplyr::select(symbol, gene_length) -> symbol_and_length

common_symbol <- intersect(symbol_and_length$symbol, all_protein_coding_genes)

symbol_and_length %>% 
  dplyr::filter(symbol %in% common_symbol) %>% 
  ggplot(aes(x = gene_length)) +
  geom_histogram(bins = 50)

symbol_and_length %>% 
  dplyr::filter(symbol %in% common_symbol) %>% 
  dplyr::mutate(bins = cut(gene_length, 50)) -> symbol_length_cut

symbol_length_cut %>% 
  dplyr::filter(symbol %in% atg_gene) %>% 
  dplyr::group_by(bins) %>% 
  dplyr::count()

fn_seed <- function(s){ 
  s %>% 
    charToRaw() %>% 
    stringr::str_extract_all(pattern = '\\d') %>% 
    unlist() %>% 
    stringr::str_c(collapse = '') %>% 
    as.numeric()
}

'liucj' %>% fn_seed %>% set.seed()

fn_sample <- function(.d){
  n <- .d$n %>% unique()
  .d %>% dplyr::sample_n(size = n)
}
  
fn_random_gene <- function(gene = gene, symbol_length_cut = symbol_length_cut){
  symbol_length_cut %>% dplyr::filter(symbol %in% gene) -> gene_cut
  gene_cut %>% dplyr::group_by(bins) %>% dplyr::count() %>% dplyr::arrange(bins) -> .te
  
  symbol_length_cut %>% 
    dplyr::filter(! symbol %in% gene) %>% 
    dplyr::inner_join(.te, by = "bins") %>% 
    dplyr::group_by(bins) %>% 
    purrrlyr::by_slice(..f = fn_sample, .collate = 'rows') %>% 
    dplyr::select(-n) %>% 
    dplyr::group_by(bins) %>% 
    dplyr::ungroup() %>% 
    dplyr::pull(symbol)
}

fn_mut_count <- function(.x , .d, gene = gene, symbol_length_cut = symbol_length_cut){
  random_symbol <- fn_random_gene(gene = gene, symbol_length_cut = symbol_length_cut)
  .d %>% dplyr::filter(symbol %in% random_symbol) %>% .$count %>% mean()
}

fn_hypermutated <- function(data, n = times, gene = atg_gene, symbol_length_cut = symbol_length_cut){

  data %>% tidyr::gather(key = symbol, value = count, -sample) -> .d
  .d %>% dplyr::filter(symbol %in% gene) %>% .$count %>% mean -> gene_val

  random_val <-
    seq(1: n) %>% 
    purrr::map_dbl(.f = fn_mut_count, .d = .d, gene = gene, symbol_length_cut = symbol_length_cut) 
  
  sum(gene_val > random_val) / length(random_val) -> p_val

}

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
sample_info_mutation_load_atg %>% 
  multidplyr::partition(cluster = cluster) %>% 
  multidplyr::cluster_library('magrittr') %>% 
  multidplyr::cluster_assign_value("fn_mut_count", fn_mut_count) %>%
  multidplyr::cluster_assign_value("fn_random_gene", fn_random_gene) %>%
  multidplyr::cluster_assign_value("fn_sample", fn_sample) %>%
  multidplyr::cluster_assign_value("fn_hypermutated", fn_hypermutated) %>%
  multidplyr::cluster_assign_value("atg_gene", atg_gene) %>%
  multidplyr::cluster_assign_value("symbol_length_cut", symbol_length_cut) %>%
  dplyr::mutate(hyper = purrr::map_dbl(.x = data, .f = fn_hypermutated, n = 1000, gene = atg_gene, symbol_length_cut = symbol_length_cut)) %>% 
  dplyr::collect() %>% 
  dplyr::as_tibble() %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-PARTITION_ID) ->
  atg_mutation_load_hypermutation

sample_info_mutation_load_lys %>% 
  multidplyr::partition(cluster = cluster) %>% 
  multidplyr::cluster_library('magrittr') %>% 
  multidplyr::cluster_assign_value("fn_mut_count", fn_mut_count) %>%
  multidplyr::cluster_assign_value("fn_random_gene", fn_random_gene) %>%
  multidplyr::cluster_assign_value("fn_sample", fn_sample) %>%
  multidplyr::cluster_assign_value("fn_hypermutated", fn_hypermutated) %>%
  multidplyr::cluster_assign_value("lys_gene", lys_gene) %>%
  multidplyr::cluster_assign_value("symbol_length_cut", symbol_length_cut) %>%
  dplyr::mutate(hyper = purrr::map_dbl(.x = data, .f = fn_hypermutated, n = 1000, gene = lys_gene, symbol_length_cut = symbol_length_cut)) %>% 
  dplyr::collect() %>% 
  dplyr::as_tibble() %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-PARTITION_ID) ->
  lys_mutation_load_hypermutation
parallel::stopCluster(cluster)

atg_mutation_load_hypermutation %>% 
  dplyr::select(-data) -> atg_mutation_load_plot_ready
atg_mutation_load_plot_ready %>%   
  readr::write_rds(path = file.path(snv_path, '.rds_02_snv_c_atg_mutation_load_hypermutation.rds.gz'), compress = 'gz')

CPCOLS <- c("#FF3030", "#C7C7C7", "#191970")

atg_mutation_load_plot_ready %>% 
  dplyr::mutate(mut_type = dplyr::case_when(hyper < 0.05 ~ 'Hypomutated', hyper > 0.95 ~ 'Hypermutated', TRUE ~ 'Non-significant')) %>% 
  ggplot(aes(x = reorder(cancer_types, - mutation_load), y = mutation_load, fill = mut_type)) +
  geom_bar(stat = 'identity') +
  # ggsci::scale_fill_npg(name = '') +
  scale_fill_manual(
    name = "Mutation Type",
    limits = c("Hypermutated", "Non-significant", "Hypomutated"),
    values = CPCOLS
  ) +
  # ggthemes::theme_gdocs() +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = 'Cancer Type', y = 'Mutation Load') -> p
ggsave(filename = 'atg_mutation_load.pdf', plot = p, device = 'pdf', path = snv_path, width = 9, height = 4)

lys_mutation_load_hypermutation %>% 
  dplyr::select(-data) -> lys_mutation_load_plot_ready
lys_mutation_load_plot_ready %>% 
  readr::write_rds(path = file.path(snv_path, '.rds_02_snv_c_lys_mutation_load_hypermutation.rds.gz'), compress = 'gz')

lys_mutation_load_plot_ready %>% 
  dplyr::mutate(mut_type = dplyr::case_when(hyper < 0.05 ~ 'Hypomutated', hyper > 0.95 ~ 'Hypermutated', TRUE ~ 'Non-significant')) %>% 
  ggplot(aes(x = reorder(cancer_types, - mutation_load), y = mutation_load, fill = mut_type)) +
  geom_bar(stat = 'identity') +
  # ggthemes::theme_gdocs() +
  # ggsci::scale_fill_npg(name = '') +
  scale_fill_manual(
    name = "Mutation Type",
    limits = c("Hypermutated", "Non-significant", "Hypomutated"),
    values = CPCOLS
  ) +
  # ggthemes::theme_gdocs() +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = 'Cancer Type', y = 'Mutation Load') ->p
ggsave(filename = 'lys_mutation_load.pdf', plot = p, device = 'pdf', path = snv_path, width = 9, height = 4)

CPCOLS <- c("#FF3030", "#C7C7C7", "#191970")
atg_mutation_load_plot_ready %>% 
  dplyr::bind_rows(lys_mutation_load_plot_ready) %>% 
  dplyr::mutate(gene_type = rep(c('Autophagy Genes', 'Lysosome Genes'), times = c(33, 28))) %>% 
  dplyr::mutate(mut_type = dplyr::case_when(hyper < 0.05 ~ 'Hypomutated', hyper > 0.95 ~ 'Hypermutated', TRUE ~ 'Non-significant')) %>% 
  dplyr::mutate(cancer_types = factor(cancer_types)) %>% 
  ggplot(aes(x = cancer_types, y = mutation_load, fill = mut_type)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(
    name = "Mutation Type",
    limits = c("Hypermutated", "Non-significant", "Hypomutated"),
    values = CPCOLS
  ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
    
    strip.background = element_rect(color = "black", fill = "white")
  ) +
  labs(x = 'Cancer Type', y = 'Mutation Load') +
  facet_grid(gene_type ~ .) -> p
ggsave(filename = 'atg_lys_combined_mutation_load.pdf', plot = p, device = 'pdf', path = snv_path, width = 7, height = 4)

#
save.image(file = file.path(snv_path, ".rda_02_snv_c_mutation_load.rda"))
load(file = file.path(snv_path, ".rda_02_snv_c_mutation_load.rda"))
