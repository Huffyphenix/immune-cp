library(magrittr)
library(ggplot2)
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")

# tcga_path = "/home/cliu18/liucj/projects/6.autophagy/TCGA"
# expr_path <- "/home/cliu18/liucj/projects/6.autophagy/02_autophagy_expr/"
# expr_path_a <- file.path(expr_path, "03_a_gene_expr")
# snv_path <- "/home/cliu18/liucj/projects/6.autophagy/04_snv"

# load cnv and gene list
snv_syn7824274 <- readr::read_rds(file.path("/data/shiny-data/GSCALite/TCGA/snv","pancan33_snv_from_syn7824274_gather.rds.gz"))
snv <- readr::read_rds(file.path(tcga_path, "pancan33_snv.rds.gz"))
#marker_file <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_marker.rds.gz"))
#gene_list <- readr::read_rds(file.path(expr_path_a, "rds_03_a_atg_lys_gene_list.rds.gz")) %>% 
#  dplyr::left_join(marker_file, by = "symbol") %>% 
#  dplyr::mutate(symbol = dplyr::recode(symbol, "ATG101" = "C12orf44"))

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)


filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

snv %>%
  dplyr::mutate(filter_snv = purrr::map(snv, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-snv) -> gene_list_snv
readr::write_rds(x = gene_list_snv, path = file.path(snv_path, ".rds_02_snv_a_gene_list.rds.gz"), compress = "gz")

snv_syn7824274 %>%
  dplyr::filter(symbol %in% gene_list$symbol) -> gene_list_snv_syn7824274
fn_merge_data <- function(.data){
  .data %>%
    tidyr::gather(-symbol,key = sample, value = count) %>%
    dplyr::left_join(gene_list_snv_syn7824274,by=c("sample","symbol")) %>%
    dplyr::select(symbol,sample,mut_n) %>%
    tidyr::spread(key = sample, value = mut_n) 
}
gene_list_snv %>%
  dplyr::mutate(filter_snv_syn7824274 = purrr::map(filter_snv,fn_merge_data)) %>%
  dplyr::select(-filter_snv) -> gene_list_snv_syn7824274.with_full_samples
readr::write_rds(x = gene_list_snv_syn7824274.with_full_samples, path = file.path(snv_path, ".rds_02_snv_a_gene_list_syn7824274.with_full_samples.rds.gz"), compress = "gz")

# filter out hypermutation samples -----
burden_path <- "/project/huff/huff/data/TCGA"
mutation_burden_class <- readr::read_rds(file.path(burden_path,"classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest()

# gene_list_snv %>%
#   tidyr::unnest() %>%
#   tidyr::gather(-cancer_types,-symbol,key="barcode",value="count") %>%
#   dplyr::left_join(mutation_burden_class,by="barcode") %>%
#   dplyr::filter(!sm_count > 1000) %>%
#   dplyr::select(-Cancer_Types) -> gene_list_snv.filter_hypermutation 

fn_get_percent <- function(cancer_types, filter_snv){
  print(cancer_types)
  n <- length(filter_snv) - 1
  filter_snv %>%
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .))) -> .d
  .d %>% 
    tidyr::gather(-symbol,key = barcode, value = count) %>% 
    dplyr::left_join(mutation_burden_class,by="barcode") %>%
    dplyr::filter(!sm_count > 1000) %>%
    dplyr::mutate(samples = ifelse(count > 0, 1, 0)) %>% 
    dplyr::group_by(symbol) %>% 
    dplyr::summarise(sm_count = sum(count), sm_sample = sum(samples)) %>% 
    dplyr::mutate(per = sm_sample / n) %>%
    dplyr::mutate(n=n)-> .d_count

  return(.d_count)
    # tibble::tibble( mut_count = list(.d_count))
}
fn_get_hyper_mut_data_class <- function(cancer_types, filter_snv){
  print(cancer_types)
  n <- length(filter_snv) - 1
  filter_snv %>%
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .))) -> .d
  .d %>% 
    tidyr::gather(key = barcode, value = count, -symbol) %>% 
    dplyr::left_join(mutation_burden_class,by="barcode") %>%
    dplyr::mutate(class = ifelse(sm_count > 1000,"hyper","non-hyper")) 
}
gene_list_snv_syn7824274.with_full_samples %>% 
  # head(1) %>% .$filter_snv_syn7824274 %>% .[[1]] -> filter_snv
  plyr::mutate(res = purrr::map2(cancer_types, filter_snv_syn7824274, fn_get_percent)) %>%
  dplyr::select(-filter_snv_syn7824274) -> gene_list_snv.count_per.filter_hypermutation

gene_list_snv_syn7824274.with_full_samples %>% 
  # head(1) %>% .$filter_snv %>% .[[1]] -> filter_snv
  plyr::mutate(res = purrr::map2(cancer_types, filter_snv_syn7824274, fn_get_hyper_mut_data_class)) %>%
  dplyr::select(-filter_snv_syn7824274) -> gene_list_snv.hypermutation_class

gene_list_snv.hypermutation_class %>%
  readr::write_rds(file.path(snv_path,"gene_list_snv_syn7824274.hypermutation_class"),compress = "gz")
# gene_list_snv.filter_hypermutation %>% 
#   dplyr::select(cancer_types,barcode,symbol,count) %>%
#   tidyr::nest(-cancer_types,.key="filter_snv") %>%
#   head(1) %>% .$filter_snv %>% .[[1]] -> filter_snv
#   plyr::mutate(res = purrr::map2(cancer_types, filter_snv, fn_get_percent)) %>%
#   dplyr::select(-filter_snv) %>%
#   dplyr::ungroup() -> gene_list_snv.count_per

# cl <- parallel::detectCores()
# cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
# gene_list_snv %>% 
#   multidplyr::partition(cluster = cluster) %>%
#   multidplyr::cluster_library("magrittr") %>%
#   multidplyr::cluster_assign_value("fn_get_percent", fn_get_percent)  %>%
#   dplyr::mutate(res = purrr::map2(cancer_types, filter_snv, fn_get_percent)) %>% 
#   dplyr::collect() %>%
#   dplyr::as_tibble() %>%
#   dplyr::ungroup() %>%
#   dplyr::select(-PARTITION_ID) %>% 
#   dplyr::select(-cancer_types, -filter_snv) %>% 
#   tidyr::unnest(res) -> gene_list_snv_count
# parallel::stopCluster(cluster)

gene_list_snv.count_per.filter_hypermutation %>% 
  readr::write_rds(path = file.path(snv_path, "syn7824274_snv_gene_list_snv_count_per.filter_hypermutation.rds.gz"), compress = "gz")

gene_list_snv.count_per.filter_hypermutation %>% 
  tidyr::unnest() %>%
  tidyr::drop_na() %>%
  dplyr::mutate(x_label = paste(cancer_types, " (n=", n,")", sep = "")) %>% 
  dplyr::mutate(sm_count = ifelse(sm_count>0, sm_count, NA)) -> plot_ready



plot_ready %>% 
  dplyr::group_by(x_label) %>% 
  dplyr::summarise(s = sum(per)) %>% 
  dplyr::arrange(dplyr::desc(s)) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(s = sum(sm_sample)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(functionWithImmune, replace = c('TwoSide' = "#EE7600", "Inhibit" = "#EE3B3B", "Activate" = "#1C86EE"))) %>% 
  dplyr::mutate(size = plyr::revalue(type,replace = c('Receptor'="bold",'Ligand'="plain",'Ligand&Receptor'="bold.italic"))) %>%
  dplyr::arrange(color,s) -> gene_rank
gene_rank$color %>% as.character() ->gene_rank$color
gene_rank$size %>% as.character() ->gene_rank$size


plot_ready %>% 
  # dplyr::filter(!symbol %in% c("TP53", "PTEN", "CDKN2A")) %>% 
  # dplyr::mutate(per = ifelse(per > 0.07, 0.07, per)) %>% 
  # dplyr::filter(per > 0.02) %>% 
  ggplot(aes(x = x_label, y = symbol, fill = per)) +
  geom_tile() +
  geom_text(aes(label = sm_count)) +
  scale_x_discrete(position = "top", limits = cancer_rank$x_label) +
  scale_y_discrete(limits = gene_rank$symbol) +
  scale_fill_gradient2(
    name = "Mutation Frequency (%)",
    limit = c(0, 0.10),
    breaks = seq(0, 0.1, 0.01),
    label = c("0", "","2","","4","","6","","8", "","10"),
    high = "red",
    na.value = "white"
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = -0.05),
    axis.text.y = element_text(color = gene_rank$color,face = gene_rank$size)
  ) +
  guides(fill = guide_legend(title = "Mutation Frequency (%)", 
                             title.position = "left", 
                             title.theme = element_text(angle = 90, vjust = 2), 
                             reverse = T, 
                             keywidth = 0.6, 
                             keyheight = 0.8 )) +
  labs(x = "", y = "") -> p;p

# ggsave(filename = "core_atg_freq.pdf", plot = p, device = "pdf", path = snv_path, width = 9, height = 7)
ggsave(filename = "ICP_mutation_freq.filter_hypermutation.pdf", plot = p, device = "pdf", path = snv_path, width = 9, height = 12)
ggsave(filename = "ICP_mutation_freq.filter_hypermutation.png", plot = p, device = "png", path = snv_path, width = 9, height = 12)


save.image(file = file.path(snv_path, ".rda_02_snv_a_gene_list.rda"))
load(file = file.path(snv_path, ".rda_02_snv_a_gene_list.rda"))









