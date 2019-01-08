library(magrittr)

out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
cnv_path <- file.path(out_path, "m_1_cnv")
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")
# load cnv and gene list
cnv <- readr::read_rds(file.path(tcga_path, "pancan34_cnv.rds.gz"))
gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

cnv %>%
  dplyr::mutate(filter_cnv = purrr::map(cnv, filter_gene_list, gene_list = gene_list)) %>%
  dplyr::select(-cnv) -> gene_list_cnv

readr::write_rds(x = gene_list_cnv, path = file.path(cnv_path, ".rds_02_cnv_a_gene_list.rds.gz"), compress = "gz")

fn_get_amplitue_threshold <- function(.x){
  ifelse(abs(.x) < log2(4 / 2) , 0, .x) -> .y #copy num variation A>1倍 OR D>-1倍
  tibble::tibble(a = sum(.y > 0) / length(.y), d = sum(.y < 0) / length(.y)) 
}
fn_get_ad <- function(.d){
  .d %>% 
    unlist(use.name = F) %>% 
    fn_get_amplitue_threshold()
}
fn_get_percent <- function(cancer_types, filter_cnv){
  filter_cnv %>%
    tidyr::nest(-symbol) %>% 
    dplyr::mutate(ad = purrr::map(data, .f = fn_get_ad)) %>% 
    dplyr::select(-data) %>% 
    tidyr::unnest(ad) %>% 
    tibble::add_column(cancer_types = cancer_types, .before = 1)
}

gene_list_cnv %>% dplyr::filter(cancer_types == "KIRC") %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>% 
  dplyr::select(-cancer_types, -filter_cnv) %>% 
  tidyr::unnest(rs)

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_cnv %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_amplitue_threshold", fn_get_amplitue_threshold)  %>%
  multidplyr::cluster_assign_value("fn_get_ad", fn_get_ad) %>% 
  multidplyr::cluster_assign_value("fn_get_percent", fn_get_percent) %>% 
  dplyr::mutate(rs = purrr::map2(cancer_types, filter_cnv, fn_get_percent)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-cancer_types, -filter_cnv) %>% 
  tidyr::unnest(rs) -> gene_list_cnv_per
on.exit(parallel::stopCluster(cluster))

library(ggplot2)
gene_list_cnv_per %>% 
  tidyr::drop_na() %>% 
  tidyr::gather(key = type, value = per, a, d) %>% 
  dplyr::mutate(type = plyr::revalue(type, replace = c("a" = "Amplification", "d" = "Deletion"))) %>% 
  dplyr::mutate(per = ifelse(per > 0.6, 0.6, per)) -> plot_ready
plot_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = type)) +
  scale_size_continuous(
    name = "CNV perl",
    breaks = c(0.1, 0.2, 0.4, 0.6),
    limits = c(0.1, 0.6),
    labels = c("10", "20", "40", "60")
  ) +
  ggthemes::scale_color_gdocs(
    name = "SCNA Type"
  ) +
  facet_wrap(~ type) -> p
ggsave(filename = "01_SCNV_all.pdf", plot = p, device = "pdf", path = cnv_path, width = 25, height = 30)




plot_ready %>% dplyr::filter(type == "Amplification") -> a_ready
a_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(s = sum(per)) %>% 
  dplyr::arrange(dplyr::desc(s)) -> cancer_rank
a_ready %>% 
  dplyr::filter(per > 0.1) %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(s = n()) %>% 
  dplyr::filter(s >= 5) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  # dplyr::filter(status %in% c("l")) %>% 
  # dplyr::filter(s > 1.8) %>% 
  dplyr::arrange(col, s) -> a_gene_rank

a_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = type)) +
  scale_size_continuous(
    name = "CNV percentage",
    breaks = c(0.1, 0.2, 0.4, 0.6),
    limits = c(0.1, 0.6),
    labels = c("10", "20", "40", "60")
  ) +
  scale_y_discrete(limit = a_gene_rank$symbol) +
  scale_x_discrete(limit = cancer_rank$cancer_types)+
  ggthemes::scale_color_gdocs(
    name = "SCNA Type"
  ) +
  ggthemes::theme_gdocs() +
  theme(axis.text.y = element_text(color = a_gene_rank$col)) -> p;p
ggsave(filename = "02_SCNV_amplification_seminar.pdf", plot = p, device = "pdf", path = cnv_path, width = 15, height = 9)


plot_ready %>% dplyr::filter(type == "Deletion") -> d_ready
d_ready %>% 
  dplyr::filter(per >= 0.1) %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::summarise(s = sum(per)) %>% 
  dplyr::arrange(dplyr::desc(s)) -> cancer_rank
d_ready %>%   
  dplyr::filter(per >= 0.1) %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(s = n()) %>% 
  dplyr::filter(s >= 5) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(status, replace = c('a' = "#e41a1c", "l" = "#377eb8", "i" = "#4daf4a", "p" = "#984ea3"))) %>% 
  # dplyr::filter(s > 5) %>%
  dplyr::arrange(status, s) -> d_gene_rank

d_ready %>% 
  ggplot(aes(y = symbol, x = cancer_types)) +
  geom_point(aes(size = per, color = type)) +
  scale_size_continuous(
    name = "CNV percentage",
    breaks = c(0.1, 0.2, 0.4, 0.6),
    limits = c(0.1, 0.6),
    labels = c("10", "20", "40", "60")
  ) +
  scale_y_discrete(limit = d_gene_rank$symbol) +
  scale_x_discrete(limit = cancer_rank$cancer_types)+
  ggthemes::scale_color_gdocs(
    name = "SCNA Type"
  ) +
  ggthemes::theme_gdocs() +
  theme(axis.text.y = element_text(color = d_gene_rank$color)) -> p
ggsave(filename = "03_SCNV_deletion_seminar.pdf", plot = p, device = "pdf", path = cnv_path, width = 15, height = 12)


save.image(file = file.path(cnv_path, ".rda_02_cnv_a_gene_list.rda"))
load(file = file.path(cnv_path, ".rda_02_cnv_a_gene_list.rda"))




