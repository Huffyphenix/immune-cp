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


fn_get_percent <- function(cancer_types, filter_snv){
  print(cancer_types)
  n <- length(filter_snv) - 1
  filter_snv %>%
    dplyr::mutate_if(.predicate = is.numeric, .fun = dplyr::funs(ifelse(is.na(.), 0, .))) -> .d
    .d %>% 
      tidyr::gather(key = barcode, value = count, -symbol) %>% 
      dplyr::mutate(samples = ifelse(count > 0, 1, 0)) %>% 
      dplyr::group_by(symbol) %>% 
      dplyr::summarise(sm_count = sum(count), sm_sample = sum(samples)) %>% 
      dplyr::mutate(per = sm_sample / n) -> .d_count

    tibble::tibble(cancer_types = cancer_types, n = n, mut_count = list(.d_count))
}
gene_list_snv %>% 
  head(1) %>% #.$filter_snv %>% .[[1]] -> filter_snv
  plyr::mutate(res = purrr::map2(cancer_types, filter_snv, fn_get_percent))

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_snv %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fn_get_percent", fn_get_percent)  %>%
  dplyr::mutate(res = purrr::map2(cancer_types, filter_snv, fn_get_percent)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-cancer_types, -filter_snv) %>% 
  tidyr::unnest(res) -> gene_list_snv_count
parallel::stopCluster(cluster)

gene_list_snv_count %>% 
  readr::write_rds(path = file.path(snv_path, ".rds_snv_a_gene_list_snv_count.rds.gz"), compress = "gz")

gene_list_snv_count %>% 
  tidyr::unnest(mut_count) %>% 
  tidyr::drop_na() %>%
  dplyr::mutate(x_label = paste(cancer_types, " (n=", n,")", sep = "")) %>% 
  dplyr::mutate(sm_count = ifelse(sm_count>0, sm_count, NA)) -> plot_ready

# core atg
# gene_list %>%
#   dplyr::filter(pathway == "autophagesome formation-core") -> gene_list_core
# plot_ready %>%
#   dplyr::semi_join(gene_list_core, by = "symbol") -> plot_ready

# lys
gene_list %>%
  dplyr::filter(status == "l") -> gene_list_lys

plot_ready %>%
  dplyr::semi_join(gene_list_lys, by = "symbol") -> plot_ready


plot_ready %>% 
  dplyr::group_by(x_label) %>% 
  dplyr::summarise(s = sum(per)) %>% 
  dplyr::arrange(dplyr::desc(s)) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(s = sum(sm_sample)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  # dplyr::filter(status %in% c("p", "i")) %>% 
  # dplyr::mutate(color = plyr::revalue(status, replace = c('a' = "#e41a1c", "l" = "#377eb8", "i" = "#4daf4a", "p" = "#984ea3"))) %>%
  # dplyr::mutate(color = ifelse(is.na(marker), "black", "red")) %>%
  dplyr::arrange(col,s) -> gene_rank
gene_rank$size %>% as.character() ->gene_rank$size
gene_rank$col %>% as.character() ->gene_rank$col


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
    limit = c(0, 0.15),
    breaks = seq(0, 0.15, 0.01),
    label = c("0", "","","3","","","6","","", "9","","", "12","","","15"),
    high = "red",
    na.value = "white"
  ) +
  ggthemes::theme_gdocs() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = -0.05),
    axis.text.y = element_text(color = gene_rank$col,face = gene_rank$size)
  ) +
  guides(fill = guide_legend(title = "Mutation Frequency (%)", 
                             title.position = "left", 
                             title.theme = element_text(angle = 90, vjust = 2), 
                             reverse = T, 
                             keywidth = 0.6, 
                             keyheight = 0.8 )) +
  labs(x = "", y = "") -> p;p

# ggsave(filename = "core_atg_freq.pdf", plot = p, device = "pdf", path = snv_path, width = 9, height = 7)
ggsave(filename = "lys_freq.pdf", plot = p, device = "pdf", path = snv_path, width = 9, height = 9)



save.image(file = file.path(snv_path, ".rda_02_snv_a_gene_list.rda"))
load(file = file.path(snv_path, ".rda_02_snv_a_gene_list.rda"))









