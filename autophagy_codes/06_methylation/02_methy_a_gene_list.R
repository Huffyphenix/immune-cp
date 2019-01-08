library(ggplot2)
library(magrittr)


# Path
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")
methy_path <- file.path(out_path, "m_4_methylation")
methy_box <- file.path(methy_path, "boxplot")


# load methylation and gene list
methy <- readr::read_rds(file.path(tcga_path, "pancan33_meth.rds.gz"))

#load gene list
gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
gene_list$col<-as.character(gene_list$col)
gene_list$size<-as.character(gene_list$size)


# functions
filter_gene_list <- function(.x, gene_list) {
  gene_list %>%
    dplyr::select(symbol) %>%
    dplyr::left_join(.x, by = "symbol")
}

methy %>% 
  # dplyr::slice(2:5) %>%  # tidyr::unnest()
  dplyr::mutate(filter_methy = purrr::map(methy, filter_gene_list, gene_list = gene_list)) %>% 
  dplyr::select(-methy) -> gene_list_methy

gene_list_methy %>% readr::write_rds(file.path(methy_path, ".rds_02_methy_a_gene_list_methy.rds.gz"), compress = 'gz')
#------------------------------------------------------------
fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
}
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
}
fun_boxplot <- function(fig_name, data, path = methy_box){
  data %>% 
    ggpubr::ggboxplot(x = 'type', y = 'meth', color = 'type', pallete = 'jco' ) +
    ggpubr::stat_compare_means(
      method = "t.test", 
      label.y = 1, 
      label.x = 1.2) + 
    ggthemes::theme_gdocs() +
    scale_color_manual(values = c("#DC3912", "#3366CC")) +
    labs(y = "B-value", x = "", title = fig_name) -> p
  ggsave(filename = paste(fig_name, "pdf", sep = "."), p, device = "pdf", path = path, width = 4, height = 3)
}
fun_compare <- function(.x, .y ){
  .x %>% 
    dplyr::mutate(gene = as.character(gene )) %>% 
    tidyr::gather(key = barcode, value = meth, -symbol, -gene) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type %in% c("01", "11")) %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(symbol, gene, barcode, meth, type) %>% 
    dplyr::mutate(type = dplyr::case_when(
      type == "01" ~ "Tumor",
      type == "11" ~ "Normal"
    )) %>% 
    dplyr::filter(!is.na(gene)) %>% 
    tidyr::drop_na() -> .d
  if(nrow(.d) < 20 || length(unique(.d$type)) != 2){return(tibble::tibble())}
  # at least 10 samples
  .d %>% 
    dplyr::select(barcode, type) %>% 
    dplyr::distinct() %>%
    dplyr::group_by(type) %>% 
    dplyr::count() %>% 
    dplyr::pull(n) -> sample_num
  if(any(sample_num < 10)){return(tibble::tibble())}
  
  .d %>% 
    dplyr::group_by(symbol, gene) %>% 
    dplyr::do(
      broom::tidy(
        t.test(meth ~ type, data = .)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    dplyr::filter(fdr < 0.05) %>%
    dplyr::mutate(
      direction = dplyr::case_when(
          estimate > 0 ~ 0, # normal high
          estimate < 0 ~ 1 # tumor high
        )) %>% 
    dplyr::select(symbol, gene, direction, p_val = p.value, fdr) -> .d_out
  
  .d %>% 
    dplyr::group_by(symbol, type) %>% 
    dplyr::summarise(m = mean(meth)) %>% 
    tidyr::spread(key = type, value = m) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(diff = Tumor - Normal) %>% 
    dplyr::select(symbol, diff) %>% 
    dplyr::inner_join(.d_out, by = "symbol") -> .d_out_diff
  
  # draw every pic
  # .d %>%
  #   dplyr::semi_join(.d_out, by = c("symbol", "gene")) %>%
  #   tidyr::nest(-symbol, -gene) %>%
  #   dplyr::mutate(fig_name = paste(.y, symbol, sep = "_")) %>% 
  #   dplyr::select(fig_name, data) %>% 
  #   purrr::pwalk(.f = fun_boxplot, path = methy_box)
  
  return(.d_out_diff)
}


gene_list_methy %>%
  dplyr::slice(1) %>%
  dplyr::mutate(methy_comparison = purrr::map(.x = filter_methy, .y = cancer_types, .f = fun_compare)) %>%
  dplyr::select(-filter_methy) %>%
  tidyr::unnest()

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
gene_list_methy %>% 
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_library("ggplot2") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type)  %>%
  multidplyr::cluster_assign_value("fun_boxplot", fun_boxplot)  %>%
  multidplyr::cluster_assign_value("fun_compare", fun_compare)  %>%
  multidplyr::cluster_assign_value("methy_box", methy_box)  %>%
  dplyr::mutate(methy_comparison = purrr::map(.x = filter_methy, .y = cancer_types, .f = fun_compare)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  dplyr::select(-filter_methy) %>% 
  tidyr::unnest() -> gene_list_methy_fdr
parallel::stopCluster(cluster)
readr::write_rds(gene_list_methy_fdr, path = file.path(methy_path, ".rds_02_gene_list_methy_fdr.rds.gz"), compress = "gz")

gene_list_methy_fdr <- readr::read_rds(path = file.path(methy_path, ".rds_02_gene_list_methy_fdr.rds.gz"))

gene_list_methy_fdr %>% 
  dplyr::filter(fdr < 0.05) %>% 
  dplyr::mutate(fdr = -log10(fdr)) %>% 
  dplyr::mutate(fdr = ifelse(fdr > 30, 30, fdr)) %>% 
  dplyr::mutate(diff = dplyr::case_when(diff < -0.2 ~ -0.2, diff > 0.2 ~ 0.2, TRUE ~ diff)) -> plot_ready
  
plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::count() %>% 
  dplyr::arrange(dplyr::desc(n)) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(v = sum(diff)) %>% 
  #dplyr::filter(symbol %in% atg_gene) %>% 
  dplyr::left_join(gene_list, by = 'symbol') %>%
  # dplyr::filter(pathway == "autophagesome formation-core") %>%
  # dplyr::filter(status == "l") %>% 
  #dplyr::mutate(color = ifelse(is.na(marker), 'black', 'red')) %>%
  dplyr::arrange(v) -> gene_rank

CPCOLS <- c("firebrick1", "gray97", "midnightblue")

plot_ready %>% 
  ggplot(aes(x = symbol, y = cancer_types)) +
  geom_point(aes(size = fdr, color = diff)) +
  scale_y_discrete(limit = cancer_rank$cancer_types) +
  scale_x_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(name = "FDR") +
  scale_color_gradient2(
    name = "Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1]
  ) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    # axis.text.y = element_text(color = gene_rank$color),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = gene_rank$col,face = gene_rank$size),
    
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    color = guide_colorbar(
      title.position = "bottom",
      title.hjust = 0.5,
      barheight = 0.5,
      barwidth = 10
    )
  ) -> p

ggsave(filename = '02_atg_meth_diff.pdf', device = 'pdf', plot = p, path = methy_path, height = 5.5, width = 15)
#------------------------------lysosome----------------

plot_ready %>% 
  dplyr::group_by(cancer_types) %>% 
  dplyr::count() %>% 
  dplyr::arrange(dplyr::desc(n)) -> cancer_rank

plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(v = sum(diff)) %>% 
  dplyr::filter(symbol %in% lys_gene) %>% 
  dplyr::left_join(gene_list, by = 'symbol') %>%
  # dplyr::filter(pathway == "autophagesome formation-core") %>%
  # dplyr::filter(status == "l") %>% 
  # dplyr::mutate(color = ifelse(is.na(marker), 'black', 'red')) %>%
  dplyr::arrange(v) -> gene_rank

CPCOLS <- c("firebrick1", "gray97", "midnightblue")

plot_ready %>% 
  ggplot(aes(x = symbol, y = cancer_types)) +
  geom_point(aes(size = fdr, color = diff)) +
  scale_y_discrete(limit = cancer_rank$cancer_types) +
  scale_x_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(name = "FDR") +
  scale_color_gradient2(
    name = "Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1]
  ) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    # axis.text.y = element_text(color = gene_rank$color),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = gene_rank$color),
    
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key = element_rect(fill = "white", colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    color = guide_colorbar(
      title.position = "bottom",
      title.hjust = 0.5,
      barheight = 0.5,
      barwidth = 10
    )
  ) -> p

ggsave(filename = '02_meth_diff.pdf', device = 'pdf', plot = p, path = methy_path, height = 5.5, width = 8)


#
save.image(file = file.path(methy_path, ".rda_02_methy_a_gene_list.rda"))
load(file = file.path(methy_path, ".rda_02_methy_a_gene_list.rda"))

