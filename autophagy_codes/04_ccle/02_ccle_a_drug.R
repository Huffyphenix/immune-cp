library(magrittr)
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")
ccle_path <- file.path(out_path, "d_1_drug")

# load data
drug_gdsc <- readr::read_rds(file.path(tcga_path, "drug_gdsc_exp_spearman.rds.gz"))
drug_ctrp <- readr::read_rds(file.path(tcga_path, "drug_ctrp_exp_spearman.rds.gz"))

#load gene list
gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

t_gdsc <- readr::read_rds(file.path(tcga_path, "drug_target_gdsc.rds.gz")) %>% 
  tidyr::unnest() %>% 
  dplyr::select(drug_name, target_pathway) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(target_pathway) %>% 
  dplyr::mutate(count = n()) %>% 
  dplyr::ungroup()
t_ctrp <- readr::read_rds(file.path(tcga_path, "drug_target_ctrp.rds.gz")) %>% 
  tidyr::unnest() %>% 
  dplyr::select(drug_name, target_pathway) %>% 
  dplyr::distinct()

# GDSC analysis
gene_list %>% 
  dplyr::inner_join(drug_gdsc, by = "symbol") -> gdsc_gene_list

fn_filter_drug <- function(.x, .cor = 0.25, .fdr = 0.05) {
  .x %>% dplyr::filter(abs(cor_sprm) > .cor, fdr < .fdr)
}

gdsc_gene_list %>% 
  # dplyr::filter(symbol == "ATG7") %>% 
  # .$drug %>% .[[1]] -> .x
  dplyr::mutate(cor_drug = purrr::map(.x = drug, .f = fn_filter_drug)) %>% 
  tidyr::unnest(cor_drug) -> gdsc_gene_list_sig_drug
gdsc_gene_list_sig_drug %>%
  write.table(file = file.path(ccle_path,"gdsc_gene_list_sig_drug.txt"),quote = F)

gdsc_gene_list_sig_drug %>%
  dplyr::filter(symbol=="TNFRSF14") %>%
  dplyr::inner_join(t_gdsc,by="drug_name") %>%
  write.table(file = file.path(ccle_path,"gdsc_TNFRSF14_drug.txt"),quote = F,sep="\t")

# gene rank
library(ggplot2)
gdsc_gene_list_sig_drug %>% 
  dplyr::mutate(
    fdr = ifelse(-log10(fdr) > 40, 40, -log10(fdr))
    #cor_sprm = ifelse(cor_sprm > 0.4, 0.4, cor_sprm),
    #cor_sprm = ifelse(cor_sprm < -0.4, -0.4, cor_sprm)
  ) %>% 
  dplyr::left_join(t_gdsc, by = "drug_name") -> gdsc_plot_ready


gdsc_plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(cor_sum = sum(cor_sprm)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  #dplyr::mutate(color = ifelse(status == "l", "black", "red")) %>% 
  dplyr::arrange(cor_sum) -> gdsc_gene_rank
gdsc_gene_rank$col <- gdsc_gene_rank$col %>% as.character()
gdsc_gene_rank$size <- gdsc_gene_rank$size %>% as.character()

gdsc_plot_ready %>% 
  dplyr::distinct(drug_name, target_pathway, count) %>% 
  dplyr::group_by(target_pathway) %>% 
  dplyr::mutate(per = n() / count) %>% print(width = Inf) %>% 
  dplyr::arrange(per) %>% 
  dplyr::ungroup() %>% dplyr::select(drug_name, target_pathway, per) -> gdsc_drug_per

gdsc_plot_ready %>% 
  dplyr::left_join(gdsc_drug_per, by = c("drug_name", "target_pathway")) %>% 
  dplyr::group_by(drug_name) %>% 
  dplyr::mutate(drug_count = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(per, target_pathway, drug_count) %>% 
  dplyr::select(drug_name, target_pathway, drug_count, per) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(target_pathway = stringr::str_to_title(target_pathway)) -> gdsc_drug_rank_pre

pathway_color <- gdsc_drug_rank_pre %>% 
  dplyr::distinct(target_pathway, per) %>% 
  dplyr::arrange(per, target_pathway) %>% 
  dplyr::mutate(color = ggthemes::gdocs_pal()(17))

gdsc_drug_rank_pre %>% 
  dplyr::select(-per) %>% 
  dplyr::left_join(pathway_color, by = "target_pathway") -> drug_rank

pathway_color %>%
  ggplot(aes(x = 1, y = per, fill = color)) +
  geom_bar(stat = "identity", position = "stack", width = 0.001) +
  scale_fill_manual(values = pathway_color$color)
  
gdsc_plot_ready %>% 
  ggplot(aes(x = symbol, y = drug_name, color = cor_sprm)) +
  geom_point(aes(size = fdr)) +
  scale_x_discrete(limits = gdsc_gene_rank$symbol, expand = c(0.012,0.012)) +
  scale_y_discrete(limits = drug_rank$drug_name, expand = c(0.012,0.012), position = "right") +
  scale_color_gradient2(
    name = "Spearman Correlation",
    high = "red",
    mid = "white",
    low = "blue"
  ) +
  scale_size_continuous(
    name = "FDR"
  ) +
  # ggthemes::theme_gdocs() +
  theme(
    panel.background = element_rect(color = "black", fill = "white", size = 0.1),
    panel.grid=element_line(colour="grey",linetype="dashed"),
    panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
    
    axis.title = element_blank(),
    axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5, color = gdsc_gene_rank$col),
    axis.text.y = element_text(size = 10, color = drug_rank$color),
    
    axis.ticks = element_line(color = "black"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    # legend.key.width = unit(1,"cm"),
    # legend.key.heigh = unit(0.3,"cm"),
    legend.key = element_rect(fill="white",colour = "black")
  ) + guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = 0.5,
      barwidth = 10
    )
  )->p;p

ggsave(filename = 'gdsc_sig_drug.pdf', plot = p, device = "pdf", path = ccle_path, width = 8, height = 12)

# CTRP analysis
gene_list %>% 
  dplyr::inner_join(drug_ctrp, by = "symbol") -> ctrp_gene_list

fn_filter_drug_ctrp <- function(.x, .cor = 0.25, .p_val = 0.05) {
  .x %>% dplyr::filter(abs(cor_sprm) > .cor, p_val < .p_val)
}
ctrp_gene_list %>% 
  dplyr::mutate(cor_drug = purrr::map(.x = drug, .f = fn_filter_drug_ctrp)) %>% 
  tidyr::unnest(cor_drug) -> ctrp_gene_list_sig_drug
ctrp_gene_list_sig_drug %>%
  write.table(file = file.path(out_path,"d_1_drug/ctrp_gene_list_sig_drug.txt"),quote = F)

ctrp_gene_list_sig_drug %>% 
  dplyr::mutate(
    p_val = ifelse(-log10(p_val) > 50, 50, -log10(p_val)),
    # cor_sprm = ifelse(cor_sprm > 0.5, 0.5, cor_sprm),
    # cor_sprm = ifelse(cor_sprm < -0.5, -0.5, cor_sprm)
  ) -> ctrp_plot_ready

ctrp_plot_ready %>% 
  dplyr::group_by(symbol) %>% 
  dplyr::summarise(cor_sum = sum(cor_sprm)) %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  #dplyr::mutate(color = ifelse(status == "l", "black", "red")) %>% 
  dplyr::arrange(cor_sum) -> ctrp_gene_rank
ctrp_gene_rank$col <- ctrp_gene_rank$col %>% as.character()
ctrp_gene_rank$size <- ctrp_gene_rank$size %>% as.character()

ctrp_plot_ready %>% 
  ggplot(aes(x = symbol, y = drug_name, color = cor_sprm)) +
  geom_point(aes(size = p_val)) +
  scale_x_discrete(limits = ctrp_gene_rank$symbol, expand = c(0.012,0.012)) +
  scale_y_discrete(
    # limits = drug_rank$drug_name, 
    expand = c(0.012,0.012), 
    position = "right") +
  scale_color_gradient2(
    name = "Spearman Correlation",
    high = "red",
    mid = "white",
    low = "blue"
  ) +
  scale_size_continuous(
    name = "P-value"
  ) +
  # ggthemes::theme_gdocs() +
  theme(
    panel.background = element_rect(color = "black", fill = "white", size = 0.1),
    panel.grid=element_line(colour="grey",linetype="dashed"),
    panel.grid.major=element_line(colour="grey",linetype="dashed",size=0.2),
    
    axis.title = element_blank(),
    axis.text.x = element_text(
      size = 9, 
      angle = 90,
      hjust = 1, 
      vjust = 0.5, 
      color = ctrp_gene_rank$col,face=ctrp_gene_rank$size),
    axis.text.y = element_text(
      # color = drug_rank$color, 
      size = 10),
    
    axis.ticks = element_line(color = "black"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    # legend.key.width = unit(1,"cm"),
    # legend.key.heigh = unit(0.3,"cm"),
    legend.key = element_rect(fill="white",colour = "black")
  ) + guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = 0.5,
      barwidth = 10
    )
  ) ->p

ggsave(filename = 'ctrp_sig_drug.pdf', plot = p, device = "pdf", path = ccle_path, width = 10, height = 20)

save.image(file = file.path(ccle_path, ".rda_02_ccle_a_gene_list.rda"))
load(file = file.path(ccle_path, ".rda_02_ccle_a_gene_list.rda"))

