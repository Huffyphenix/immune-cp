
# pie plot for gene pairs correlation statistic among all cancers ---------

# load data 
data_path <- "/home/huff/project/immune_checkpoint/result_20171025/e_2_DE"

cor_res <- readr::read_tsv(file.path(data_path,"tsv_06_receptor-ligand.correlation.tsv"))



# draw plot ---------------------------------------------------------------
cor_res %>%
  dplyr::mutate(cor_class = ifelse(Cor >0 & p.value < 0.05, "Pos", "Not_significant")) %>%
  dplyr::mutate(cor_class = ifelse(Cor <0 & p.value < 0.05, "Neg", cor_class)) %>%
  tidyr::nest(-symbol,-symbol.x,-cor_class,-Recepter_pairs) %>%
  dplyr::mutate(count = purrr::map(data,.f=function(.x){
    nrow(.x)
  })) %>%
  dplyr::select(symbol,symbol.x,cor_class,count,Recepter_pairs) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() %>%
  tidyr::spread(key="cor_class",value="count") %>%
  dplyr::mutate(Neg=ifelse(is.na(Neg),0,Neg)) %>%
  dplyr::mutate(Not_significant=ifelse(is.na(Not_significant),0,Not_significant)) %>%
  dplyr::mutate(Pos=ifelse(is.na(Pos),0,Pos)) %>%
  dplyr::mutate(N_all = Neg + Not_significant + Pos) %>%
  dplyr::mutate(Neg=Neg/N_all, Not_significant=Not_significant/N_all, Pos=Pos/N_all) %>%
  dplyr::select(-N_all) %>%
  tidyr::gather(key="type", value="per", -c(symbol, symbol.x, Recepter_pairs )) -> cor_res.percentage
  
cor_res.percentage %>%
  dplyr::arrange(Recepter_pairs) -> symbol_rank
cor_res.percentage %>%
  dplyr::mutate(
    symbol = factor(x = symbol, levels = unique(symbol_rank$symbol)),
    symbol.x = factor(x = symbol.x, levels = unique(symbol_rank$symbol.x))
  ) -> pie_plot_ready

pie_plot_ready %>%
  ggplot(aes(x = factor(1), y = per, fill = type)) +
  geom_bar(stat = "identity", position = "stack", color = NA) +
  # scale_y_continuous(limits = c(0,1))
  coord_polar("y") +
  facet_grid(as.formula(symbol ~ symbol.x)) + # cancer_types ~ symbol
  # scale_x_discrete(limits = cnv_gene_rank$symbol) +
  # scale_x_discrete(expand=c(0,0)) +
  # scale_y_discrete(expand=c(0,0)) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    
    strip.text.y = element_text(angle = 0, hjust = 0, size = 10),
    strip.text.x = element_text(size = 10, angle = 0, vjust = 0),
    strip.background = element_blank(),
    
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    legend.key.size = unit(0.25, "cm"),
    
    panel.background = element_blank(),
    panel.spacing = unit(0, "null"), # unit(0.01, "lines"),
    panel.spacing.x = unit(0, "null"),
    
    plot.margin = rep(unit(0, "null"), 4),
    axis.ticks.length = unit(0, "cm")
  ) +
  scale_fill_manual(
    limits = c("Neg", "Not_significant", "Pos"),
    label = c("Negative correlation", "Not significant", "Positive correlation"),
    # Amp RColorBrewer name = "Spectral"
    # Del RColorBrewer name = "BrBG"
    values = c("aquamarine3", "grey", "brown1")
  ) -> p;p

ggsave(file.path(data_path,"fig_09_receptor-ligand.correlation-pie.pdf"),device = "pdf",height = 10, width = 10)
