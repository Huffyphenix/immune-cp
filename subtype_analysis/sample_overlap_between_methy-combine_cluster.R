library(methods)
library(magrittr)
library(CancerSubtypes)
library(ConsensusClusterPlus)
library(SNFtool)
library(survival)
library(survminer)
library(ggplot2)

# load methy cluster  -----------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

methy_results <- readr::read_rds(file.path("/project/huff/huff/immune_checkpoint/genelist_data","genelist_methy_CC_20.rds.gz"))

C=6
methy_cluster <- methy_results[[C]][['consensusClass']] 
methy_group <- data.frame(barcode = names(methy_cluster),Cluster=methy_cluster) %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(M_Group = ifelse(Cluster %in% c(1,2,3), "Group1","Group2")) %>% # for 6 clusters
  dplyr::mutate(M_Group = ifelse(Cluster == 6, "Group3", M_Group))

combine_cluster <- readr::read_tsv(file.path("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/combine/Get_best_clutser_20",paste(C,"clusters_group_info.tsv",sep="_")))
combine_cluster %>%
  dplyr::rename("barcode" = "sample","Cluster" = "group") %>%
  dplyr::mutate(C_Group = ifelse(Cluster %in% c(1,3,5,6), "Group1","Group2")) %>% # for 6 clusters
  dplyr::mutate(C_Group = ifelse(Cluster == 2, "Group3", C_Group)) -> combine_group

combine_group %>%
  dplyr::inner_join(methy_group,by="barcode") %>%
  dplyr::select(barcode,C_Group,M_Group) %>%
  dplyr::mutate(same_group = ifelse(M_Group==C_Group,"overlap","diff")) -> methy_combine_group

methy_combine_group %>%
  dplyr::select(same_group) %>%
  table()

methy_combine_group %>%
  dplyr::group_by(C_Group, M_Group) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::select(C_Group, M_Group,n) %>%
  unique()

methy_combine_group %>%
  dplyr::select(-same_group) %>%
  dplyr::mutate(line_group = paste(C_Group,M_Group,sep=".")) %>%
  dplyr::arrange(C_Group,line_group) %>%
  dplyr::mutate(C_n=1:nrow(methy_combine_group)) %>%
  dplyr::arrange(M_Group,line_group) %>%
  dplyr::mutate(M_n=1:nrow(methy_combine_group)) %>%
  dplyr::arrange(C_Group,line_group) %>%
  dplyr::mutate(ID=1:nrow(methy_combine_group)) -> methy_combine_group.draw
methy_combine_group.draw %>%
  dplyr::select(barcode,C_Group,C_n,ID,line_group) %>%
  dplyr::mutate(`Data Platform`= "C_Group") %>%
  dplyr::rename("Group"="C_Group","n" = "C_n") -> methy_combine_group.combine
methy_combine_group.draw %>%
  dplyr::select(barcode,M_Group,M_n,ID,line_group) %>%
  dplyr::mutate(`Data Platform`= "M_Group") %>%
  dplyr::rename("Group"="M_Group","n" = "M_n") -> methy_combine_group.methy
methy_combine_group.combine %>%
  rbind(methy_combine_group.methy) -> plot_ready
# plot 参考
# https://stackoverflow.com/questions/44351127/how-to-plot-parallel-coordinates-with-multiple-categorical-variables-in-r
plot_ready %>%
  ggplot(aes(x=`Data Platform`,y=n,group=ID)) +
  geom_line(aes(color=line_group),alpha=0.05) +
  geom_point(aes(color = Group),shape=0) +
  scale_color_manual(
    # breaks = c(plot_ready$Group %>% unique(),plot_ready$line_group %>% unique()),
    values =  c("#CD2626", "#CDC0B0", "#838B8B", "#000000", "#CD9B1D","#0000FF", "#00008B", "#8A2BE2",  "#00B2EE","#98F5FF", "#53868B", "#EEAD0E", "#458B00", "#EEA2AD", "#E066FF", "#EE3A8C", "#00FF00", "#FFFF00", "#5CACEE", "#8B6914", "#FF7F24")
  ) +
  theme(
    axis.text.y = element_blank(),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.background = element_rect(fill = "white"),
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12,color = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_blank()
  )
result_path <- file.path("/project/huff/huff/immune_checkpoint/result_20171025/subtype_analysis/combine")
ggsave(file.path(result_path,"Overlap_between.Methy.3Combined.result.png"),device = "png",width = 4,height = 4)
ggsave(file.path(result_path,"Overlap_between.Methy.3Combined.result.pdf"),device = "pdf",width = 4,height = 4)

