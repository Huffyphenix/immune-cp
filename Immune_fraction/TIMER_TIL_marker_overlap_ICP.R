#########
# get the overlapping between ICGs and immune cell markers
#########


# load data ---------------------------------------------------------------

data_path <- file.path("/home/huff/project/immune_checkpoint/data/immunity")

gene_marker <- readr::read_tsv(file.path(data_path,"cancer-specific_genes_markers.txt"))
  
ICP_list <- readr::read_tsv(file.path("/home/huff/project/immune_checkpoint/checkpoint/20171021_checkpoint","ICPs_all_info_class-new.tsv"))  


# TIMER gene marker processing: genes as marker in more than 7 cancers --------------------------------------------

gene_marker %>%
  tidyr::nest(- `Immune Cell`) %>%
  dplyr::mutate(marker_filter = purrr::map(data,.f=function(.x){
    .x$`Disease Name` %>% unique() %>% length() -> cancers_number
    .x$`Gene Symbol` %>%
      table() %>%
      as.data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::filter(Freq>=(cancers_number/3))%>%
      .$.
  })) %>%
  dplyr::select(-data) -> cell_specific_markers


# overlap analysis --------------------------------------------------------

library(VennDiagram)
setwd(file.path("/home/huff/project/immune_checkpoint/result_20171025/e_5_immune_infiltration/TIMERmarkers_overlap_ICGs"))
cell_specific_markers %>%
  dplyr::mutate(run = purrr::map2(`Immune Cell`,marker_filter,.f=function(.x,.y){
    vennplot <- venn.diagram(list(TIMER_markers=.y,
                                  ICGs=ICP_list$symbol),
                             # file.path(out_path,"Venn_EHZ2_CBX2_targets_overlap.pdf"),
                             filename = NULL,
                             imagetype = "svg",
                             col = RColorBrewer::brewer.pal(6,"Paired")[c(6,2)],
                             cat.col = RColorBrewer::brewer.pal(6,"Paired")[c(6,2)],
                             # cat.pos=c(0,-5),
                             scaled = TRUE,
                             ext.text = TRUE,
                             ext.line.lwd = 2,
                             ext.dist = -0.15,
                             ext.length = 0.9,
                             ext.pos = -4,
                             inverted = TRUE,
                             # rotation.degree = 45,
                             cex = 1,
                             cat.cex = 1,
                             main = paste(.x),
                             main.cex = 2)
    ggsave(vennplot,file=paste(.x,"TIMER-moreThan7Cancers_overlap_ICGs.svg",sep="."),device = "svg", width = 4, height = 4)
    1
  }))

# TIMER gene marker processing: genes as marker in more than 2 cancers --------------------------------------------

gene_marker %>%
  tidyr::nest(- `Immune Cell`) %>%
  dplyr::mutate(marker_filter = purrr::map(data,.f=function(.x){
    .x$`Disease Name` %>% unique() %>% length() -> cancers_number
    .x$`Gene Symbol` %>%
      table() %>%
      as.data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::filter(Freq>=2)%>%
      .$.
  })) %>%
  dplyr::select(-data) -> cell_specific_markers.2cancers


# overlap analysis --------------------------------------------------------

library(VennDiagram)
cell_specific_markers.2cancers %>%
  dplyr::mutate(run = purrr::map2(`Immune Cell`,marker_filter,.f=function(.x,.y){
    vennplot <- venn.diagram(list(TIMER_markers=.y,
                                  ICGs=ICP_list$symbol),
                             # file.path(out_path,"Venn_EHZ2_CBX2_targets_overlap.pdf"),
                             filename = NULL,
                             imagetype = "svg",
                             col = RColorBrewer::brewer.pal(6,"Paired")[c(6,2)],
                             cat.col = RColorBrewer::brewer.pal(6,"Paired")[c(6,2)],
                             # cat.pos=c(0,-5),
                             scaled = TRUE,
                             ext.text = TRUE,
                             ext.line.lwd = 2,
                             ext.dist = -0.15,
                             ext.length = 0.9,
                             ext.pos = -4,
                             inverted = TRUE,
                             # rotation.degree = 45,
                             cex = 1,
                             cat.cex = 1,
                             main = paste(.x),
                             main.cex = 2)
    ggsave(vennplot,file=paste(.x,"TIMER-moreThan2Cancers_overlap_ICGs.svg",sep="."),device = "svg", width = 4, height = 4)
    1
  }))
