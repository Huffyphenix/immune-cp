### TIMER data
### TIL distribution

TIMER_data <- readr::read_tsv(file.path("/project/huff/huff/immune_checkpoint/data/immunity","immuneEstimation.txt"))

sample_cacner <- readr::read_tsv(file.path("/home/huff/project/data/TCGA","TCGA_RNAseq_sample_info.tsv")) 

TIMER_data %>%
  dplyr::inner_join(sample_cacner,by="barcode") %>%
  dplyr::filter(substr(barcode,14,14)=="0") %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+ CD8_Tcell+ Neutrophil+ Macrophage+ Dendritic) -> TIMER_data_cancer

# get data peak =----------
fn_density_peak <-function(.name,.x){
  print(.name)
  # get peak and secondary peak
  set.seed(42)
  .x <- .x$Leukocyte_Fraction
  if(length(unique(as.vector(.x)))<3){
    tibble::tibble(peak.x = NA, peak.y = NA)
  }else{
    d1 <- density(.x)
    d1.ym <- which.max(d1$y)
    d1.ym.x <- d1$x[d1.ym]
    d1.ym.y <- d1$y[d1.ym]
    
    tibble::tibble(peak.x = d1.ym.x, peak.y = d1.ym.y)
  }
}

TIMER_data_cancer %>%
  tidyr::gather(key="cell_type",value="Leukocyte_Fraction",-barcode,-cancer_types) %>%
  tidyr::nest(-cancer_types,-cell_type) %>%
  dplyr::mutate(LF_peak = purrr::map2(cancer_types,data,fn_density_peak)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> timer.peak

TIMER_data_cancer %>%
  tidyr::gather(key="cell_type",value="Leukocyte_Fraction",-barcode,-cancer_types) %>%
  dplyr::inner_join(timer.peak,by=c("cancer_types","cell_type")) -> immune_lanscape_immunity.plot

res_path <- "/home/huff/project/immune_checkpoint/result_20171025/e_5_immune_infiltration/TIMER_cancer_TIL_density"

cancer_color <- readr::read_tsv(file.path("/home/huff/project","data/TCGA","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(TIMER_data_cancer$cancer_types)) -> cancer21_color

for (celltype in c("B_cell", "CD4_Tcell", "CD8_Tcell", "Neutrophil", "Macrophage", "Dendritic","TIL")){
  immune_lanscape_immunity.plot %>%
    dplyr::filter(cell_type %in% celltype) -> data_fr_plot
  timer.peak %>%
    dplyr::filter(cell_type %in% celltype) -> peak_fr_plot
  data_fr_plot %>%
    dplyr::arrange(peak.x) %>%
    .$cancer_types %>%
    unique() -> cancer.rank
  
  data_fr_plot <- within(data_fr_plot,cancer_types <- factor(cancer_types,levels = cancer.rank))
  with(data_fr_plot, levels(cancer_types))
  
  peak_fr_plot <- within(peak_fr_plot,cancer_types <- factor(cancer_types,levels = cancer.rank))
  with(peak_fr_plot, levels(cancer_types))
  
  # plot --------------
  
  data_fr_plot %>%
    ggdensity(x="Leukocyte_Fraction",fill="cancer_types",
              palette = "jco") +
    facet_wrap(~cancer_types,scales = "free_y") +
    geom_vline(data=peak_fr_plot,aes(xintercept = peak.x),linetype = "dotted") +
    geom_text(data=peak_fr_plot,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
    scale_fill_manual(
      values = cancer21_color$color,
      limits = cancer21_color$cancer_types
    ) +
    ggtitle(paste("Density of", celltype, "in each cancers (all tumor samples)")) +
    ylab("Density") +
    theme(
      legend.position = "none",
      legend.key = element_rect(fill = "white", colour = "black"),
      axis.text = element_text(colour = "black",size = 12),
      strip.background = element_rect(fill="white",color="black"),
      strip.text = element_text(color="black",size=12)
    )
  ggsave(file.path(res_path,paste(celltype,"TIMER_cancer_density",".png",sep="")),device = "png",width = 10,height = 8)
  ggsave(file.path(res_path,paste(celltype,"TIMER_cancer_density",".pdf",sep="")),device = "pdf",width = 10,height = 8)
}

