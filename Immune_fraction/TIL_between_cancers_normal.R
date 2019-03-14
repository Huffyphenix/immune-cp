####################
library(magrittr)
library(ggplot2)
library(ggpubr)
# immune infitration of each cancers --------------------------------------


# data path ---------------------------------------------------------------
# 1. data from : The Immune Landscape of Cancer. Immunity 48, 812-830.e14.
immunity_path_1 <- "/project/huff/huff/data/TCGA/immune_infiltration/immune_lanscape.immunity"
immune_lanscape_immunity <- readr::read_tsv(file.path(immunity_path_1,"mmc2.txt")) 

# 2. data from TIMER
immunity_path_2 <- "/project/huff/huff/immune_checkpoint/data/immunity"
TIMER_immunity <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic)


# result path -------------------------------------------------------------
result_path <- "/project/huff/huff/immune_checkpoint/result_20171025/e_5_immune_infiltration"

# load data ---------------------------------------------------------------
TCGA_cancer_info <- readr::read_tsv(file.path("/project/huff/huff/data/TCGA", "TCGA_sample_cancer_info.txt")) %>%
  dplyr::rename("barcode" = "TCGA Participant Barcode","cancer_types" = "TCGA Study")


# density plot ------------------------------------------------------------
cancer_color <- readr::read_tsv(file.path("/data/shiny-data/GSCALite","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(immune_lanscape_immunity$`TCGA Study`)) -> cancer21_color

### All cancer samples ----
# 1. for data from : The Immune Landscape of Cancer. Immunity 48, 812-830.e14.
fn_density_peak <-function(.name,.x){
  print(.name)
  # get peak and secondary peak
  set.seed(42)
  .x <- .x$Leukocyte_Fraction
  if(length(unique(as.vector(.x)))<3){
    .x <- .x + runif(length(.x),0.001,0.002)
  }
  d1 <- density(.x)
  d1.ym <- which.max(d1$y)
  d1.ym.x <- d1$x[d1.ym]
  d1.ym.y <- d1$y[d1.ym]
  
  tibble::tibble(peak.x = d1.ym.x, peak.y = d1.ym.y)
}

immune_lanscape_immunity %>%
  dplyr::filter(!is.na(`Leukocyte Fraction`)) %>%
  dplyr::rename("Leukocyte_Fraction"="Leukocyte Fraction","cancer_types" = "TCGA Study") %>%
  dplyr::select(`TCGA Participant Barcode`,cancer_types,Leukocyte_Fraction) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(LF_peak = purrr::map2(cancer_types,data,fn_density_peak)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> immune_lanscape_immunity.peak

immune_lanscape_immunity %>%
  dplyr::filter(!is.na(`Leukocyte Fraction`)) %>%
  dplyr::rename("Leukocyte_Fraction"="Leukocyte Fraction","cancer_types" = "TCGA Study") %>%
  ggdensity(x="Leukocyte_Fraction",fill="cancer_types",
            palette = "jco") +
  facet_wrap(~cancer_types) +
  geom_vline(data=immune_lanscape_immunity.peak,aes(xintercept = peak.x)) +
  geom_text(data=immune_lanscape_immunity.peak,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
  scale_fill_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )
  
ggsave(file.path(result_path,"Immunity.48.812-830.e14.Leukocyte_Fraction_allsamples_each_cancers.pdf"),device = "pdf",width = 10,height = 8)
ggsave(file.path(result_path,"Immunity.48.812-830.e14.Leukocyte_Fraction_allsamples_each_cancers.png"),device = "png",width = 10,height = 8)


# 2. for data from TIMER, Sherly Liu.
# function get immune cells' density peak of tumor and normal samples of each cancers, 
fn_density_peak_timer <-function(.x){
  # get peak and secondary peak
  set.seed(42)
  .x <- .x$Immune_Score
  if(length(.x)<10){
    tibble::tibble(peak.x = NA, peak.y = NA)
  }else{
    d1 <- density(.x)
    d1.ym <- which.max(d1$y)
    d1.ym.x <- d1$x[d1.ym]
    d1.ym.y <- d1$y[d1.ym]
    tibble::tibble(peak.x = d1.ym.x, peak.y = d1.ym.y)
  }
}

# data prepare
TIMER_immunity %>%
  dplyr::rename("barcode_1"="barcode") %>%
  dplyr::mutate(barcode = substr(barcode_1,1,12)) %>%
  dplyr::inner_join(TCGA_cancer_info,by="barcode") %>%
  dplyr::mutate(Type = ifelse(substr(barcode_1,14,14)=="1","Normal","Tumor")) -> TIMER_data_dealed

# get peak
TIMER_data_dealed %>%
  tidyr::gather(-cancer_types,-Type,-barcode_1,-barcode,key="cell_type",value="Immune_Score") %>%
  tidyr::nest(-cancer_types,-Type,-cell_type) %>%
  dplyr::mutate(peak = purrr::map(data,fn_density_peak_timer)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> TIMER_data_dealed.peak

# function to draw density plot 
fn_draw_density <- function(value){
  TIMER_data_dealed.peak %>%
    dplyr::filter(cell_type==value) -> anno
  
  cancer_color %>%
    dplyr::filter(cancer_types %in% unique(TIMER_data_dealed$cancer_types)) -> cancer_color.timer
  
  TIMER_data_dealed %>%
    ggdensity(x=value,fill="Type",
              palette = "jco") +
    facet_wrap(~cancer_types,scales = "free") +
    geom_vline(data=anno,aes(xintercept = peak.x,color=Type),linetype="dotted") +
    scale_color_manual(
      values = c("#104E8B", "#8B6914")
    )+
    # geom_text(data=anno,
    #           aes(x=peak.x+0.1,y=peak.y+0.2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep="")),
    #           size=5) +
    theme(
      # legend.position = "none",
      legend.key = element_rect(fill = "white", colour = "black"),
      axis.text = element_text(colour = "black",size = 12),
      strip.background = element_rect(fill="white",color="black"),
      strip.text = element_text(color="black",size=12)
    )
  filename <- paste("TIMER",value,"density.normal_tumor_compare",sep=".")
  ggsave(file.path(result_path,paste(filename,"png",sep=".")),device = "png",width = 10,height = 8)
  ggsave(file.path(result_path,paste(filename,"pdf",sep=".")),device = "pdf",width = 10,height = 8)
}

# draw density plot and output
fn_draw_density("TIL")
fn_draw_density("B_cell")
fn_draw_density("CD4_Tcell")
fn_draw_density("CD8_Tcell")
fn_draw_density("Neutrophil")
fn_draw_density("Macrophage")
fn_draw_density("Dendritic")

# function to draw density plot (by cancer types)
fn_draw_density_cancer <- function(value){
  TIMER_data_dealed.peak %>%
    dplyr::filter(cancer_types==value) -> anno
  
  TIMER_data_dealed %>%
    dplyr::filter(cancer_types==value) %>%
    tidyr::gather(-cancer_types,-Type,-barcode_1,-barcode,key="cell_type",value = "TIL") %>%
    ggdensity(x="TIL",fill="Type",
              palette = "jco") +
    facet_wrap(~cell_type,scales = "free") +
    geom_vline(data=anno,aes(xintercept = peak.x,color=Type),linetype="dotted") +
    scale_color_manual(
      values = c("#104E8B", "#8B6914")
    )+
    # geom_text(data=anno,
    #           aes(x=peak.x+0.1,y=peak.y+0.2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep="")),
    #           size=5) +
    ylab("Density") +
    theme(
      # legend.position = "none",
      legend.key = element_rect(fill = "white", colour = "black"),
      axis.text = element_text(colour = "black",size = 12),
      strip.background = element_rect(fill="white",color="black"),
      strip.text = element_text(color="black",size=12),
      axis.title.x = element_blank()
    )
  filename <- paste("TIMER",value,"density.normal_tumor_compare",sep=".")
  ggsave(file.path(result_path,paste(filename,"png",sep=".")),device = "png",width = 6,height = 4)
  ggsave(file.path(result_path,paste(filename,"pdf",sep=".")),device = "pdf",width = 6,height = 4)
  "Done"
}

# draw density plot and output
TIMER_data_dealed.peak.compareTN %>%
  dplyr::filter(peak.pos=="peak.x") %>%
  dplyr::filter(Relationship_between_T.N!="Not Applicable") %>%
  dplyr::select(cancer_types) %>%
  unique() %>%
  dplyr::mutate(x=purrr::map(cancer_types,fn_draw_density_cancer)) %>%
  tidyr::unnest()

# compare the density peak between tumor and normal
TIMER_data_dealed.peak %>%
  tidyr::gather(-cancer_types,-Type,-cell_type,key="peak.pos",value="peak.value") %>%
  tidyr::spread(key="Type",value="peak.value") %>%
  dplyr::mutate(Relationship_between_T.N = purrr::map2(Normal,Tumor,.f=function(.x,.y){
      if(is.na(.x)){
        "Not Applicable"
      }else if(.x > .y){
        "N>T"
      }else if(.x < .y){
        "T>N"
      }else{
        "T=N"
      }
    })) %>%  
  tidyr::unnest(Relationship_between_T.N) -> TIMER_data_dealed.peak.compareTN
  
TIMER_data_dealed.peak.compareTN %>%
  readr::write_tsv(file.path(result_path,"compareTN.TIMER_cell_density.peak.tsv")


# draw density compare peak
TIMER_data_dealed.peak.compareTN %>%
  dplyr::filter(peak.pos=="peak.x") %>%
  dplyr::filter(Relationship_between_T.N!="Not Applicable") -> TIMER_data_dealed.peak.compareTN.peak.x

TIMER_data_dealed.peak.compareTN.peak.x %>%
  dplyr::mutate(n = ifelse(Relationship_between_T.N=="N>T",1,-1)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(N = sum(n)) %>%
  dplyr::select(cancer_types,N) %>%
  unique() %>%
  dplyr::arrange(N) -> cancer_rank

TIMER_data_dealed.peak.compareTN.peak.x %>%
  ggplot(aes(x=cancer_types,y=cell_type)) +
  geom_tile(aes(fill=Relationship_between_T.N, width=0.9, height=0.9),color="grey",size=0.5) +
  scale_x_discrete(
    limits = cancer_rank$cancer_types
  ) +
  scale_fill_manual(
    name = "TIL (T:N) ",
    values = c("#1E90FF", "#FF8C00")
  ) +
  theme(
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12),
    axis.title = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
  )
ggsave(file.path(result_path,paste("compareTN.TIMER_cell_density.peak","pdf",sep=".")),device = "pdf",width = 6,height = 3)
ggsave(file.path(result_path,paste("compareTN.TIMER_cell_density.peak","png",sep=".")),device = "png",width = 6,height = 3)

  