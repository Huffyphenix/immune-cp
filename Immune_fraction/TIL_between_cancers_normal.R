####################
library(magrittr)
library(ggplot2)
library(ggpubr)
# immune infitration of each cancers --------------------------------------
basic_path <- "/home/huff/project"


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


# density compare ------------------------------------------------------------
cancer_color <- readr::read_tsv(file.path("/data/shiny-data/GSCALite","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(immune_lanscape_immunity$`TCGA Study`)) -> cancer21_color

### All cancer samples ------------------------------------------------------
# 1. for data from : The Immune Landscape of Cancer. Immunity 48, 812-830.e14.
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
  geom_vline(data=immune_lanscape_immunity.peak,aes(xintercept = peak.x),linetype = "dotted") +
  geom_text(data=immune_lanscape_immunity.peak,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
  scale_fill_manual(
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  ) +
  ggtitle("Density of LF in each cancers (all tumor samples)") +
  ylab("Density") +
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
  readr::write_tsv(file.path(result_path,"compareTN.TIMER_cell_density.peak.tsv"))


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

### tumor-normal paired samples ------------------------------------------------------

# 1. for data from TIMER, Sherly Liu.
# function get immune cells' density peak of tumor and normal samples of each cancers, 
# fn_density_peak_timer 

# data prepare
TIMER_data_dealed %>%
  dplyr::filter(Type=="Normal") %>%
  dplyr::select(barcode) %>%
  unique() -> TIMER_T.N.paired_sample

# get peak
TIMER_data_dealed %>%
  dplyr::filter(barcode %in% TIMER_T.N.paired_sample$barcode) %>%
  tidyr::gather(-cancer_types,-Type,-barcode_1,-barcode,key="cell_type",value="Immune_Score") %>%
  tidyr::nest(-cancer_types,-Type,-cell_type) %>%
  dplyr::mutate(peak = purrr::map(data,fn_density_peak_timer)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> TIMER_data_T.N.paired.peak

# function to draw density plot 
fn_draw_density_TN_paired <- function(data,peak_data,value){
  data %>%
    .$cancer_types %>% unique() -> cancers
  peak_data %>%
    dplyr::filter(cell_type==value) %>%
    dplyr::filter(cancer_types %in% cancers)-> anno
  
  data %>%
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
    ylab("Density") +
    ggtitle(paste(value,"density","in each cancers (T-N paired)")) +
    theme(
      # legend.position = "none",
      legend.key = element_rect(fill = "white", colour = "black"),
      axis.text = element_text(colour = "black",size = 12),
      strip.background = element_rect(fill="white",color="black"),
      strip.text = element_text(color="black",size=12)
    )
  filename <- paste("TIMER",value,"density.TN-paired_compare",sep=".")
  ggsave(file.path(result_path,paste(filename,"png",sep=".")),device = "png",width = 10,height = 8)
  ggsave(file.path(result_path,paste(filename,"pdf",sep=".")),device = "pdf",width = 10,height = 8)
}

# draw density plot and output
TIMER_data_dealed %>%
  dplyr::filter(barcode %in% TIMER_T.N.paired_sample$barcode) %>%
  dplyr::group_by(cancer_types,Type) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::select(-n) %>%
  dplyr::ungroup() -> TIMER_ready_draw
fn_draw_density_TN_paired(TIMER_ready_draw,TIMER_data_T.N.paired.peak,"TIL")
fn_draw_density_TN_paired(TIMER_ready_draw,TIMER_data_T.N.paired.peak,"CD4_Tcell")
fn_draw_density_TN_paired(TIMER_ready_draw,TIMER_data_T.N.paired.peak,"CD8_Tcell")
fn_draw_density_TN_paired(TIMER_ready_draw,TIMER_data_T.N.paired.peak,"Neutrophil")
fn_draw_density_TN_paired(TIMER_ready_draw,TIMER_data_T.N.paired.peak,"Macrophage")
fn_draw_density_TN_paired(TIMER_ready_draw,TIMER_data_T.N.paired.peak,"Dendritic")

# function to draw density plot (by cancer types)
fn_draw_density_cancer_paired <- function(value){
  TIMER_ready_draw %>%
    .$cancer_types %>% unique() -> cancers
  TIMER_data_T.N.paired.peak %>%
    dplyr::filter(cancer_types==value) %>%
    dplyr::filter(cancer_types %in% cancers) -> anno
  
  TIMER_ready_draw %>%
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
    ggtitle(paste("Immune cells density of paired T-N samples in", value)) +
    theme(
      # legend.position = "none",
      legend.key = element_rect(fill = "white", colour = "black"),
      axis.text = element_text(colour = "black",size = 12),
      strip.background = element_rect(fill="white",color="black"),
      strip.text = element_text(color="black",size=12),
      axis.title.x = element_blank()
    )
  filename <- paste("TIMER",value,"density.T-N.paired_compare",sep=".")
  ggsave(file.path(result_path,paste(filename,"png",sep=".")),device = "png",width = 6,height = 5)
  ggsave(file.path(result_path,paste(filename,"pdf",sep=".")),device = "pdf",width = 6,height = 5)
  "Done"
}
# compare the density peak between tumor and normal
TIMER_data_T.N.paired.peak %>%
  tidyr::gather(-cancer_types,-Type,-cell_type,key="peak.pos",value="peak.value") %>%
  tidyr::spread(key="Type",value="peak.value") %>%
  dplyr::mutate(Relationship_between_T.N = purrr::map2(Normal,Tumor,.f=function(.x,.y){
    if(is.na(.x) | is.na(.y)){
      "Not Applicable"
    }else if(.x > .y){
      "N>T"
    }else if(.x < .y){
      "T>N"
    }else{
      "T=N"
    }
  })) %>%  
  tidyr::unnest(Relationship_between_T.N) -> TIMER_data_dealed.peak.compareTN.paired

# draw density plot and output
TIMER_data_dealed.peak.compareTN.paired %>%
  dplyr::filter(peak.pos=="peak.x") %>%
  dplyr::filter(Relationship_between_T.N!="Not Applicable") %>%
  dplyr::select(cancer_types) %>%
  unique() %>%
  dplyr::mutate(x=purrr::map(cancer_types,fn_draw_density_cancer_paired)) %>%
  tidyr::unnest()

TIMER_data_dealed.peak.compareTN.paired %>%
  readr::write_tsv(file.path(result_path,"compareTN.paired.TIMER_cell_density.peak.tsv"))


# draw density compare peak
TIMER_data_dealed.peak.compareTN.paired %>%
  dplyr::filter(peak.pos=="peak.x") %>%
  dplyr::filter(Relationship_between_T.N!="Not Applicable") -> TIMER_data_dealed.peak.compareTN.paired.peak.x

TIMER_data_dealed.peak.compareTN.paired.peak.x %>%
  dplyr::mutate(n = ifelse(Relationship_between_T.N=="N>T",1,-1)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(N = sum(n)) %>%
  dplyr::select(cancer_types,N) %>%
  unique() %>%
  dplyr::arrange(N) -> cancer_rank.paired

TIMER_data_dealed.peak.compareTN.paired.peak.x %>%
  ggplot(aes(x=cancer_types,y=cell_type)) +
  geom_tile(aes(fill=Relationship_between_T.N, width=0.9, height=0.9),color="grey",size=0.5) +
  scale_x_discrete(
    limits = cancer_rank.paired$cancer_types
  ) +
  scale_fill_manual(
    name = NULL,
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
ggsave(file.path(result_path,paste("compareTN.paired.TIMER_cell_density.peak","pdf",sep=".")),device = "pdf",width = 6,height = 3)
ggsave(file.path(result_path,paste("compareTN.paired.TIMER_cell_density.peak","png",sep=".")),device = "png",width = 6,height = 3)

# 1. for data from : The Immune Landscape of Cancer. Immunity 48, 812-830.e14.
fn_density_peak_paired <- function(.name,.x){
  print(.name)
  # get peak and secondary peak
  set.seed(42)
  .x <- .x$Leukocyte_Fraction
  if(length(as.vector(.x))<10){
    tibble::tibble(peak.x = NA, peak.y = NA)
  }else{
    d1 <- density(.x)
    d1.ym <- which.max(d1$y)
    d1.ym.x <- d1$x[d1.ym]
    d1.ym.y <- d1$y[d1.ym]
    
    tibble::tibble(peak.x = d1.ym.x, peak.y = d1.ym.y)
  }
}
#  -----------------
immune_lanscape_immunity %>%
  dplyr::filter(`TCGA Participant Barcode` %in% TIMER_T.N.paired_sample$barcode) %>%
  dplyr::filter(!is.na(`Leukocyte Fraction`)) %>%
  dplyr::rename("Leukocyte_Fraction"="Leukocyte Fraction","cancer_types" = "TCGA Study") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::select(-n) %>%
  dplyr::ungroup() -> immune_lanscape_immunity.withPairedNormal

immune_lanscape_immunity.withPairedNormal %>%
  dplyr::select(`TCGA Participant Barcode`,cancer_types,Leukocyte_Fraction) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(LF_peak = purrr::map2(cancer_types,data,fn_density_peak_paired)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> immune_lanscape_immunity.peak.withPairedNormal

cancer_color %>%
  dplyr::filter(cancer_types %in% unique(immune_lanscape_immunity.withPairedNormal$cancer_types)) -> cancer_color.withPairedNormal

immune_lanscape_immunity.withPairedNormal %>%
  ggdensity(x="Leukocyte_Fraction",fill="cancer_types",
            palette = "jco",alpha = 0.5) +
  facet_wrap(~cancer_types) +
  geom_vline(data=immune_lanscape_immunity.peak.withPairedNormal,aes(xintercept = peak.x),linetype="dotted") +
  geom_text(data=immune_lanscape_immunity.peak.withPairedNormal,aes(x=peak.x+0.2,y=peak.y+2,label=paste("(",signif(peak.x,2),",",signif(peak.y,2),")",sep=""))) +
  scale_fill_manual(
    values = cancer_color.withPairedNormal$color,
    limits = cancer_color.withPairedNormal$cancer_types
  ) +
  ylab("Density") +
  ggtitle(paste("Density of LF in each cancers (with paired normal samples)")) +
  theme(
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )

ggsave(file.path(result_path,"Immunity.48.812-830.e14.Leukocyte_Fraction_Tumor-samples(with paired normal samples)_each_cancers.pdf"),device = "pdf",width = 10,height = 8)
ggsave(file.path(result_path,"Immunity.48.812-830.e14.Leukocyte_Fraction_Tumor-samples(with paired normal samples)_each_cancers.png"),device = "png",width = 10,height = 8)

save.image(file.path(result_path,"Immunity_diff_in_cancers.rda"))
load(file.path(result_path,"Immunity_diff_in_cancers.rda"))