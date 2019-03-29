##########################
# use xCell data from cancers with normal-tumor paired
# to confirm the conclusion that the TIL difference between tumor and normal influence the ICP expression between them.

####################
library(magrittr)
library(ggplot2)
library(ggpubr)
# immune infitration of each cancers --------------------------------------
basic_path <- "/home/huff/project"
xcell_path <- file.path(basic_path,"data/TCGA/immune_infiltration/xCell")


# load data ---------------------------------------------------------------
xCell_immunity <- readr::read_rds(file.path(xcell_path,"Pan14_xCell_results.rds.gz"))

xCell_cell_type <- readr::read_tsv(file.path(xcell_path,"cell_type_classification.13059_2017_1349_MOESM1_ESM.txt")) %>%
  dplyr::select(`Full name`,`Cell types`,Subgroup) %>%
  dplyr::mutate(cell_set = ifelse(Subgroup == "Lymphoid" | Subgroup == "Myeloid","Immune Cell","Stromal Cell")) %>%  # fantom5 and xCell overlapped cells
  dplyr::select(`Cell types`,cell_set) %>%
  unique()

# TCGA data
TCGA_cancer_info <- readr::read_tsv(file.path("/project/huff/huff/data/TCGA", "TCGA_sample_cancer_info.txt")) %>%
  dplyr::rename("barcode" = "TCGA Participant Barcode","cancer_types" = "TCGA Study")


# density compare ------------------------------------------------------------
cancer_color <- readr::read_tsv(file.path(basic_path,"data/TCGA","02_pcc.tsv"))
# cancer_color %>%
#   dplyr::filter(cancer_types %in% unique(immune_lanscape_immunity$`TCGA Study`)) -> cancer21_color

# result path -------------------------------------------------------------
result_path <- file.path(basic_path,"immune_checkpoint/result_20171025/e_5_immune_infiltration","xCell")

# xCell data process ------------------------------------------------------
# to get section score of immune and stromal cells 

xCell_immunity %>%
  dplyr::mutate(xCell_gather = purrr::map(xCell,.f=function(.x){
    .x %>%
      as.data.frame() %>%
      dplyr::mutate(`Cell types` = rownames(.x)) %>%
      tidyr::gather(-`Cell types`,key="barcode",value="score")
  })) %>%
  dplyr::select(-xCell) %>%
  tidyr::unnest() %>%
  dplyr::inner_join(xCell_cell_type,by="Cell types") %>%
  dplyr::group_by(barcode,cell_set) %>%
  dplyr::mutate(score_sum = sum(score)) %>%
  dplyr::select(cancer_types,barcode,cell_set,score_sum) %>%
  unique() %>%
  dplyr::ungroup() -> xCell_immunty.cell_set

# function get immune cells' density peak of tumor and normal samples of each cancers, 
fn_density_peak_timer <-function(.x){
  # get peak and secondary peak
  set.seed(42)
  .x <- .x$score_sum
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
xCell_immunty.cell_set %>%
  dplyr::mutate(Type = ifelse(substr(barcode,14,14)=="1","Normal","Tumor")) %>%
  dplyr::mutate(cell_set = ifelse(cell_set=="Immune Cell","Immune_Cell","Stromal_Cell")) %>%
  tidyr::spread(key="cell_set",value="score_sum")  -> xCell_immunty.cell_set_dealed

# get peak
xCell_immunty.cell_set_dealed %>%
  tidyr::gather(-cancer_types,-barcode,-Type,key="cell_set",value="score_sum") %>%
  tidyr::nest(-cancer_types,-Type,-cell_set) %>%
  dplyr::mutate(peak = purrr::map(data,fn_density_peak_timer)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> xCell_immunty.cell_set.peak

# function to draw density plot 
fn_draw_density <- function(value){
  xCell_immunty.cell_set.peak %>%
    dplyr::filter(cell_set==value) -> anno
  
    cancer_color %>%
    dplyr::filter(cancer_types %in% unique(xCell_immunty.cell_set_dealed$cancer_types)) -> cancer_color.timer
  
  xCell_immunty.cell_set_dealed %>%
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
  filename <- paste("xCell",value,"density.normal_tumor_compare",sep=".")
  ggsave(file.path(result_path,paste(filename,"png",sep=".")),device = "png",width = 10,height = 8)
  ggsave(file.path(result_path,paste(filename,"pdf",sep=".")),device = "pdf",width = 10,height = 8)
}

fn_draw_density("Immune_Cell")
fn_draw_density("Stromal_Cell")

### tumor-normal paired samples ------------------------------------------------------

# 1. for data from TIMER, Sherly Liu.
# function get immune cells' density peak of tumor and normal samples of each cancers, 
# fn_density_peak_timer 

# data prepare
xCell_immunty.cell_set_dealed %>%
  dplyr::filter(Type=="Normal") %>%
  dplyr::select(barcode) %>%
  dplyr::mutate(barcode_1 = substr(barcode,1,12)) %>%
  unique() -> xCell_T.N.paired_sample

# get peak
xCell_immunty.cell_set_dealed %>%
  dplyr::mutate(barcode_1 = substr(barcode,1,12)) %>%
  dplyr::filter(barcode_1 %in% xCell_T.N.paired_sample$barcode_1) %>%
  tidyr::gather(-cancer_types,-Type,-barcode,-barcode_1,key="cell_set",value="score_sum") %>%
  tidyr::nest(-cancer_types,-Type,-cell_set) %>%
  dplyr::mutate(peak = purrr::map(data,fn_density_peak_timer)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> xCell_data_T.N.paired.peak

# function to draw density plot 
fn_draw_density_TN_paired <- function(data,peak_data,value){
  data %>%
    .$cancer_types %>% unique() -> cancers
  peak_data %>%
    dplyr::filter(cell_set==value) %>%
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
  filename <- paste("xCell",value,"density.TN-paired_compare",sep=".")
  ggsave(file.path(result_path,paste(filename,"png",sep=".")),device = "png",width = 10,height = 8)
  ggsave(file.path(result_path,paste(filename,"pdf",sep=".")),device = "pdf",width = 10,height = 8)
}

# draw density plot and output
xCell_immunty.cell_set_dealed %>%
  dplyr::mutate(barcode_1 = substr(barcode,1,12)) %>%
  dplyr::filter(barcode_1 %in% xCell_T.N.paired_sample$barcode_1) %>%
  dplyr::group_by(cancer_types,Type) %>%
  dplyr::mutate(n=n()) %>%
  dplyr::filter(n>=10) %>%
  dplyr::select(-n) %>%
  dplyr::ungroup() -> xCell_ready_draw
fn_draw_density_TN_paired(xCell_ready_draw,xCell_data_T.N.paired.peak,"Immune_Cell")
fn_draw_density_TN_paired(xCell_ready_draw,xCell_data_T.N.paired.peak,"Stromal_Cell")


# function to draw density plot (by cancer types)
fn_draw_density_cancer_paired <- function(value){
  xCell_ready_draw %>%
    .$cancer_types %>% unique() -> cancers
  xCell_data_T.N.paired.peak %>%
    dplyr::filter(cancer_types==value) %>%
    dplyr::filter(cancer_types %in% cancers) -> anno
  
  xCell_ready_draw %>%
    dplyr::filter(cancer_types==value) %>%
    tidyr::gather(-cancer_types,-Type,-barcode_1,-barcode,key="cell_set",value = "score_sum") %>%
    ggdensity(x="score_sum",fill="Type",
              palette = "jco") +
    facet_wrap(~cell_set,scales = "free") +
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
  filename <- paste("xCell",value,"density.T-N.paired_compare",sep=".")
  ggsave(file.path(result_path,paste(filename,"png",sep=".")),device = "png",width = 6,height = 5)
  ggsave(file.path(result_path,paste(filename,"pdf",sep=".")),device = "pdf",width = 6,height = 5)
  "Done"
}
# compare the density peak between tumor and normal
xCell_data_T.N.paired.peak %>%
  tidyr::gather(-cancer_types,-Type,-cell_set,key="peak.pos",value="peak.value") %>%
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
  tidyr::unnest(Relationship_between_T.N) -> xCell_data.peak.compareTN.paired

# draw density plot and output
xCell_data.peak.compareTN.paired %>%
  dplyr::filter(peak.pos=="peak.x") %>%
  dplyr::filter(Relationship_between_T.N!="Not Applicable") %>%
  dplyr::select(cancer_types) %>%
  unique() %>%
  dplyr::mutate(x=purrr::map(cancer_types,fn_draw_density_cancer_paired)) %>%
  tidyr::unnest()

xCell_data.peak.compareTN.paired %>%
  readr::write_tsv(file.path(result_path,"compareTN.paired.xCell_cell_density.peak.tsv"))


# draw density compare peak
xCell_data.peak.compareTN.paired %>%
  dplyr::filter(peak.pos=="peak.x") %>%
  dplyr::filter(Relationship_between_T.N!="Not Applicable") -> xCell_data.peak.compareTN.paired.peak.x

xCell_data.peak.compareTN.paired.peak.x %>%
  dplyr::mutate(n = ifelse(cell_set=="Immune_Cell" & Relationship_between_T.N=="N>T",1,-1)) %>%
  dplyr::mutate(n = ifelse(cell_set=="Stromal_Cell" & Relationship_between_T.N=="T>N",1,n)) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(N = sum(n)) %>%
  dplyr::select(cancer_types,N) %>%
  unique() %>%
  dplyr::arrange(N) -> cancer_rank.paired

xCell_data.peak.compareTN.paired.peak.x %>%
  ggplot(aes(x=cancer_types,y=cell_set)) +
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
ggsave(file.path(result_path,paste("compareTN.paired.xCell_cell_density.peak","pdf",sep=".")),device = "pdf",width = 6,height = 2)
ggsave(file.path(result_path,paste("compareTN.paired.xCell_cell_density.peak","png",sep=".")),device = "png",width = 6,height = 2)
