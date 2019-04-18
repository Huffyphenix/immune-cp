############################################
# get association between all ICPs mutation and TIL.
############################################

library(magrittr)
library(ggplot2)
out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")


# LOAD DATA 
ICP_highGene_snv <- readr::read_rds(file.path(snv_path,"ICP_highGene_snv.rds.gz")) %>%
  dplyr::select(-highGene_SNV)

# data from TIMER
immunity_path_2 <- "/home/huff/project/immune_checkpoint/data/immunity"
TIMER_immunity_onlyTumor <- readr::read_tsv(file.path(immunity_path_2,"immuneEstimation.txt")) %>%
  dplyr::filter(substr(barcode,14,14) == 0) %>%
  dplyr::mutate(barcode = substr(barcode,1,12)) %>%
  dplyr::mutate(TIL = B_cell+CD4_Tcell+CD8_Tcell+Neutrophil+Macrophage+Dendritic) 


# GET TIL DIFFERENCE BETWEEN GROUPS DEFINED BY ICPS MUTATION OR NOT -----------------------

library(showtext)
font_add("Arial","ARIAL.TTF") # sudo进入container，将/home/huff/fonts/中的字体拷贝到/usr/share/fonts/,then do:fc-cache

library(export)
source("/home/huff/project/github/immune-cp/autophagy_codes/03_mutation/fn_get_mutation_exclusive.R")
library(survival)

fn_mutation_burden_all <- function(data,group,filter,value,facet,color,xlab,title,comp_list,m_a_name,result_path,w=4,h=3){
  data %>%
    dplyr::filter(Cell_type %in% filter$Cell_type) -> data
  
  color %>%
    dplyr::inner_join(data,by = "group") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = paste(group,", n=",n,sep="")) %>%
    dplyr::select(group,color) %>%
    unique() -> color_paired
  
  data %>%
    dplyr::arrange(group) %>%
    ggpubr::ggboxplot(x = group, y = value,fill = "white",alpha = 0,width = 0.1,
                      color = group #add = "jitter",#, palette = "npg"
    ) +
    geom_violin(aes(fill = group),alpha = 0.5) +
    # geom_jitter(aes(color = group),alpha = 0.2,width = 0.1,size=0.1) +
    geom_boxplot(fill = "white",alpha = 0,width = 0.1) +
    scale_x_discrete(#breaks = c(1:10),
      limits = sort(unique(data$group))
      # labels = unique(data$group)
      # expand = c(0.2,0.2,0.2)
    ) +
    facet_wrap(as.formula(facet), strip.position = "bottom", scales = "free") +
    scale_color_manual(
      values = color_paired$color
    ) +
    scale_fill_manual(
      values = color_paired$color
    ) +
    # ylim(4,12) +
    labs(y = xlab, title = title) +
    theme(legend.position = "none",
          axis.text = element_text(colour = "black"),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill = "white",colour = "white"),
          text = element_text(size = 12, colour = "black"),
          strip.text = element_text(size = 12)) +
    ggpubr::stat_compare_means()
    # ggpubr::stat_compare_means(comparisons = list(c("ICPs_Mut", "WT")),method = "wilcox.test",label = "p.signif")
    
  ggsave(filename = paste(m_a_name,"png",sep = "."), path = result_path,device = "png",width = w,height = h)
  ggsave(filename = paste(m_a_name,"pdf",sep = "."), path = result_path,device = "pdf",width = w,height = h)
  print("end draw --------------")
}


fn_diff_test <- function(data){
  data$group %>% table() -> x
  which(x > 3) -> index
  names(index) -> valid_group
  data %>%
    dplyr::filter(group %in% valid_group) -> data
  if (length(valid_group) == 2) {
    broom::tidy(
      wilcox.test(value ~ group, data = data, var.equal = F)
    )
  } else {
    tibble::tibble()
  }
  
}

fn_TIL_on <- function(V1, V2, .data1, .data2, cancer_types){
  print(paste(cancer_types,sep = "_"))
  
  .data1 %>% 
    dplyr::mutate(ICPs = ifelse(ICPs == 0, "WT", "ICPs_Mut")) %>%
    dplyr::inner_join(TIMER_immunity_onlyTumor, by = "barcode") %>%
    tidyr::gather(-barcode, -ICPs, key = "Cell_type",value = "value") %>%
    dplyr::rename("group" = "ICPs") -> plot_ready
  
  # get TIL diff between groups, get pvalue 
  plot_ready %>%
    tidyr::nest(-Cell_type) %>%
    dplyr::mutate(pvalue = purrr::map(data,.f = fn_diff_test)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() -> result
  
  # plot if pvalue is significant
  result %>%
    dplyr::filter(p.value <= 0.05) -> result.sig
  if (nrow(result.sig) >= 1) {
    color_list <- tibble::tibble(group = c( "WT", "ICPs_Mut"), 
                                 color = c("#00B2EE", "#CD2626"))
    if (nrow(result.sig) <= 3 ) {
      w = 4
      h = 3} else if (nrow(result.sig) > 3) {
        w = 6
        h = 5
      } else {
        w = 6
        h = 7
      } 
    fn_mutation_burden_all(data = plot_ready,group = "group",filter = result.sig,title = cancer_types,
                           value = "value", facet = "~ Cell_type",color = color_list, 
                           xlab = "TIL",m_a_name = cancer_types,
                           result_path = file.path(snv_path,"TIL"),
                           w = w,h = h)
  }
  
  result %>%
    dplyr::mutate(cancer_types = cancer_types)
}

fn_TIL_pre <- function(cancer_types, ICP_SNV){
  ICP_SNV %>%
    tidyr::gather(-symbol,key = "barcode",value = "mut") %>%
    dplyr::mutate(mut = ifelse(is.na(mut),0,mut)) %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(mut_overall = sum(mut)) %>%
    dplyr::mutate(mut_overall = ifelse(mut_overall > 0,1,mut_overall)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-symbol,-mut) %>%
    dplyr::rename("ICPs" = "mut_overall") %>%
    unique() -> ICP_overall_mut
  
  fn_TIL_on(.data1 = ICP_overall_mut, cancer_types = cancer_types) -> .pval  

  .pval 
}

# RUNNING
ICP_highGene_snv %>% 
  # tail(20) %>%
  purrr::pmap(.f = fn_TIL_pre) %>% 
  dplyr::bind_rows()  -> ICPs_muation_groups_TIL

ICPs_muation_groups_TIL %>%
  dplyr::select(cancer_types, 1:5) %>%
  readr::write_tsv(file.path(snv_path, "TIL", "allICPs_mutation_groups-TIL.tsv"))

ICPs_muation_groups_TIL %>%
  dplyr::select(cancer_types, 1:5) %>%
  dplyr::filter(p.value <= 0.05) %>%
  readr::write_tsv(file.path(snv_path, "TIL", "allICPs_mutation_groups-TIL-sigP.tsv"))
