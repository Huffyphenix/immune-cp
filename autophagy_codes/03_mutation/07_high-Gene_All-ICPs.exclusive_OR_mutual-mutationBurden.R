############################################
# get mutual exclusive mutation profile between all ICPs.
############################################

library(magrittr)
library(ggplot2)
out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")


# LOAD DATA 
ICP_highGene_snv <- readr::read_rds(file.path(snv_path,"ICP_highGene_snv.rds.gz")) 

# data from TIMER
mutation_burden_class <- readr::read_rds(file.path("/home/huff/project/data/TCGA","classfication_of_26_cancers_by_mutation_burden192.rds.gz")) %>%
  tidyr::unnest() %>%
  dplyr::rename("cancer_types"="Cancer_Types")

# survival results from 07_high-Gene_All-ICPs.exclusive_OR_mutual-survival.R
survival_sig_res <- readr::read_tsv(file.path(snv_path,"mutual_exclusive_allICPs-highGene","survival","allICPs_mutual_exclusive_highGeneSNV-survival.tsv")) %>%
  dplyr::filter(kmp <= 0.05)

# GET TIL BETWEEN GROUPS DEFINED BY ICPS AND HIGHGENE MUTATION IN EACH CANCERS -----------------------
# TIL DIFFERENCE BETWEEN EACH MUTATION GROUPS
library(showtext)
font_add("Arial","ARIAL.TTF") # sudo进入container，将/home/huff/fonts/中的字体拷贝到/usr/share/fonts/,then do:fc-cache

library(export)
source("/home/huff/project/github/immune-cp/autophagy_codes/03_mutation/fn_get_mutation_exclusive.R")
library(survival)

fn_mutation_burden_all <- function(data,group,value,color,title,xlab,m_a_name,result_path,w=4,h=3){
  
  color %>%
    dplyr::inner_join(data,by = "group") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group = paste(group,", n=",n,sep="")) %>%
    dplyr::select(group,color) %>%
    unique() -> color_paired
  comp_mat <- expand.grid("0_0",sort(unique(data$group)))
  comp_list <- list()
  for (col in 2:nrow(comp_mat)) {
    comp_list[[col - 1]] <- c(as.character(comp_mat[col,1]),as.character(comp_mat[col,2]))
  }
  
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
    # facet_wrap(as.formula(facet), strip.position = "bottom", scales = "free") +
    scale_color_manual(
      values = color_paired$color,
      labels = sort(unique(data$group))
    ) +
    scale_fill_manual(
      values = color_paired$color,
      labels = sort(unique(data$group))
    ) +
    # ylim(4,12) +
    labs(y = xlab, x = "Group", title  = title) +
    theme(legend.position = "none",
          axis.title = element_text(colour = "black"),
          strip.background = element_rect(fill = "white",colour = "white"),
          text = element_text(size = 8, colour = "black")) +
    # ggpubr::stat_compare_means(label.y = 14,paired = TRUE) +
    ggpubr::stat_compare_means(comparisons = comp_list,method = "wilcox.test",label = "p.signif")
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
  if (length(valid_group) >= 2) {
    broom::tidy(
      oneway.test(value ~ group, data = data, var.equal = F)
    ) 
  } else {
    tibble::tibble()
  }
  
}

fn_MB_on <- function(V1, V2, .data1, .data2, cancer_types){
  # .data <- filter_cnv
  # V1 <- 'TP53'
  # V2 <- 'HLA-E'
  # print("start fn_survival_on --------------")
  V1 <- as.character(V1)
  V2 <- as.character(V2)
  print(paste(cancer_types,V1,V2,sep = "_"))
  
  .data2 %>%
    dplyr::filter(symbol %in% c(V1)) %>% 
    dplyr::select(-symbol) %>%
    tidyr::gather(key = barcode, value = mut) %>%
    dplyr::mutate(mut = ifelse(is.na(mut),0,1)) %>%
    dplyr::mutate(mut = as.integer(mut)) %>%
    dplyr::inner_join(.data1, by = "barcode") %>%
    tidyr::unite(group, c("mut", V2)) %>%
    dplyr::inner_join(mutation_burden_class, by = "barcode") %>%
    dplyr::mutate(value = ifelse(sm_count == 0, log2(0.5),log2(sm_count))) -> plot_ready
  
  # get TIL diff between groups, get pvalue 
  plot_ready %>%
    dplyr::mutate(n =1) %>%
    tidyr::nest(-n) %>%
    dplyr::mutate(pvalue = purrr::map(data,fn_diff_test)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() -> result
  
  # plot if pvalue is significant
  result %>%
    dplyr::filter(p.value <= 0.1) -> result.sig
  if (nrow(result.sig) >= 1) {
    color_list <- tibble::tibble(group = c( "0_0", "0_1", "1_0", "1_1"), 
                                 color = c("#00B2EE", "#CDAD00", "pink1","#CD2626"))
    w = 4
    h = 3
    fn_mutation_burden_all(data = plot_ready,group = "group",
                           value = "value",color = color_list, 
                           title = paste(cancer_types, V1, V2),
                           xlab = "log2 (mutation burden)",m_a_name = paste(cancer_types, V1, V2, sep = "_"),
                           result_path = file.path(snv_path,"mutual_exclusive_allICPs-highGene","mutation_burden"),
                           w = w,h = h)
  }
  
  
  name <- paste(c(cancer_types, V1, V2), collapse = "_")
  # print(name)
  
  # print("end fn_survival_on --------------")
  
  result %>%
    dplyr::mutate(te = name)
}

fn_MB_pre <- function(cancer_types, ICP_SNV, highGene_SNV, survival, cluster){
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
  
  expand.grid(highGene_SNV$symbol,"ICPs") -> .gene_pairs
  ct <- cancer_types
  survival_sig_res %>%
    dplyr::filter(cancer_types %in% ct) %>%
    .$g1 -> g1
  if (length(g1) >= 1) {
    .gene_pairs %>% 
      dplyr::filter(Var1 %in% g1) %>%
      dplyr::mutate(rs = purrr::map2(Var1, Var2, .f = fn_MB_on, .data1 = ICP_overall_mut, .data2 = highGene_SNV,cancer_types = cancer_types)) %>% 
      dplyr::as_tibble() %>%
      dplyr::ungroup() %>%
      dplyr::select(rs) %>% 
      tidyr::unnest() -> .gene_pairs_pval  
  } else{
    tibble::tibble() -> .gene_pairs_pval
  }
  .gene_pairs_pval 
}

# RUNNING
ICP_highGene_snv %>% 
  # tail(20) %>%
  purrr::pmap(.f = fn_MB_pre, cluster = cluster) %>% 
  dplyr::bind_rows() %>%
  tidyr::separate(col = te, into = c('cancer_types', 'g1', 'g2')) -> mutual_exclusive_MB

mutual_exclusive_MB %>%
  readr::write_tsv(file.path(snv_path,"mutual_exclusive_allICPs-highGene", "mutation_burden", "allICPs_mutual_exclusive_highGeneSNV-mutation_burden.tsv"))

