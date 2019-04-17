############################################
# get mutual exclusive mutation profile between all ICPs.
############################################

library(magrittr)
library(ggplot2)
out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")


# LOAD DATA 
ICP_highGene_snv <- readr::read_rds(file.path(snv_path,"ICP_highGene_snv.rds.gz")) 

# gene_list_snv.count_per.filter_hypermutation <-
#   readr::read_rds(path = file.path(snv_path, "syn7824274_snv_gene_list_snv_count_per.filter_hypermutation.rds.gz"))
# 
# ICP_highGene_snv %>%
#   dplyr::inner_join(gene_list_snv.count_per.filter_hypermutation,by="cancer_types") %>%
#   dplyr::rename("sm_count"="res") -> ICP_SNV_combine_data

# GET EXCLUSIVE MUTATION PROFILE BETWEEN ICPS IN EACH CANCERS -----------------------
library(showtext)
font_add("Arial","ARIAL.TTF") # sudo进入container，将/home/huff/fonts/中的字体拷贝到/usr/share/fonts/,then do:fc-cache

library(export)
source("/home/huff/project/github/immune-cp/autophagy_codes/03_mutation/fn_get_mutation_exclusive.R")

fn_draw_exclusive <- function(plot_ready, V1, V2, cancer_types,p_val,type){
  print("start fn_draw_exclusive --------------")
  p_val <- signif(p_val,3)
  plot_ready %>% 
    dplyr::mutate(mut = ifelse(mut >= 1,as.integer(1),mut)) %>%
    tidyr::spread(key = symbol, value = mut) -> rank_ready
  
  rank_ready %>%
    dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::any_vars(. != 0)) %>% 
    nrow() / nrow(rank_ready) -> diff
  
  fig_name <- paste(type, cancer_types, V1, V2, p_val, signif(diff, digits = 2),sep = "_")
  
  rank_ready %>%
    dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::any_vars(. >= 1)) %>%
    dplyr::arrange(dplyr::desc(get(V1, .)), dplyr::desc(get(V2, .))) %>%
    dplyr::pull(barcode) -> plot_rank
  
  rank_ready %>%
    dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::all_vars(. < 1)) %>%
    dplyr::arrange(dplyr::desc(get(V2, .)), dplyr::desc(get(V1, .))) %>%
    dplyr::pull(barcode) %>%
    c(plot_rank, .) -> plot_rank
  
  
  plot_ready %>%
    dplyr::mutate(mut = ifelse(mut >= 1,as.integer(1),mut)) %>%
    ggplot(aes(x = barcode, y = symbol, fill = as.factor(mut))) +
    geom_tile(color = "white") +
    scale_x_discrete(limits = plot_rank) +
    scale_y_discrete(limits = c(V1,V2)) +
    # facet_grid(symbol ~ .) +
    scale_fill_manual(
      name = "Mutation_type",
      limits = c("0", "1"),
      label = c("NA", "Mut"),
      values = c( "#D6D6D6", "#FF6A6A")) +
    labs(x = "", y = "", title = paste(V1, "and", V2, "SNV mutual exclusive in", cancer_types, ", p=", p_val)) +
    theme(
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      
      # strip.text.y = element_text(angle = 0,hjust = 0,size = 8),
      # strip.text.x = element_text(size = 8,angle = 90,vjust = 0),
      # strip.background = element_blank(),
      
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      
      panel.background = element_blank(),
      panel.spacing.y  = unit(0, "lines"),
      
      text = element_text(size = 8, color = "black")
    )  -> p
  
  # pdf(file.path(snv_path,"mutual_exclusive_allICPs-highGene",paste(fig_name,"pdf",sep = ".")),width = 4,height = 2)
  # showtext_begin()
  # p
  # showtext_end()
  # dev.off()
  ggsave(file.path(snv_path,"mutual_exclusive_allICPs-highGene",paste(fig_name,"pdf",sep=".")),device = "pdf",width=4,height=2)
  ggsave(file.path(snv_path,"mutual_exclusive_allICPs-highGene",paste(fig_name,"png",sep = ".")),device = "png",width = 4,height = 2)
  
  # graph2ppt(x = p,file = file.path(snv_path,"mutual_exclusive_allICPs-highGene",paste(fig_name,"pptx",sep = ".")),width = 4,height = 2)
}
fn_ex <- function(V1, V2, .data, cancer_types){
  # .data <- filter_cnv
  # V1 <- 'TP53'
  # V2 <- 'HLA-E'
  print("start fn_ex --------------")
  V1 <- as.character(V1)
  V2 <- as.character(V2)
  # print(V1)
  #  print(V2)
  .data %>% 
    dplyr::filter(symbol %in% c(V1, V2)) %>% 
    tidyr::gather(key = barcode, value = mut, -symbol) %>%
    dplyr::mutate(mut = ifelse(is.na(mut),0,mut)) %>%
    dplyr::mutate(mut = as.integer(mut)) -> plot_ready
  
  # get exclusive significance by combetExactTest
  plot_ready %>%
    # dplyr::filter(symbol %in% gene_rank$symbol) %>%
    tidyr::spread(key="barcode",value = "mut") %>%
    as.matrix() -> data.matrix
  rownames(data.matrix) <- data.matrix[,1]
  data.matrix <- data.matrix[,-1]
  data.matrix <- apply(data.matrix,2,as.numeric)
  
  result <- fn_get_mutation_pattern(data.matrix)
  
  # plot if pvalue is significant
  plot_ready %>%
    tidyr::spread(key = symbol, value = mut) %>% 
    dplyr::select(-barcode) -> .d
  .g_name <- colnames(.d)
  
  name <- paste(c(cancer_types, .g_name), collapse = "_")
  print(name)
  
  if (result$exclusive_pvalue <= 0.05) {
    fn_draw_exclusive(plot_ready, V1, V2, cancer_types,result$exclusive_pvalue,"Exclusive")
  }
  if (result$mutual_pvalue <= 0.05){
    fn_draw_exclusive(plot_ready, V1, V2, cancer_types,result$mutual_pvalue,"Mutual")
  }
  print("end fn_ex --------------")
  
  result %>%
    dplyr::mutate(te = name)
}
fn_mutal_exclusive <- function(cancer_types, ICP_SNV, highGene_SNV, cluster){
  ICP_SNV %>%
    tidyr::gather(-symbol,key="barcode",value="mut") %>%
    dplyr::mutate(mut = ifelse(is.na(mut),0,mut)) %>%
    dplyr::group_by(barcode) %>%
    dplyr::mutate(mut_overall = sum(mut)) %>%
    dplyr::mutate(mut_overall = ifelse(mut_overall >0,1,mut_overall)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-symbol,-mut) %>%
    dplyr::mutate(symbol="ICPs") %>%
    dplyr::rename("mut"="mut_overall") %>%
    unique() %>%
    tidyr::spread(key="barcode",value="mut") -> ICP_overall_mut
  
  expand.grid(highGene_SNV$symbol,ICP_overall_mut$symbol) -> .gene_pairs
  ICP_overall_mut %>%
    rbind(highGene_SNV) -> combine_data
  .gene_pairs %>% 
    # multidplyr::partition(cluster = cluster) %>%
    # multidplyr::cluster_library("magrittr") %>%
    # multidplyr::cluster_library("showtext") %>%
    # multidplyr::cluster_library("ggplot2") %>%
    # multidplyr::cluster_library("export") %>%
    # multidplyr::cluster_assign_value("fn_ex", fn_ex) %>%
    # multidplyr::cluster_assign_value("combine_data", combine_data) %>%
    # multidplyr::cluster_assign_value("cancer_types", cancer_types) %>%
    # multidplyr::cluster_assign_value("fn_draw_exclusive",fn_draw_exclusive) %>%
    # multidplyr::cluster_assign_value("snv_path",snv_path) %>%
    dplyr::mutate(rs = purrr::map2(Var1, Var2, .f = fn_ex, .data = combine_data, cancer_types = cancer_types)) %>% 
    # dplyr::collect() %>%
    dplyr::as_tibble() %>%
    dplyr::ungroup() %>%
    # dplyr::select(-PARTITION_ID) %>%
    dplyr::select(rs) %>% 
    tidyr::unnest() -> .gene_pairs_pval
  
  .gene_pairs_pval 
}

# cluster <- multidplyr::create_cluster(5)
ICP_highGene_snv %>% 
  # dplyr::filter(cancer_types == "UVM") %>%
  purrr::pmap(.f = fn_mutal_exclusive, cluster = cluster) %>% 
  dplyr::bind_rows() -> mutual_exclusive

mutual_exclusive %>%
  readr::write_tsv(file.path(snv_path,"mutual_exclusive_allICPs-highGene","allICPs_mutual_exclusive_highGeneSNV.tsv"))

mutual_exclusive %>% 
  tidyr::separate(col = te, into = c('cancer_types', 'g1', 'g2')) %>% 
  dplyr::mutate(exclusive_fdr = p.adjust(exclusive_pvalue, method = 'fdr')) %>%
  dplyr::mutate(mutual_fdr = p.adjust(mutual_pvalue, method = 'fdr')) -> mutual_exclusive.fdr
mutual_exclusive.fdr %>%
  readr::write_tsv(file.path(snv_path,"mutual_exclusive_allICPs-highGene","allICPs_mutual_exclusive_highGeneSNV.fdr.tsv"))
