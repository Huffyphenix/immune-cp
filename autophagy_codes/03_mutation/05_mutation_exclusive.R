############################################
# get mutual exclusive mutation between ICPs and genes with high mutation rate in cancers.
############################################

library(magrittr)
library(ggplot2)
out_path <- "/home/huff/project/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")


# LOAD DATA 
ICP_highGene_snv <- readr::write_rds(file.path(snv_path,"ICP_highGene_snv.rds.gz"))

# GET EXCLUSIVE MUTATION BETWEEN ICP AND HIGH MUTATED GENES IN EACH CANCERS ------
library(showtext)
font_add("Arial","ARIAL.TTF") # sudo进入container，将/home/huff/fonts/中的字体拷贝到/usr/share/fonts/,then do:fc-cache

library(export)
# library(Cairo)
fn_draw_exclusive <- function(plot_ready, V1, V2, cancer_types,p_val){
  plot_ready %>% 
    tidyr::spread(key = symbol, value = mut) -> rank_ready
  
  rank_ready %>%
    dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::any_vars(. != 0)) %>% 
    nrow() / nrow(rank_ready) -> diff
  
  fig_name <- paste(cancer_types, p_val, signif(diff, digits = 2), V1, V2, sep = "_")
  
  rank_ready %>%
    dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::any_vars(. >= 1)) %>%
    dplyr::arrange(dplyr::desc(get(V1, .)), dplyr::desc(get(V2, .))) %>%
    dplyr::pull(barcode) -> plot_rank
  
  rank_ready %>%
    dplyr::filter_if(.predicate = is.integer, .vars_predicate = dplyr::all_vars(. < 1)) %>%
    dplyr::arrange(dplyr::desc(get(V2, .)), dplyr::desc(get(V1, .))) %>%
    dplyr::pull(barcode) %>%
    c(plot_rank, .) -> plot_rank
  
  showtext.begin()
  plot_ready %>%
    dplyr::mutate(mut = ifelse(mut>=2,2,mut)) %>%
    ggplot(aes(x = barcode, y = factor(1), fill = as.factor(mut))) +
    geom_tile(color="white",width=0.9) +
    scale_x_discrete(limits = plot_rank) +
    facet_grid(symbol ~ .) +
    scale_fill_manual(
      name="Mutation_type",
      limits = c("0", "1", "2"),
      label=c("NA", "1", ">=2"),
      values= c( "#D6D6D6", "#FF6A6A","#FF0000")) +
    labs(x = "", y = "", title = paste(V1, "and", V2, " SNV mutual exclusive in", cancer_types, ", p=", p_val))+
    theme(
      axis.text=element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      
      strip.text.y = element_text(angle =0,hjust=0,size=11),
      strip.text.x = element_text(size=11,angle = 90,vjust = 0),
      strip.background = element_blank(),
      
      legend.title= element_blank(),
      legend.text = element_text(size=11),
      legend.position = "bottom",
      
      panel.background = element_blank(),
      panel.spacing.y  = unit(0, "lines"),
      
      text = element_text(family = "Arial")
    )  -> p
  showtext.end()
  ggsave(file.path(snv_path,"mutual_exclusive",paste(fig_name,"pdf",sep=".")),device = "pdf",width=6,height=2)
  ggsave(file.path(snv_path,"mutual_exclusive",paste(fig_name,"png",sep=".")),device = "png",width=6,height=2)
  
  graph2ppt(x=p,file=file.path(snv_path,"mutual_exclusive",paste(fig_name,"pptx",sep=".")),width=6,height=2)
}
fn_ex <- function(V1, V2, .data, cancer_types){
  # .data <- filter_cnv
  # V1 <- 'GNAQ'
  # V2 <- 'BTN2A1'
  V1<-as.character(V1)
  V2<-as.character(V2)
 print(V1)
  print(V2)
  .data %>% 
    dplyr::filter(symbol %in% c(V1, V2)) %>% 
    tidyr::gather(key = barcode, value = mut, -symbol) %>%
    dplyr::mutate(mut=ifelse(is.na(mut),0,mut)) %>%
    dplyr::mutate(mut=as.integer(mut)) -> plot_ready
  plot_ready %>%
    tidyr::spread(key = symbol, value = mut) %>% 
    dplyr::select(-barcode) -> .d
  .g_name <- colnames(.d)
  
  name <- paste(c(cancer_types, .g_name), collapse = "_")
  print(name)
  
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(. == 0)) %>% 
    nrow() -> nn
  .d %>%
    dplyr::filter_all(.vars_predicate = dplyr::all_vars(. != 0)) %>% 
    nrow()-> aa
  .d %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. != 0)) %>% 
    dplyr::filter_all(.vars_predicate = dplyr::any_vars(. == 0)) -> .d_an
  
  sum(.d_an %>% dplyr::pull(1) != 0) -> an
  sum(.d_an %>% dplyr::pull(2) != 0) -> na
  c(nn = nn, an = an, na = na, aa = aa) %>% 
    cometExactTest::comet_exact_test(mutmatplot = F) -> p_val
  
  if(p_val <= 0.05){
    fn_draw_exclusive(plot_ready, V1, V2, cancer_types,p_val)
  }
  
  tibble::tibble(te = name, nn = nn, an = an ,na = na, aa = aa, p_val = p_val)
}
fn_mutal_exclusive <- function(cancer_types, ICP_SNV, highGene_SNV, cluster){
  expand.grid(highGene_SNV$symbol,ICP_SNV$symbol) -> .gene_pairs
  ICP_SNV %>%
    rbind(highGene_SNV) -> combine_data
  .gene_pairs %>% 
    multidplyr::partition(cluster = cluster) %>%
    multidplyr::cluster_library("magrittr") %>%
    multidplyr::cluster_assign_value("fn_ex", fn_ex) %>%
    multidplyr::cluster_assign_value("combine_data", combine_data) %>%
    multidplyr::cluster_assign_value("cancer_types", cancer_types) %>%
    dplyr::mutate(rs = purrr::map2(Var1, Var2, .f = fn_ex, .data = combine_data, cancer_types = cancer_types)) %>% 
    dplyr::collect() %>%
    dplyr::as_tibble() %>%
    dplyr::ungroup() %>%
    dplyr::select(-PARTITION_ID) %>%
    dplyr::select(rs) %>% 
    tidyr::unnest() %>% 
    tidyr::separate(col = te, into = c('cancer_types', 'g1', 'g2')) -> .gene_pairs_pval
  
  .gene_pairs_pval %>% 
    dplyr::mutate(fdr = p.adjust(p_val, method = 'fdr'))
}

cluster <- multidplyr::create_cluster(5)
ICP_highGene_snv %>% 
  # dplyr::filter(cancer_types == "UVM") %>%
  purrr::pmap(.f = fn_mutal_exclusive, cluster = cluster) %>% 
  dplyr::bind_rows() -> mutual_exclusive
parallel::stopCluster(cluster)

mutual_exclusive %>%
  readr::write_tsv(file.path(snv_path,"mutual_exclusive","ICP_mutual_exclusive_highGeneSNV.tsv"))