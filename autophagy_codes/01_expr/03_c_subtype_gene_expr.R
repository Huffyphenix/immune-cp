library(magrittr)
library(ggplot2)
basic_path <- "/home/huff/project/immune_checkpoint"
out_path <- file.path(basic_path,"result_20171025")
subtype_path <- file.path(out_path, "c_2_subtype/")
tcga_path <- file.path(basic_path,"data/TCGA_data")
expr_path <- file.path(out_path,"expr_rds")

clinical_subtype <- 
  readr::read_rds(path = file.path(tcga_path,"pancan34_clinical_subtype.rds.gz")) %>% 
  dplyr::select(-n)

gene_list_path <- file.path(basic_path,"checkpoint/20171021_checkpoint")
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))


# merge clinical and expr
expr_subtype <- 
  gene_list_expr %>%
  dplyr::inner_join(clinical_subtype, by = "cancer_types")

fun_barcode <- function(.b){
  stringr::str_sub(
    string = .b,
    start = 1,
    end = 12
  )
} #get short barcode from long barcode
fun_tn_type <- function(.b){
  type <- .b %>% 
    stringr::str_split(pattern = "-", simplify = T) %>% 
    .[, 4] %>% 
    stringr::str_sub(1, 2)
} # get tumor and normal info
fun_expr_subtype_merge <- function(filter_expr, subtype){
  # merge clinical and expr data
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type == "01") %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(symbol, barcode, expr)  %>% 
    dplyr::inner_join(subtype, by = "barcode") -> expr_subtype_ready
} 
fun_subtype_test <- function(expr_subtype_ready){
  expr_subtype_ready %>% 
    tidyr::drop_na(expr) %>%
    dplyr::group_by(symbol) -> d
  
  d %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(symbol, subtype) %>% 
    dplyr::mutate(l = n() > 5) -> tl
  
  if(! all(tl$l)){
    return(tibble::tibble())
  } else{
    d %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct(subtype) %>% 
      .$subtype %>% 
      length() -> n_subtype
    # d %>%
    #   #dplyr::ungroup() %>%
    #   dplyr::group_by(symbol,subtype) %>%
    #   dplyr::do(
    #     me_exp=mean(.$expr) 
    #   ) %>%
    #   dplyr::ungroup()
    #n_subtype == 2 t.test
    #n_subtype >3 anova
    if(n_subtype == 2){
      d %>% 
        dplyr::do(broom::tidy(t.test(log2(expr + 1) ~ subtype, data = .))) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
        dplyr::select(symbol, p.value, fdr)  #%>% 
        #dplyr::filter(p.value < 0.01, fdr < 0.1)
    } else{
      d %>% 
        dplyr::do(broom::tidy(oneway.test(log2(expr + 1) ~ subtype, data = .))) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
        dplyr::select(symbol, p.value, fdr)  #%>% 
        #dplyr::filter(p.value < 0.01, fdr < 0.1) 
      }
  }
}

# expr_subtype %>% 
#   dplyr::mutate(merged_clean = purrr::map2(filter_expr, subtype, fun_expr_subtype_merge)) %>% 
#   dplyr::select(-filter_expr, -subtype) %>% 
#   dplyr::mutate(diff_pval = purrr::map(merged_clean, fun_subtype_test)) %>% 
#   tidyr::unnest(diff_pval, .drop = F) -> expr_subtype_sig_pval

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
expr_subtype %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type) %>% 
  multidplyr::cluster_assign_value("fun_expr_subtype_merge", fun_expr_subtype_merge) %>% 
  multidplyr::cluster_assign_value("fun_subtype_test", fun_subtype_test) %>% 
  dplyr::mutate(merged_clean = purrr::map2(filter_expr, subtype, fun_expr_subtype_merge)) %>% 
  dplyr::select(-filter_expr, -subtype) %>% 
  dplyr::mutate(diff_pval = purrr::map(merged_clean, fun_subtype_test)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  tidyr::unnest(diff_pval, .drop = F) -> expr_subtype_sig_pval
on.exit(parallel::stopCluster(cluster))

expr_subtype %>%
  dplyr::mutate(merged_clean = purrr::map2(filter_expr, subtype, fun_expr_subtype_merge)) %>% 
  dplyr::select(-filter_expr, -subtype) ->test

expr_subtype_sig_pval %>% 
  readr::write_rds(path = file.path(subtype_path, ".rds_03_c_subtype_gene_expr.rds.gz"), compress = "gz")
expr_subtype_sig_pval <- readr::read_rds(path = file.path(subtype_path, ".rds_03_c_subtype_gene_expr.rds.gz"))

expr_subtype_sig_pval %>% 
  #filter fdr<=0.05
  dplyr::filter(fdr<=0.05) %>%
  dplyr::select(cancer_types,symbol,p.value,fdr) %>%
  readr::write_csv(path = file.path(subtype_path, "03_c_subtype_gene_fdr0.05.csv"))

#----------------------------------------------------

fun_rank_cancer <- function(pattern){
  pattern %>% 
    dplyr::summarise_if(.predicate = is.numeric, dplyr::funs(sum(., na.rm = T))) %>%
    tidyr::gather(key = cancer_types, value = rank) %>%
    dplyr::arrange(dplyr::desc(rank))
} #get cancer rank
fun_rank_gene <- function(pattern){
  pattern %>% 
    dplyr::rowwise() %>%
    dplyr::do(
      symbol = .$symbol,
      rank =  unlist(.[-1], use.names = F) %>% sum(na.rm = T)
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest() %>%
    dplyr::arrange(rank)
} # get gene rank

expr_subtype_sig_pval %>% 
  dplyr::select(cancer_types, symbol) %>% 
  dplyr::mutate(n = 1) %>% 
  tidyr::spread(key = cancer_types, value = n) -> pattern

cancer_rank <- pattern %>% fun_rank_cancer()
gene_rank <- pattern %>% 
  fun_rank_gene() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::mutate(color = plyr::revalue(type, replace = c("Lysosome" = "black", "Autophagy" = "red")))


expr_subtype_sig_pval %>% 
  ggplot(aes(x = cancer_types, y = symbol, color = cancer_types)) +
  geom_point(aes(size = -log10(p.value))) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(name = "P-value") +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed"),
    panel.grid.major = element_line(
      colour = "grey",
      linetype = "dashed",
      size = 0.2
    ),
    
    axis.title = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(color = gene_rank$color),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) -> p
ggsave(
  filename = "fig_03_c_subtype_sig_genes.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height = 25,
  path = subtype_path
)
# readr::write_rds(
#   p,
#   path = file.path(subtype_path, ".fig_03_c_subtype_sig_genes.pdf.rds.gz"),
#   compress = "gz"
# )

#-----------------------------------------------------------
library(export)
compare_list<-function(x){
  res<-list()
  for(i in 2:length(x)-1){
    res<-c(res,list(c(x[i],x[i+1])))
  }
  return(res)
}
library("ggthemes","ggpubr")

fun_draw_boxplot <- function(cancer_types, merged_clean, symbol, p.value, fdr){
  # print(cancer_types)
  p_val <- signif(-log10(p.value), digits = 3)
  gene <- symbol
  fig_name <- paste(cancer_types, gene, p_val, sep = "_")
  # comp_list <- list(c("Stage I", "Stage II"), c("Stage II", "Stage III"), c("Stage III", "Stage IV"))

  merged_clean %>% 
    dplyr::select(subtype) %>% 
    unique() %>% .[,1] %>% sort() %>%
    compare_list() ->comp_list
  
  merged_clean %>% 
    dplyr::filter(symbol == gene) %>% 
    dplyr::mutate( expr = log2(expr)) %>% 
    dplyr::filter(expr != "-Inf") %>%
    dplyr::select(expr) %>%
    max() ->max_exp
  merged_clean %>% 
    dplyr::filter(symbol == gene) %>% 
    dplyr::mutate( expr = log2(expr)) %>% 
    dplyr::filter(expr != "-Inf") %>%
    dplyr::select(expr) %>%
    min() ->min_exp
  merged_clean %>% 
    dplyr::filter(symbol == gene) %>% 
    dplyr::mutate(expr = log2(expr)) %>% 
    dplyr::arrange(subtype) %>% 
    ggpubr::ggviolin(x = "subtype", y = "expr",  fill = "subtype", pallete = "jco",alpha = 0.5) +
    theme(legend.position = "none") +
    ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test") + 
    # ylim(min_exp-1,max_exp+1) +
    #ggpubr::stat_compare_means(method = "kruskal.test",label.y = max_exp+6,label.x=1,label.sep = "\n") +
    # geom_text(label = paste("Oneway.test:" ,"\n","p=",signif(p.value,3),sep=""),x=1,y=max_exp+5,size=3)+
    labs(x  = "", y = "Expression (log2 RSEM)",title=paste(cancer_types,gene,sep=", ")) +
    # geom_text(label=paste(cancer_types,gene,sep="\n"),x =title.pos.x,y=max_exp+6,size=4,colour="black")+
    ggthemes::scale_fill_gdocs() -> p
  ggsave(filename = paste(fig_name,".pdf",sep = ""), plot = p, path = file.path(subtype_path, "boxplot"), width = 4, height = 3,  device = "pdf")
  ggsave(filename = paste(fig_name,".png",sep = ""), plot = p, path = file.path(subtype_path, "boxplot"), width = 4, height = 3,  device = "png")
  
  export::graph2ppt(x = p,file = file.path(subtype_path, "boxplot",paste(fig_name,".pptx",sep = "")), width = 4, height = 3,font = "Times New Roman")
}

expr_subtype_sig_pval %>% 
  dplyr::filter(fdr<=0.05) %>%
  purrr::pwalk(.f = fun_draw_boxplot)


save.image(file = file.path(subtype_path, ".rda_03_c_subtype_gene_expr.rda"))
load(file = file.path(subtype_path, ".rda_03_c_subtype_gene_expr.rda"))
