library(magrittr)
library(ggplot2)
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
stage_path <- file.path(out_path, "c_1_stage")
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")


clinical_stage <- 
  readr::read_rds(path = file.path(tcga_path,"pancan34_clinical_stage.rds.gz")) %>% 
  dplyr::filter(n >= 40) %>% 
  dplyr::select(-n)

gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))


# merge clinical and expr
expr_stage <- 
  gene_list_expr %>%
  dplyr::inner_join(clinical_stage, by = "cancer_types")

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
fun_expr_stage_merge <- function(filter_expr, stage){
  # merge clinical and expr data
  filter_expr %>% 
    dplyr::select(-entrez_id) %>% 
    tidyr::gather(key = barcode, value = expr, -symbol) %>% 
    dplyr::mutate(type = fun_tn_type(barcode)) %>% 
    dplyr::filter(type == "01") %>% 
    dplyr::mutate(barcode = fun_barcode(barcode)) %>% 
    dplyr::select(symbol, barcode, expr) -> expr_clean
  expr_clean %>% dplyr::inner_join(stage, by = "barcode") -> expr_stage_ready
} # merge stage and expr to clean
fn_get_order <- function(.d){
  .d %>% 
    dplyr::group_by(stage) %>% 
    dplyr::summarise(me = mean(expr)) %>% 
    .$me %>% rank() -> .d_m
  
  if(identical(.d_m, c(1,2,3,4))){
    return(1)
  } else if(identical(.d_m, c(4,3,2,1))){
    return(2)
  } else{
    return(3)
  }
}
fun_stage_test <- function(expr_stage_ready){
  expr_stage_ready %>% 
    tidyr::drop_na(expr) %>%
    dplyr::group_by(symbol, stage) %>% 
    dplyr::mutate(l = n() > 5) -> tl
  

  if(! all(tl$l)){ #filter out cancers with less than 10 samples in one of stage.
    return(tibble::tibble())
  } else{
    expr_stage_ready %>% 
      tidyr::drop_na(expr) %>%
      dplyr::group_by(symbol) %>% 
      dplyr::do(broom::tidy(oneway.test(log2(expr + 1) ~ stage, data = .))) %>%
      dplyr::ungroup() %>% 
      dplyr::mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
      dplyr::select(symbol, p.value, fdr) -> diff_pval
    
    expr_stage_ready %>% 
      tidyr::drop_na(expr) %>%
      tidyr::nest(-symbol) %>% 
      dplyr::mutate(order = purrr::map_dbl(data, .f = fn_get_order)) %>% 
      dplyr::select(- data) -> symbol_order

    diff_pval %>% 
      dplyr::inner_join(symbol_order, by = "symbol")
  }
} # stage_test


# expr_stage %>%
#   dplyr::filter(cancer_types == "BRCA") %>% 
#   dplyr::mutate(merged_clean = purrr::map2(filter_expr, stage, fun_expr_stage_merge)) %>%
#   dplyr::select(-filter_expr, -stage) %>%
#   dplyr::mutate(diff_pval = purrr::map(merged_clean, fun_stage_test)) %>%
#   tidyr::unnest(diff_pval, .drop = F) -> te

cl <- parallel::detectCores()
cluster <- multidplyr::create_cluster(floor(cl * 5 / 6))
expr_stage %>%
  multidplyr::partition(cluster = cluster) %>%
  multidplyr::cluster_library("magrittr") %>%
  multidplyr::cluster_assign_value("fun_barcode", fun_barcode)  %>%
  multidplyr::cluster_assign_value("fun_tn_type", fun_tn_type) %>% 
  multidplyr::cluster_assign_value("fun_expr_stage_merge", fun_expr_stage_merge) %>% 
  multidplyr::cluster_assign_value("fun_stage_test", fun_stage_test) %>% 
  multidplyr::cluster_assign_value("fn_get_order", fn_get_order) %>% 
  dplyr::mutate(merged_clean = purrr::map2(filter_expr, stage, fun_expr_stage_merge)) %>% 
  dplyr::select(-filter_expr, -stage) %>% 
  dplyr::mutate(diff_pval = purrr::map(merged_clean, fun_stage_test)) %>% 
  dplyr::collect() %>%
  dplyr::as_tibble() %>%
  dplyr::ungroup() %>%
  dplyr::select(-PARTITION_ID) %>% 
  tidyr::unnest(diff_pval, .drop = F) -> expr_stage_sig_pval
on.exit(parallel::stopCluster(cluster))

expr_stage_sig_pval %>% 
  readr::write_rds(path = file.path(stage_path, ".rds_03_b_stage_gene_expr.rds.gz"), compress = "gz")

expr_stage_sig_pval %>%
  dplyr::select(cancer_types,symbol,p.value,fdr) %>%
  readr::write_csv(file.path(stage_path,"03_b_stage_gene_expr.csv"))

expr_stage_sig_pval %>%
  dplyr::select(cancer_types,symbol,p.value,fdr) %>%
  dplyr::filter(fdr<=0.05) %>%
  readr::write_csv(file.path(stage_path,"03_b_stage_gene_fdr0.05.csv"))

# expr_stage_sig_pval <- readr::read_rds(path = file.path(stage_path, ".rds_03_b_stage_gene_expr.rds.gz"))
#--------------------------------------------------------

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
get_pattern <- function(p.value) {
  if(!is.na(p.value)){
    if(p.value < 0.05) {
      return(1)
    }else{
      return(0)
    }}else{
      return(0)
    }
}# get pattern

expr_stage_sig_pval %>% 
  dplyr::mutate(n = mapply(get_pattern,p.value)) %>% 
  dplyr::select(cancer_types, symbol,n) %>% 
  tidyr::spread(key = cancer_types, value = n) -> pattern

cancer_rank <- pattern %>% fun_rank_cancer()
gene_rank <- pattern %>% 
  fun_rank_gene() %>% 
  dplyr::left_join(gene_list, by = "symbol") %>% 
  dplyr::arrange(col, rank)
gene_rank$col<-gene_rank$col %>% as.character()
gene_rank$size<-gene_rank$size %>% as.character()

expr_stage_sig_pval %>% 
  ggplot(aes(x = cancer_types, y = symbol, color = cancer_types)) +
  geom_point(aes(size = -log10(p.value))) +
  scale_x_discrete(limit = cancer_rank$cancer_types) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_size_continuous(
    limit = c(-log10(0.05), 15),
    range = c(1, 6),
    breaks = c(-log10(0.05), 5, 10, 15),
    labels = c("0.05", latex2exp::TeX("$10^{-5}$"), latex2exp::TeX("$10^{-10}$"), latex2exp::TeX("$< 10^{-15}$")),
    name = "P-value"
  ) +
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
    axis.text.y = element_text(color = gene_rank$col,face=gene_rank$size),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(fill = "white", colour = "black")
  ) ->p;p
ggsave(
  filename = "fig_03_b_stage_sig_genes.pdf",
  plot = p,
  device = "pdf",
  width = 8,
  height =10,
  path = stage_path
)
readr::write_rds(
  p,
  path = file.path(stage_path, ".fig_03_b_stage_sig_genes.pdf.rds.gz"),
  compress = "gz"
)
#-------------------------------------
library("ggthemes","ggpubr")
fun_draw_boxplot <- function(cancer_types, merged_clean, symbol, p.value, fdr){
  # print(cancer_types)
  p_val <- signif(-log10(p.value), digits = 3)
  gene <- symbol
  fig_name <- paste(cancer_types, gene, p_val, "pdf", sep = ".")
  comp_list <- list(c("Stage I", "Stage IV"), c("Stage II", "Stage IV"), c("Stage III", "Stage IV"))
  merged_clean %>% 
    dplyr::filter(symbol == gene) %>% 
    dplyr::mutate(expr = log2(expr)) %>% 
    dplyr::arrange(stage) %>% 
    ggpubr::ggboxplot(x = "stage", y = "expr",  color = "stage", pallete = "jco"  ) +
    ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test") + 
    ggpubr::stat_compare_means(method = "anova", label.y = 14) +
    labs(x  = "", y = "Expression (log2 RSEM)", title = paste(gene, "expression stage change in", cancer_types)) +
    ggthemes::scale_color_gdocs() -> p
    ggsave(filename = fig_name, plot = p, path = file.path(stage_path, "boxplot"), width = 6, height = 6,  device = "pdf")
}
fun_draw_boxplot_filter <- function(cancer_types, merged_clean, symbol, p.value, fdr){
  p_val <- signif(-log10(p.value), digits = 3)
  gene <- symbol
  fig_name <- paste(cancer_types, gene, p_val, "pdf", sep = ".")
  comp_list <- list(c("Stage I", "Stage II"), c("Stage II", "Stage III"), c("Stage III", "Stage IV"))
  
  d <- merged_clean %>% 
    dplyr::filter(symbol == gene) %>% 
    dplyr::mutate( expr = log2(expr)) %>% 
    dplyr::arrange(stage)
    
  p_list <- ggpubr::compare_means(expr ~ stage, data = d, method = "t.test")
  
  if(p_list %>% .$p %>% .[3] < 0.05){
    d %>% 
      ggpubr::ggboxplot(x = "stage", y = "expr",  color = "stage", pallete = "jco"  ) +
      ggpubr::stat_compare_means(comparisons = comp_list, method = "t.test") + 
      ggpubr::stat_compare_means(method = "anova", label.y = 14) +
      labs(x  = "", y = "Expression (log2 RSEM)", title = paste(gene, "expression subtype change in", cancer_types)) +
      ggthemes::scale_color_gdocs() -> p
      ggsave(filename = fig_name, plot = p, path = file.path(stage_path, "boxplot_stageI_vs_IV"), width = 6, height = 6,  device = "pdf")
  } else{
    print(fig_name)
  }
}


expr_stage_sig_pval %>% dplyr::select(-order)%>% purrr::pwalk(.f = fun_draw_boxplot)
# expr_stage_sig_pval %>% purrr::pwalk(.f = fun_draw_boxplot_filter)


save.image(file = file.path(stage_path, ".rda_03_b_stage_gene_expr.rda"))
load(file = file.path(stage_path, ".rda_03_b_stage_gene_expr.rda"))
