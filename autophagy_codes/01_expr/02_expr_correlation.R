###########################################################################
# expression correlation between ligand and receptor in cancers ------------------------------

library(magrittr)
library(ggplot2)

# data path config
immune_path <- "/project/huff/huff/immune_checkpoint"
tcga_path <- file.path(immune_path,"/data/TCGA_data")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"

# load data
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
gene_list_expr <- readr::read_rds(path = file.path(expr_path, ".rds_03_a_gene_list_expr.rds.gz"))

# function to get correlation -----
fn_cor_a <- function(.data){
  .data %>%
    tidyr::gather(-symbol,-entrez_id,key="barcode",value="expr") %>%
    dplyr::as.tbl() %>%
    dplyr::filter(!is.na(expr)) %>%
    tidyr::nest(-symbol,-entrez_id) -> .tmp
  .tmp %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(spm = purrr::map(data,data_b=.tmp,fn_cor_b)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::ungroup() %>%
    dplyr::rename(Cor=estimate)-> .out
  return(.out)
}

fn_cor_b <- function(.data_a,data_b){
  res <- data.frame(estimate=NULL,statistic=NULL,p.value=NULL,symbol.x=NULL)
  for (i in 1:nrow(data_b)) {
    data_b$data[[i]] -> .tmp_b
    set.seed(12345)
    chrun <- runif(nrow(.tmp_b),min=0,max=0.001)
    .data_a %>%
      dplyr::inner_join(.tmp_b,by="barcode") %>%
      dplyr::mutate(expr.x = expr.x + chrun, expr.y = expr.y + chrun) -> .tmp_c
    broom::tidy(cor.test(.tmp_c$expr.x,.tmp_c$expr.y,method = c("spearman")),
                warning =function(e) 2 ,
                error=function(e) 1) %>%
      dplyr::mutate(symbol.x = data_b$symbol[i])-> tmp.spm
    rbind(res,tmp.spm) -> res
  }
  return(res)
}

gene_list_expr %>%
  # dplyr::filter(cancer_types=="KIRC") %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(spm = purrr::map(filter_expr,fn_cor_a)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-filter_expr) -> gene_list_expr_correlation

# trigle heatmap ----
CPCOLS <- c("red", "white", "blue")
gene_list_expr_correlation %>%
  tidyr::unnest() %>%
  dplyr::filter(p.value<=0.05 & abs(Cor)>=0.5) %>%
  dplyr::select(symbol,symbol.x,Cor) %>%
  tidyr::spread(key="symbol.x",value="Cor") -> draw_ready

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
draw_ready_up_tri <- get_upper_tri(draw_ready[,-1])

draw_ready_up_tri %>%
  dplyr::mutate(symbol = draw_ready$symbol) %>%
  tidyr::gather(-symbol,key="symbol.x",value="Cor") %>%
  dplyr::filter(!is.na(Cor)) %>%
  ggplot(aes(x=symbol,y=symbol.x)) +
  geom_tile(aes(fill=Cor),color = "white") +
  scale_fill_gradient2(
    name = "Correlation", # "Methylation diff (T - N)",
    low = CPCOLS[3],
    mid = CPCOLS[2],
    high = CPCOLS[1],
    limit = c(-1,1), space = "Lab"
  ) +
  theme_minimal()+ 
  theme(
    axis.text.x = element_text(angle = 30)
  ) +
  coord_fixed()

# correlation point plot ----
fn_cor_c <- function(.data) {
  .data %>%
    dplyr::select(barcode,symbol,expr) %>%
    tidyr::spread(key="symbol",value="expr") %>%
    as.data.frame() -> .tmp
  broom::tidy(cor.test(.tmp[,2], .tmp[,3],data = .data, method = "spearman"))
}
gene_list_expr %>%
  tidyr::unnest() %>%
  dplyr::filter(symbol %in% c("TNFSF4","TNFRSF4")) %>%
  tidyr::gather(-cancer_types,-symbol,-entrez_id,key="barcode",value="expr") %>%
  dplyr::filter(! is.na(expr)) -> ligand_receptor_data

ligand_receptor_data %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::group_by(cancer_types) %>%
  dplyr::mutate(spm = purrr::map(data,fn_cor_c)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> ligand_receptor_correlation

ligand_receptor_correlation %>%
  dplyr::filter(abs(estimate)>=0.2,p.value <=0.05) %>%
  dplyr::mutate(cancer_label = paste(cancer_types,"r=",signif(estimate,2),"p=",signif(p.value,2))) %>%
  dplyr::select(cancer_types,cancer_label) -> sig_cancers

ligand_receptor_data %>%
  dplyr::inner_join(sig_cancers,by="cancer_types") %>%
  dplyr::select(-entrez_id) %>%
  dplyr::mutate(expr = log2(expr)) %>%
  tidyr::spread(key="symbol",value="expr") %>%
  ggplot(aes(x=TNFSF4,y=TNFRSF4)) +
  geom_point(aes(color = cancer_label),size=0.5) +
  geom_smooth(aes(color = cancer_label),method = "lm") +
  facet_wrap(~cancer_label) +
  theme_classic() +
  theme(
    legend.position = "none"
  )
ggsave(file.path(out_path,"e_6_exp_profile","TNFSF4_TNFRSF4_cor.pdf"),width = 8,height = 4)
ggsave(file.path(out_path,"e_6_exp_profile","TNFSF4_TNFRSF4_cor.png"),width = 8,height = 4)
