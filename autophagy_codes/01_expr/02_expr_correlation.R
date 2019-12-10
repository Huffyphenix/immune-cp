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
gene_list <- readr::read_tsv(file.path(gene_list_path, "ICPs_all_info_class-new.tsv"))

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
  # dplyr::filter(p.value<=0.05 & abs(Cor)>=0.5) %>%
  dplyr::select(cancer_types,symbol,symbol.x,Cor) %>%
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
draw_ready_up_tri <- get_upper_tri(draw_ready[,-c(1:2)])

draw_ready_up_tri %>%
  dplyr::mutate(symbol = draw_ready$symbol,cancer_types = draw_ready$cancer_types) %>%
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

my_theme <-   theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_line(colour=NA),
  axis.text.y = element_text(size = 10,colour = "black"),
  axis.text.x = element_text(size = 10,colour = "black"),
  # legend.position = "none",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.background = element_blank(),
  legend.key = element_rect(fill = "white", colour = "black"),
  plot.title = element_text(size = 20),
  axis.text = element_text(colour = "black"),
  strip.background = element_rect(fill = "white",colour = "black"),
  strip.text = element_text(size = 10),
  text = element_text(color = "black")
)
# all correlation in one plot -----
# get receptor and ligand pairs 
gene_list %>%
  dplyr::filter(!is.na(Recepter_pairs)) %>%
  dplyr::filter(Recepter_pairs!="Pair 10") %>%
  dplyr::select(symbol,Recepter_pairs) %>%
  tidyr::nest(-Recepter_pairs) %>%
  dplyr::mutate(pairs = purrr::map(data,.f=function(.x){
    combn(.x$symbol,2) %>%
      t() %>%
      as.data.frame() %>%
      dplyr::as.tbl() %>%
      dplyr::rename("Gene1"="V1","Gene2" ="V2")
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> recepor_ligand.pairs

# filter correlation of receptor and ligand pairs 
gene_list_expr_correlation %>%
  dplyr::mutate(filter_cor = purrr::map(spm,.f=function(.x){
    .x %>%
      tidyr::nest(-symbol) %>%
      dplyr::mutate(filter = purrr::map2(symbol,data,.f=function(.x,.y){
        recepor_ligand.pairs %>%
          dplyr::filter(Gene1 == .x) %>%
          .$Gene2 -> .gene2
        recepor_ligand.pairs %>%
          dplyr::filter(Gene1 == .x) %>%
          .$Recepter_pairs %>%
          unique() -> .Recepter_pairs
        .y %>%
          dplyr::filter(symbol.x %in% .gene2) %>%
          dplyr::mutate(Recepter_pairs = .Recepter_pairs)
      })) %>%
      dplyr::select(-data) %>%
      tidyr::unnest()
  })) %>%
  dplyr::select(-spm) %>%
  tidyr::unnest() -> recepor_ligand.pairs.correlation
recepor_ligand.pairs.correlation %>%
  readr::write_tsv(file.path(out_path,"e_2_DE","tsv_06_receptor-ligand.correlation.tsv"))
## point plot
gene_list %>%
  dplyr::filter(!is.na(Recepter_pairs)) %>%
  dplyr::filter(Recepter_pairs!="Pair 10") %>%
  dplyr::select(symbol,Recepter_pairs) %>%
  tidyr::nest(-Recepter_pairs) %>%
  dplyr::mutate(Feature_symbol = purrr::map2(Recepter_pairs,data,.f=function(.x,.y){
    paste(.x,"(",paste(.y$symbol,collapse = "/"),")",sep="")
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest()  -> label.data

recepor_ligand.pairs.correlation %>%
  dplyr::mutate(log10P = ifelse(-log10(p.value)>=10,10,-log10(p.value))) %>% 
  dplyr::inner_join(label.data,by="Recepter_pairs") %>%
  ggplot(aes(y=log10P,x=Cor)) +
  geom_jitter(aes(color=Feature_symbol)) +
  facet_wrap(.~cancer_types) +
  my_theme +
  ggpubr::color_palette(palette = c("#FFB6C1", "#FFBBFF", "#0000FF", "#8B2323", "#CDAA7D", "#8EE5EE", "#76EE00",
                                    "#D2691E", "#8B008B", "#6E8B3D", "#FF3030", "#006400", "#FFD700", "#EE00EE")) +
  labs(x="Spearman correlation of expression between genes within receptor-ligand pairs",
       y=latex2exp::TeX("-log_{10} (P value)")) +
  theme(legend.title = element_blank())
ggsave(file.path(out_path,"e_2_DE","fig_06_receptor-ligand.correlation.pdf"),device = "pdf",height = 8,width = 12)
ggsave(file.path(out_path,"e_2_DE","fig_06_receptor-ligand.correlation.png"),device = "png",height = 8,width = 12)

## heat map plot
recepor_ligand.pairs.correlation -> recepor_ligand.pairs.correlation.heatmap
recepor_ligand.pairs.correlation.heatmap %>%
  dplyr::arrange(Recepter_pairs) %>%
  dplyr::select(symbol,symbol.x,Recepter_pairs) -> receptor_ligand.rank

recepor_ligand.pairs.correlation.heatmap <- within(recepor_ligand.pairs.correlation.heatmap,symbol<-factor(symbol,levels=unique(receptor_ligand.rank$symbol)))
with(recepor_ligand.pairs.correlation.heatmap,levels(symbol))

recepor_ligand.pairs.correlation.heatmap <- within(recepor_ligand.pairs.correlation.heatmap,symbol.x<-factor(symbol.x,levels=unique(receptor_ligand.rank$symbol.x)))
with(recepor_ligand.pairs.correlation.heatmap,levels(symbol.x))

recepor_ligand.pairs.correlation.heatmap %>%
  # dplyr::filter(cancer_types=="SKCM") %>%
  ggplot(aes(x=symbol,y=symbol.x)) +
  geom_tile(aes(fill=Cor),color="grey",size=0.2) +
  facet_wrap(.~cancer_types) +
  scale_fill_gradient2(
    name = "Spearman cor.",
    low = "skyblue3",
    mid = "white",
    high = "tomato3",
    limit = c(-1,1), 
    space = "Lab"
  ) +
  guides(fill=guide_colorbar(title.position="left")) +
  my_theme +
  theme(
    legend.title = element_text(angle = 90),
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
    axis.title = element_blank()
  )
ggsave(file.path(out_path,"e_2_DE","fig_07_receptor-ligand.correlation-heatmap.pdf"),device = "pdf",height = 15,width = 15)
ggsave(file.path(out_path,"e_2_DE","fig_07_receptor-ligand.correlation-heatmap.png"),device = "png",height = 15,width = 15)

# TNFRSF4 correlation point plot ----
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

save.image(file.path(out_path,"e_2_DE","receptor-ligand.correlation.rdata"))
