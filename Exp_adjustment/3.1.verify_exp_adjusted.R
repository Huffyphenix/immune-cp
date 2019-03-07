#####################################################
# verify the expression adjusted by our method
#####################################################
library(magrittr)
library(ggplot2)

immune_path <- "/project/huff/huff/immune_checkpoint"
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
# load data ---------------------------------------------------------------

gene_list_adjust_expr <-  readr::read_rds(file.path(gene_list_path,"pancan21.ICP.exp_adjust.by_cell_ratio.rds.gz"))

# global distribution, the peak value of density -----------------------

library(mixtools)
require(graphics)
require(pastecs)
fn_density_peak_pastecs <-function(.name,.x){
  print(.name)
  # get peak and secondary peak
  set.seed(42)
  
  D1 <- .x$T_T
  d1 <- density(D1, bw = "sj")
  ts_y1<-ts(d1$y)
  tp1<-turnpoints(ts_y1)
  peaks1 <- grep("TRUE",tp1$peaks)
  tp1$pos[peaks1]
  d1$x[peaks1]
  
  set.seed(42)
  
  D2 <- .x$I_T
  
  d2 <- density(D2, bw = "sj")
  ts_y<-ts(d2$y)
  tp<-turnpoints(ts_y)
  peaks <- grep("TRUE",tp$peaks)
  tp$pos[peaks]
  d2$x[peaks]
}
ggplot(.x, aes(I_T)) + geom_density() + geom_vline(xintercept = d2$x[peaks])
ggplot(.x, aes(T_T)) + geom_density() + geom_vline(xintercept = d1$x[peaks1])

gene_list_adjust_expr %>%
  tidyr::unnest() %>%
  dplyr::filter(!is.na(T_T)) %>%
  dplyr::filter(!is.na(I_T)) %>%
  dplyr::select(cancer_types,symbol,barcode,T_T,I_T) %>%
  dplyr::mutate(T_T=ifelse(T_T==0,0.01,T_T)) %>%
  dplyr::mutate(I_T=ifelse(I_T==0,0.01,I_T)) %>%
  dplyr::mutate(T_T=log2(T_T)) %>%
  dplyr::mutate(I_T=log2(I_T)) %>%
  tidyr::nest(-symbol,-cancer_types) %>%
  tidyr::unite(sy_can,symbol,cancer_types,sep = "_") %>%
  dplyr::mutate(run = purrr::map2(sy_can,data,fn_density_peak)) %>% 太差，还是只看第一个峰吧！
  dplyr::select(-data) %>%
  tidyr::unnest() -> gene_exp_density_peak_in_each_cancer


# distribution of each genes in each cancer -----------------------

cancer_color <- readr::read_tsv(file.path("/data/shiny-data/GSCALite","02_pcc.tsv"))
cancer_color %>%
  dplyr::filter(cancer_types %in% unique(gene_list_adjust_expr$cancer_types)) -> cancer21_color

gene_list_adjust_expr %>%
  tidyr::unnest() %>%
  dplyr::filter(!is.na(T_T)) %>%
  dplyr::filter(!is.na(I_T)) %>%
  dplyr::select(cancer_types,symbol,barcode,T_T,I_T) %>%
  dplyr::mutate(T_T=ifelse(T_T==0,0.01,T_T)) %>%
  dplyr::mutate(I_T=ifelse(I_T==0,0.01,I_T)) %>%
  dplyr::mutate(T_T=log2(T_T)) %>%
  dplyr::mutate(I_T=log2(I_T)) %>%
  tidyr::nest(-symbol) %>%
  dplyr::mutate(run = purrr::map2(symbol,data,fun_draw_distribution))

fun_draw_distribution <- function(symbol,plot_ready){
  # cancer annotation on scatter plot
  plot_ready %>% 
    dplyr::group_by(cancer_types) %>%
    dplyr::mutate(middle_I = mean(I_T)) %>%
    dplyr::mutate(middle_T = mean(T_T)) %>%
    dplyr::select(cancer_types,middle_T,middle_I) %>%
    unique() %>%
    dplyr::ungroup() -> anno
  
  sp<- ggscatter(plot_ready, x = "T_T", y = "I_T",
                 color = "cancer_types", palette = "jco",
                 size = 1, alpha = 0.6)+
    border() +
    scale_color_manual(
      values = cancer21_color$color,
      limits = cancer21_color$cancer_types
    )+
    geom_text(data=anno,aes(x=middle_T,y=middle_I,label=cancer_types)) +
    ylab(paste(symbol,"expression in Immune cell set")) +
    xlab(paste(symbol,"expression in Tumor cell set")) +
    theme(
      axis.text = element_text(colour = "black")
    )
  # Marginal density plot of x (top panel) and y (right panel)
  xplot <- ggdensity(plot_ready, "T_T", fill = "cancer_types",
                     palette = "jco")+
    scale_fill_manual(
      values = cancer21_color$color,
      limits = cancer21_color$cancer_types
    )
  yplot <- ggdensity(plot_ready, "I_T", fill = "cancer_types", 
                     palette = "jco")+
    rotate()+
    scale_fill_manual(
      values = cancer21_color$color,
      limits = cancer21_color$cancer_types
    )
  # Cleaning the plots
  yplot <- yplot + clean_theme() 
  xplot <- xplot + clean_theme()
  # Arranging the plot
  ggarrange(xplot, NULL, sp, yplot, 
            ncol = 2, nrow = 2,  align = "hv", 
            widths = c(4, 1), heights = c(1, 3),
            legend = "top",
            common.legend = TRUE)
  filename <- paste(symbol,"exp_distribution_in_cancers",sep="_")
  ggsave(filename = paste(filename,"png",sep="."),device = "png",path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_adjusted_exp","distribution"),width = 10,height = 8)
  ggsave(filename = paste(filename,"pdf",sep="."),device = "pdf",path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_adjusted_exp","distribution"),width = 10,height = 8)
  return(1)
}


# adjusted expression comfirm ---------------------------------------------


