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
fn_density_peak <-function(.name,.x){
  print(.name)
  # get peak and secondary peak
  set.seed(42)
  D1 <- .x$T_T
  if(length(unique(as.vector(D1)))<3){
    D1 <- D1 + runif(length(D1),0.001,0.002)
  }
  d1 <- density(D1)
  d1.ym <- which.max(d1$y)
  d1.ym.x <- d1$x[d1.ym]
  
  set.seed(42)
  D2 <- .x$I_T
  if(length(unique(as.vector(D2)))<3){
    D2 <- D2 + runif(length(D2),0.001,0.002)
  }
  d2 <- density(D2)
  d2.ym <- which.max(d2$y)
  d2.ym.x <- d2$x[d2.ym]
  tibble::tibble(Immune_peak = d2.ym.x, Tumor_peak = d1.ym.x)
}
# ggplot(.x, aes(I_T)) + geom_density() + geom_vline(xintercept = d2$x[peaks])
# ggplot(.x, aes(T_T)) + geom_density() + geom_vline(xintercept = d1$x[peaks1])

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
  dplyr::mutate(run = purrr::map2(sy_can,data,fn_density_peak)) %>% 
  dplyr::select(-data) %>%
  tidyr::unnest() -> gene_exp_density_peak_in_each_cancer

# draw peak value
gene_exp_density_peak_in_each_cancer %>%
  tidyr::separate(sy_can,c("symbol","cancer_types"),"_") %>%
  ggplot(aes(x=Tumor_peak,y=Immune_peak)) +
  geom_jitter(aes(color = cancer_types),size=1) +
  facet_wrap(~symbol) +
  scale_color_manual(
    name = "Cancer types",
    values = cancer21_color$color,
    limits = cancer21_color$cancer_types
  )+
  ylab("Expression peak in Immune cell set") +
  xlab("Expression peak in Tumor cell set") +
  theme_bw() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key.width = unit(0.15,"inches"),
    legend.key.height=unit(0.15,"inches"),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black",size = 12),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=12)
  )
ggsave(filename = paste("gene_expression_peak_value_in_cancers","png",sep="."),device = "png",path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_adjusted_exp","global_distribution"),width = 10,height = 8)
ggsave(filename = paste("gene_expression_peak_value_in_cancers","pdf",sep="."),device = "pdf",path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_adjusted_exp","global_distribution"),width = 10,height = 8)

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
fn_correlation <- function(.x){
  broom::tidy(cor.test(.x$exp,.x$predic_exp_from_celldata,method = "pearson")) 
}
gene_list_adjust_expr %>%
  tidyr::unnest() %>%
  dplyr::select(cancer_types,symbol,barcode,exp,predic_exp_from_celldata) %>%
  tidyr::nest(-cancer_types,-symbol) %>%
  dplyr::group_by(symbol,cancer_types) %>%
  dplyr::mutate(cor = purrr::map(data,fn_correlation)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::ungroup() -> tcgaExp_cellE_correlation

gene_list_adjust_expr %>%
  tidyr::unnest() %>%
  dplyr::select(cancer_types,symbol,barcode,exp,predic_exp_from_celldata) %>%
  tidyr::nest(-symbol) %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(cor = purrr::map(data,fn_correlation)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(cancer_types="Total") %>%
  dplyr::ungroup()-> tcgaExp_cellE_correlation.total

fn_cor_class <- function(.x){
  if(.x<(-0.7)){
    class <-"1" 
  }else if(.x>=(-0.7) & .x<(-0.5)){
    class <-"2" 
  }else if(.x>=(-0.5) & .x<(-0.3)){
    class <-"3" 
  }else if(.x>=(-0.3) & .x<0){
    class <-"4" 
  }else if(.x>=0 & .x<0.3){
    class <-"5" 
  }else if(.x>=0.3 & .x<0.5){
    class <-"6" 
  }else if(.x>=0.5 & .x<0.7){
    class <-"7" 
  }else if(.x>=0.7 & .x<1){
    class <-"8" 
  }
  class
}

tcgaExp_cellE_correlation %>%
  rbind(tcgaExp_cellE_correlation.total) %>%
  dplyr::filter(!is.na(estimate)) %>%
  dplyr::mutate(Cor_class = purrr::map(estimate,fn_cor_class)) %>%
  tidyr::unnest() %>%
  dplyr::mutate(p = ifelse(p.value>0.05,"1","2")) %>%
  ggplot(aes(x=cancer_types,y=symbol)) +
  geom_point(aes(size=p,fill=Cor_class),color="grey",shape=22) +
  scale_fill_manual(
    values=c("#104E8B","#0000EE","#1874CD", "#87CEFF","pink","pink2", "palevioletred2","#FF0000"),
    labels = c("[-1,-0.7)","[-0.7,-0.5)","[-0.5,-0.3)","[-0.3,0)","[0,0.3)","[0.3,0.5)","[0.5,0.7)","[0.7,1]")
  ) +
  scale_size_discrete(
    labels = c("p>0.05","p<=0.05"),
    name = NULL
  ) +
  scale_x_discrete(limit = c("Total",sort(tcgaExp_cellE_correlation$cancer_types %>% unique()))) +
  theme(
    panel.background = element_rect(fill = "white", 
                                    colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(vjust = 1, hjust = 1, size = 10,angle = 90),
    axis.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.text = element_text(colour = "black")
  )
ggsave(filename = paste("tcgaExp_cellE_correlation","pdf",sep="."),device = "pdf",path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_adjusted_exp","correlation"),width = 8,height = 10)
ggsave(filename = paste("tcgaExp_cellE_correlation","png",sep="."),device = "png",path = file.path(immune_path,"result_20171025/Exp_adjustment/verify_adjusted_exp","correlation"),width = 8,height = 10)
