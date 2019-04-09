# merge stage, subtype and survival results into one


# path --------------------------------------------------------------------

data_path <- "/home/huff/project/immune_checkpoint/result_20171025"
gene_list_path <- "/home/huff/project/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
ICP_expr_pattern <- readr::read_tsv(file.path(data_path,"ICP_exp_patthern","manual_edit_2_ICP_exp_pattern_in_immune_tumor_cell.tsv"))
fn_site_color <- function(.n,.x){
  print(.n)
  if(.x=="Mainly_Tumor"){
    "red"
  }else if(.x=="Mainly_Immune"){
    "Blue"
  }else if(.x=="Both"){
    c("#9A32CD")
  }else{
    "grey"
  }
}
gene_list %>%
  dplyr::inner_join(ICP_expr_pattern,by="symbol") %>%
  dplyr::rename("Exp_site"="Exp site") %>%
  dplyr::mutate(Exp_site=ifelse(is.na(Exp_site),"N",Exp_site)) %>%
  dplyr::mutate(site_col = purrr::map2(symbol,Exp_site,fn_site_color)) %>%
  tidyr::unnest() -> gene_list
# load result -------------------------------------------------------------

stage_res <- readr::read_csv(file.path(data_path,"c_1_stage","03_b_stage_gene_fdr0.05.csv")) %>%
  dplyr::rename("stage_fdr"="fdr")
subtype_res <- readr::read_csv(file.path(data_path,"c_2_subtype","03_c_subtype_gene_fdr0.05.csv")) %>%
  dplyr::rename("subtype_fdr"="fdr")
survival_res <- readr::read_csv(file.path(data_path,"c_3_survival","PFS_c_3_survival_genelist_pval.csv"))


# merge three data into one -----------------------------------------------

survival_res %>%
  dplyr::select(-estimate,-p.value) %>%
  dplyr::left_join(stage_res,by=c("cancer_types","symbol")) %>%
  dplyr::left_join(subtype_res,by=c("cancer_types","symbol")) %>%
  dplyr::select(-p.value.x,-p.value.y) %>%
  dplyr::mutate(kmp = ifelse(kmp < 0.05, status, NA)) %>%
  dplyr::mutate(stage_fdr = ifelse(!is.na(stage_fdr),"Stage",NA)) %>%
  dplyr::mutate(subtype_fdr = ifelse(!is.na(subtype_fdr),"Subtype",NA)) %>%
  dplyr::select(-status) %>%
  tidyr::gather(-cancer_types,-symbol,key="key",value="value") -> ready_draw

ready_draw %>%
  # dplyr::filter(cancer_types %in% c("BLCA","BRCA")) %>%
  dplyr::filter(!is.na(value)) %>% 
  dplyr::select(symbol) %>%
  unique() %>%
  dplyr::inner_join(gene_list,by="symbol") %>%
  dplyr::arrange(`Exp_site`,functionWithImmune) -> gene_rank
  
library(ggplot2)
ready_draw %>%
  # dplyr::filter(cancer_types %in% c("BLCA","BRCA")) %>%
  dplyr::filter(!is.na(value)) %>%
  ggplot(aes(x=key,y=symbol,fill = value))+
  geom_tile(color="grey",aes(height=0.9,width=0.9)) +
  facet_wrap(~cancer_types,nrow=1) +
  scale_fill_manual(
    values = c("red3","dodgerblue2","violet","gold","white")
  ) +
  scale_y_discrete(limit = gene_rank$symbol) +
  theme(
    legend.key = element_rect(fill = "white", colour = "black"),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text = element_text(colour = "black",size = 8),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = gene_rank$site_col),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill="white",color="black"),
    strip.text = element_text(color="black",size=8),
    axis.title = element_blank(),
    panel.border = element_rect(colour = "black",fill=NA,size=0.5),
    panel.background = element_blank(),
    panel.spacing.x = unit(0,"lines")
  ) -> p1

gene_rank %>%
  dplyr::mutate(fun = "functionWithImmune") %>%
  ggplot(aes(y=symbol,x=fun)) +
  geom_tile(aes(fill = functionWithImmune),color="grey",size=1) +
  scale_y_discrete(limit = gene_rank$symbol) +
  scale_fill_manual(
    name = "Immune Checkpoint",
    values = c("#1C86EE", "#EE3B3B", "#EE7600")
  ) +
  theme(
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_blank(),
    # panel.grid = element_line(colour = "grey", linetype = "dashed"),
    # panel.grid.major = element_line(
    #   colour = "grey",
    #   linetype = "dashed",
    #   size = 0.2
    # ),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 10),
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    plot.margin=unit(c(0,0,0,-0), "cm")
  ) -> p2
ggarrange(p2,p1,
          ncol = 2, nrow = 1,  align = "hv", 
          widths = c(1, 20), legend = "bottom",
          common.legend = FALSE) -> p;p

ggsave(file.path(data_path,"clinical_result","clinical_merge.pdf"),device = "pdf",height = 8,width = 12)
ggsave(file.path(data_path,"clinical_result","clinical_merge.png"),device = "png",height = 8,width = 12)
