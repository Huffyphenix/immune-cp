############## ICP exp in inhibitory immune cell ###############
#################### FANTOM data ###############################
library(curl)
library(ggbeeswarm)
library(magrittr)

# data path ---------------------------------------------------------------
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
result_path <- file.path(immune_path,"result_20171025","ICP_exp_patthern-byratio")
fantom_path <- file.path(basic_path,"data/FANTOM5/extra")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")


# load image --------------------------------------------------------------

load(file.path(result_path,"pattern_validation","FANTOM5.validation.Rdata"))

# load data ---------------------------------------------------------------

ICP_fantom.gene_exp.Immune_cell.combine <-
  readr::read_rds(file.path(immune_path,"genelist_data","FANTOM5","ICP_fantom.gene_exp.cell_line.Immune_cell.raw.exp.rds.gz")) %>%
  dplyr::filter(Group != "Stromal Cell") 
  # dplyr::filter(`Characteristics[Tissue]` != "blood")

fantom_sample_info <- readr::read_tsv(file.path(fantom_path,"HumanSamples2.0.classification.txt")) %>%
  # dplyr::filter(`Characteristics [Category]` %in% c("primary cells")) %>%
  dplyr::mutate(sample = toupper(paste(curl_escape(`Charateristics [description]`),`Source Name`,sep="."))) %>%
  dplyr::select(sample,`Charateristics [description]`,`Characteristics[Tissue]`,`Characteristics [Cell type]`,`Characteristics [Category]`)  %>%
  dplyr::filter(sample %in% unique(ICP_fantom.gene_exp.Immune_cell.combine$sample))


# ICP exp in each cell types ----------------------------------------------

gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)

ICP_exp_site <- readr::read_tsv(file.path(result_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv"))
gene_list_exp_site <- readr::read_tsv(file.path(result_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T)`) %>%
  dplyr::inner_join(gene_list,by="symbol") 

inhibitory_immune_cell<- c("monocyte","neutrophil","monoblast")
ICP_fantom.gene_exp.Immune_cell.combine %>%
  dplyr::inner_join(fantom_sample_info,by="sample") %>%
  dplyr::mutate(`Characteristics[Tissue].x` = ifelse(is.na(`Characteristics[Tissue].x`),"notsure",`Characteristics[Tissue].x`)) %>% # cell lines tissue type, not include immune cell. Cause replaced by "PrimaryCell" in ICP_fantom.gene_exp.Immune_cell.combine
  dplyr::filter(`Characteristics[Tissue].x`!="notsure") %>%
  dplyr::mutate(Group = ifelse(`Characteristics[Tissue].x`=="blood","blood",Group)) %>%
  # dplyr::mutate(Group = ifelse(`Characteristics [Category]`=="cell lines" & Group!="blood","Tumor Cell",Group)) %>%
  dplyr::mutate(Role = ifelse(`Characteristics [Cell type]` %in% inhibitory_immune_cell, "Immune cells(Inhibitory)", "Immune cells(Activate)")) %>%
  dplyr::mutate(Role = ifelse(Group=="Tumor Cell","Tumor cells(non-blood)",Role)) %>%
  dplyr::mutate(Role = ifelse(Group=="blood",paste("Tumor cells(blood)",Role,sep="-"),Role)) %>%
  dplyr::mutate(Role = ifelse(Role =="Tumor cells(blood)-Immune cells(Activate)","Immune cells of blood tumor(Activate)",Role)) %>%
  dplyr::mutate(Role = ifelse(Role =="Tumor cells(blood)-Immune cells(Inhibitory)","Immune cells of blood tumor(Inhibitory)",Role)) %>%
  dplyr::inner_join(ICP_exp_site,by="symbol") %>%
  dplyr::mutate(log2gene_tpm = log2(gene_tpm+1) ) %>%
  dplyr::filter(Exp_site!="Not_sure") -> ready_for_analysis

# exp profile between immune and tumor, blood activate and inhibit immne cells---------------------------

ready_for_draw <- ready_for_analysis
ICP_exp_site %>%
  dplyr::inner_join(data.frame(Exp_site = c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor" ),
                               rank = c(5,4,3,2,1)), by = "Exp_site") %>%
  dplyr::filter(Exp_site!="Not_sure") %>%
  dplyr::arrange(rank,`log2FC(I/T)`) %>%
  .$symbol -> symbol_rank

ready_for_draw <- within(ready_for_draw, symbol <- factor(symbol, levels = unique(symbol_rank))) # change symbol's factor level to rank the facets on plot
with(ready_for_draw, levels(symbol))

ready_for_draw %>%
  dplyr::select(symbol,Exp_site) %>%
  # dplyr::filter(symbol %in% c("CD276","CTLA4")) %>%
  unique() -> color_bac
color_bac$Role <- color_bac$log2gene_tpm <- 1

ready_for_draw %>%
  ggplot(aes(x=Role,y=log2gene_tpm)) +
  geom_boxplot() +
  rotate() +
  geom_rect(data=color_bac,aes(fill = Exp_site),xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf,alpha = 0.1) +
  facet_wrap(~symbol,scales = "free_x")  +
  ylab("log2(TPM)") +
  scale_fill_manual(
    # values = site_cplor,
    values = c("yellow",  "green","pink","blue", "red"),
    # values = c("#008B00", "#00EE00", "#CD8500", "#FF4500"),
    breaks = c("Only_exp_on_Immune", "Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor")
  ) +
  theme(
    panel.background = element_rect(fill = "white",colour = "black"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank(),
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
ggsave(file.path(result_path,"pattern_validation","1.FANTOM5.ICP_exp_in_tumor(no-blood-tissue)_Activate-Inhibitory-immune.png"),device = "png",height = 10,width = 20)  
ggsave(file.path(result_path,"pattern_validation","1.FANTOM5.ICP_exp_in_tumor(no-blood-tissue)_Activate-Inhibitory-immune.pdf"),device = "pdf",height = 10,width = 20)  


# tSNE --------------------------------------------------------------------

library("Rtsne")
ready_for_analysis %>%
  dplyr::select(sample,gene_tpm,symbol,Group,Role) %>%
  dplyr::mutate(Role = ifelse(Role %in% c("Immune cells of blood tumor(Activate)","Immune cells of blood tumor(Inhibitory)"),"Blood tumor cells",Role)) %>%
  dplyr::mutate(Role = ifelse(Role %in% c("Tumor cells(non-blood)"),"Solid tumor cells",Role)) %>%
  dplyr::filter(!Role %in% c("Blood tumor cells")) %>%
  tidyr::spread(key="symbol",value="gene_tpm") -> data_for_tSNE

# normalization
data_for_tSNE.matx.normalize <- normalize_input(data_for_tSNE[,-c(1:3)] %>% as.matrix())

# do tSNE
tsne_res <- Rtsne(data_for_tSNE.matx.normalize,dims=2,pca=FALSE,theta=0.0)

# Show the objects in the 2D tsne representation
plot(tsne_res$Y,col=factor(data_for_tSNE$Role), asp=1)
tsne_res$Y %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  # dplyr::rename("tSNE 1"="V1","tSNE 2"="V2") %>%
  dplyr::mutate(sample = data_for_tSNE$sample,
                Cell_type = data_for_tSNE$Role) %>%
  ggplot(aes(x=V1,y= V2)) +
  geom_jitter(aes(color=Cell_type),size=1) +
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  ggpubr::color_palette("jco") +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) 
ggsave(filename = file.path(result_path,"pattern_validation","1.3.FANTOME5.ICP_exp-T-I_tSNE.png"),device = "png",width = 6,height = 3)
ggsave(filename = file.path(result_path,"pattern_validation","1.3.FANTOME5.ICP_exp-T-I_tSNE.pdf"),device = "pdf",width = 6,height = 3)

# PCA analysis of blood tumor immune cells and primary immune cells -------

library("FactoMineR")
library("factoextra")
ready_for_analysis %>%
  dplyr::select(sample,Role,gene_tpm,symbol) %>%
  tidyr::spread(key="symbol",value="gene_tpm") -> ready_for_PCA

ready_for_PCA.df <- as.data.frame(ready_for_PCA[,-c(1)])
rownames(ready_for_PCA.df) <- ready_for_PCA$sample

res.pca <- PCA(ready_for_PCA.df[,-1], scale.unit = T,graph = FALSE)

# As described in previous sections, the eigenvalues measure the amount of variation retained by each principal component. Eigenvalues are large for the first PCs and small for the subsequent PCs. That is, the first PCs corresponds to the directions with the maximum amount of variation in the data set. We examine the eigenvalues to determine the number of principal components to be considered. The eigenvalues and the proportion of variances (i.e., information) retained by the principal components (PCs) can be extracted using the function get_eigenvalue() [factoextra package].
eig.val <- get_eigenvalue(res.pca)
var <- get_pca_var(res.pca)

# PCA biplot
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
fviz_pca_biplot(res.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                labelsize = 2,
                fill.ind = ready_for_PCA.df$Role,
                col.ind = "black",
                # Color variable by groups
                col.var = ICP_exp_site %>%
                  dplyr::inner_join(data.frame(Exp_site = c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor" ),
                                               rank = c(5,4,3,2,1)), by = "Exp_site") %>%
                  dplyr::filter(Exp_site!="Not_sure") %>%
                  dplyr::arrange(rank,`log2FC(I/T)`) %>%
                  .$Exp_site,
                # gradient.cols =c("yellow",  "green","pink","blue", "red"),
                legend.title = list(fill = "Sample class", color = "Exp_site"),
                repel = TRUE        # Avoid label overplotting
)+
  ggpubr::fill_palette("npg")+      # Indiviual fill color
  ggpubr::color_palette("jco")     # Variable colors

ggsave(file.path(result_path,"pattern_validation","2.PCA_analysis_of_ICPs_in_different_cell_types.png"),device = "png",height = 6,width = 8)
ggsave(file.path(result_path,"pattern_validation","2.PCA_analysis_of_ICPs_in_different_cell_types.pdf"),device = "pdf",height = 6,width = 8)


# vertical comparison of ICPs in tumor and immune -------------------------
strip_color <- data.frame(Exp_site = c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor" ),
                          site_cplor = c("blue", "green", "orange", "pink","red"),
                          rank = c(5,4,3,2,1))

my_theme <-   theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_blank(),
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

fn_plot_ICP_exp_in_dataset <- function(.data,ylab,facet,title,filename){
  .data %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mid = quantile(Exp,0.5)) %>%
    dplyr::arrange(mid) %>%
    dplyr::select(symbol,mid) %>%
    unique() %>%
    dplyr::inner_join(gene_list_exp_site,by="symbol") %>%
    dplyr::inner_join(strip_color,by="Exp_site") -> .symbol_rank
  .data %>%
    dplyr::mutate(Exp = log2(Exp+1)) %>%
    ggplot(aes(x=symbol,y=Exp)) +
    geom_boxplot(outlier.colour = "grey",outlier.size = 0.5) +
    facet_wrap(as.formula(facet)) +
    rotate() +
    ggtitle(title) +
    ylab(ylab) +
    xlab("Symbol") +
    scale_x_discrete(limits = .symbol_rank$symbol) +
    my_theme +
    theme(
      axis.text.y = element_text(size = 10,colour = .symbol_rank$site_cplor),
      plot.title = element_text(size=12)
    )
  ggsave(file.path(result_path,"pattern_validation",paste(filename,"pdf",sep=".")),device = "pdf",width = 4, height = 10)
  ggsave(file.path(result_path,"pattern_validation",paste(filename,"png",sep=".")),device = "png", width = 4, height = 10)
}
# immune and tumor cell together
ready_for_analysis %>%
  dplyr::mutate(Exp = gene_tpm) %>%
  dplyr::mutate(cell_source = ifelse(Role %in% c("Immune cells(Activate)","Immune cells(Inhibitory)"),"Immune cell","Tumor cell")) %>%
  dplyr::filter(! Role %in% c("Tumor cells(blood)-Immune cells(Activate)","Tumor cells(blood)-Immune cells(Inhibitory)")) %>%
  fn_plot_ICP_exp_in_dataset(ylab="log2(TPM)",facet="~cell_source",title="FANTOM 5, ICP expression",filename="1.1.FANTOM5.ICP_exp.T-I")

# only tumor cell
ready_for_analysis %>%
  dplyr::mutate(Exp = gene_tpm) %>%
  dplyr::mutate(cell_source = ifelse(Role %in% c("Immune cells(Activate)","Immune cells(Inhibitory)"),"Immune cell","Tumor cell")) %>%
  dplyr::filter(! Role %in% c("Tumor cells(blood)-Immune cells(Activate)","Tumor cells(blood)-Immune cells(Inhibitory)")) %>%
  dplyr::filter(cell_source == "Tumor cell") %>%
  fn_plot_ICP_exp_in_dataset(ylab="log2(TPM)",facet="~cell_source",title="FANTOM 5, ICP expression",filename="1.1.1.FANTOM5.ICP_exp.Tumor")

ready_for_analysis %>%
  dplyr::mutate(Exp = gene_tpm) %>%
  dplyr::mutate(cell_source = ifelse(Role %in% c("Immune cells(Activate)","Immune cells(Inhibitory)"),"Immune cell","Tumor cell")) %>%
  dplyr::filter(! Role %in% c("Tumor cells(blood)-Immune cells(Activate)","Tumor cells(blood)-Immune cells(Inhibitory)")) %>%
  dplyr::filter(cell_source == "Immune cell") %>%
  fn_plot_ICP_exp_in_dataset(ylab="log2(TPM)",facet="~cell_source",title="FANTOM 5, ICP expression",filename="1.1.2.FANTOM5.ICP_exp.Immune")

# correlation of ICP mean exp between immune and tumor --------------------
ready_for_analysis %>%
  dplyr::select(symbol,mean_cell_line,mean_immune_exp,Exp_site) %>%
  unique() %>%
  dplyr::mutate(mean_immune_exp=log2(mean_immune_exp+1),mean_cell_line=log2(mean_cell_line+1))-> ready_for_cor

broom::tidy(
  cor.test(ready_for_cor$mean_immune_exp,ready_for_cor$mean_cell_line,method = "spearman")
) %>%
  dplyr::mutate(y=(max(ready_for_cor$mean_immune_exp)-min(ready_for_cor$mean_immune_exp))*0.75+min(ready_for_cor$mean_immune_exp),x=min(ready_for_cor$mean_cell_line)+(max(ready_for_cor$mean_cell_line)-min(ready_for_cor$mean_cell_line))*0.4) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2))) -> cor_anno

ready_for_cor %>%
  tidyr::nest(-Exp_site) %>%
  dplyr::mutate(spm = purrr::map(data,.f=function(.x){
    if(nrow(.x)>5){
      broom::tidy(
        cor.test(.x$mean_immune_exp,.x$mean_cell_line,method = "spearman")) %>%
        dplyr::mutate(y=(max(.x$mean_immune_exp)-min(.x$mean_immune_exp))*0.85+min(.x$mean_immune_exp),x=cor_anno$x) %>%
        dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2)))
    }else{
      tibble::tibble()
    }
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> cor_anno.specific

ready_for_cor %>%
  dplyr::filter(Exp_site!="Not_sure") %>%
  ggplot(aes(x=mean_cell_line,y=mean_immune_exp)) +
  geom_jitter(aes(color=Exp_site),size=0.5) +
  geom_smooth(se= F, method = "lm", aes(color=Exp_site,group=Exp_site)) +
  geom_smooth(se = F, method = "lm",color= "black") +
  geom_text(aes(x=x,y=y,label = label),data=cor_anno,color="black") +
  geom_text(aes(x=x,y=y,label = label,color=Exp_site),data = cor_anno.specific) +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme +
  xlab("log2(Mean exppression in tumor cells)") +
  ylab("log2(Mean exppression in immune cells)") 
ggsave(file.path(result_path,"pattern_validation","1.2.FANTOM5-T-I-meanExp.correlation.pdf"),device = "pdf",height = 6,width = 8)
ggsave(file.path(result_path,"pattern_validation","1.2.FANTOM5-T-I-meanExp.correlation.png"),device = "png",height = 6,width = 8)

# correlation of ICP rank between immune and tumor --------------------
ready_for_cor %>%
  dplyr::arrange(mean_cell_line) %>%
  dplyr::mutate(Tumor_rank=rank(mean_cell_line)) %>%
  dplyr::arrange(mean_immune_exp) %>%
  dplyr::mutate(Immune_rank=1:nrow(ready_for_cor)) -> ready_for_rank.cor

broom::tidy(
  cor.test(ready_for_rank.cor$Immune_rank,ready_for_rank.cor$Tumor_rank,method = "spearman")
) %>%
  dplyr::mutate(y=(max(ready_for_rank.cor$Immune_rank)-min(ready_for_rank.cor$Immune_rank))*0.75+min(ready_for_rank.cor$Immune_rank),x=min(ready_for_rank.cor$Tumor_rank)+(max(ready_for_rank.cor$Tumor_rank)-min(ready_for_rank.cor$Tumor_rank))*0.4) %>%
  dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2))) -> cor_anno

ready_for_rank.cor %>%
  tidyr::nest(-Exp_site) %>%
  dplyr::mutate(spm = purrr::map(data,.f=function(.x){
    if(nrow(.x)>5){
      broom::tidy(
        cor.test(.x$Immune_rank,.x$Tumor_rank,method = "spearman")) %>%
        dplyr::mutate(y=(max(.x$Immune_rank)-min(.x$Immune_rank))*0.85+min(.x$Immune_rank),x=cor_anno$x) %>%
        dplyr::mutate(label = paste("r = ",signif(estimate,2),", p = ",signif(p.value,2)))
    }else{
      tibble::tibble()
    }
  })) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> cor_anno.specific

ready_for_rank.cor %>%
  dplyr::filter(Exp_site!="Not_sure") %>%
  ggplot(aes(x=Tumor_rank,y=Immune_rank)) +
  geom_jitter(aes(color=Exp_site),size=0.5) +
  geom_smooth(se= F, method = "lm", aes(color=Exp_site,group=Exp_site)) +
  geom_smooth(se = F, method = "lm",color= "black") +
  geom_text(aes(x=x,y=y,label = label),data=cor_anno,color="black") +
  geom_text(aes(x=x,y=y,label = label,color=Exp_site),data = cor_anno.specific) +
  scale_color_manual(values=c("#CD661D",  "#008B00", "#FF69B4", "#1874CD","#CD3333")) +
  my_theme +
  xlab("log2(Mean exppression in tumor cells)") +
  ylab("log2(Mean exppression in immune cells)") 
ggsave(file.path(result_path,"pattern_validation","1.2.FANTOM5-T-I-meanExp.correlation.pdf"),device = "pdf",height = 6,width = 8)
ggsave(file.path(result_path,"pattern_validation","1.2.FANTOM5-T-I-meanExp.correlation.png"),device = "png",height = 6,width = 8)
# save image --------------------------------------------------------------

save.image(file.path(result_path,"pattern_validation","FANTOM5.validation.Rdata"))
