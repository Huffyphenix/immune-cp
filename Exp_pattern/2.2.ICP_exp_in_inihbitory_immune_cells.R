############## ICP exp in inhibitory immune cell ###############
#################### FANTOM data ###############################
library(curl)
library(ggbeeswarm)

# data path ---------------------------------------------------------------
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
result_path <- file.path(immune_path,"result_20171025","ICP_exp_patthern-byratio")
fantom_path <- file.path(basic_path,"data/FANTOM5/extra")

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


ICP_exp_site <- readr::read_tsv(file.path(result_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv"))

inhibitory_immune_cell<- c("monocyte","neutrophil","monoblast")
ICP_fantom.gene_exp.Immune_cell.combine %>%
  dplyr::inner_join(fantom_sample_info,by="sample") %>%
  dplyr::mutate(`Characteristics[Tissue].x` = ifelse(is.na(`Characteristics[Tissue].x`),"notsure",`Characteristics[Tissue].x`)) %>%
  dplyr::mutate(Group = ifelse(`Characteristics[Tissue].x`=="blood","blood",Group)) %>%
  # dplyr::mutate(Group = ifelse(`Characteristics [Category]`=="cell lines" & Group!="blood","Tumor Cell",Group)) %>%
  dplyr::mutate(Role = ifelse(`Characteristics [Cell type]` %in% inhibitory_immune_cell, "Immune cells(Inhibitory)", "Immune cells(Activate)")) %>%
  dplyr::mutate(Role = ifelse(Group=="Tumor Cell","Tumor cells(non-blood)",Role)) %>%
  dplyr::mutate(Role = ifelse(Group=="blood",paste("Tumor cells(blood)",Role,sep="-"),Role)) %>%
  dplyr::inner_join(ICP_exp_site,by="symbol") %>%
  dplyr::mutate(gene_tpm = log2(gene_tpm+1) ) %>%
  dplyr::filter(Exp_site!="Not_sure") -> ready_for_draw


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
color_bac$Role <- color_bac$gene_tpm <- 1

ready_for_draw %>%
  ggplot(aes(x=Role,y=gene_tpm)) +
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


# PCA analysis of blood tumor immune cells and primary immune cells -------

library("FactoMineR")
library("factoextra")
ready_for_draw %>%
  # dplyr::filter(Group %in% c("Immune Cell","blood")) %>%
  dplyr::select(sample,Role,gene_tpm,symbol) %>%
  dplyr::mutate(gene_tpm=2^gene_tpm) %>%
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
  ggpubr::fill_palette("jco")+      # Indiviual fill color
  ggpubr::color_palette("npg")     # Variable colors

ggsave(file.path(result_path,"pattern_validation","2.PCA_analysis_of_ICPs_in_different_cell_types.png"),device = "png",height = 8,width = 10)
ggsave(file.path(result_path,"pattern_validation","2.PCA_analysis_of_ICPs_in_different_cell_types.pdf"),device = "pdf",height = 8,width = 10)

