############## verify the gene expression site in multiple data set ############
library(magrittr)

basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio")

load(file.path(res_path,"pattern_validation","ICP_exp.Rdata"))

#### gene list ----
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_exp_site <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T)`) %>%
  dplyr::inner_join(gene_list,by="symbol")

#### Zhang ZM Single cell deep sequence of T cells (n=5063)------
liver_Tcell_ICP_exp <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/liver_single_cell_TPM-ZhangZM","GSE98638_HCC.TCell.S5063.TPM.txt")) %>%
  dplyr::filter(geneID %in% gene_list_exp_site$entrez_ID)
dim(liver_Tcell_ICP_exp) # 68 rows

# sample info
liver_Tcell_sample <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/liver_single_cell_TPM-ZhangZM","patient_info.txt"),col_names = F) %>%
  dplyr::select("X1","X3","X4") %>%
  dplyr::rename("sample"="X1","cell_cluster"="X3","source_name"="X4") %>%
  dplyr::mutate(sample = gsub("\\.","-",sample))
liver_Tcell_ICP_exp %>%
  tidyr::gather(-geneID,-symbol,key="sample",value="Exp")%>%
  dplyr::inner_join(liver_Tcell_sample,by="sample") -> liver_Tcell_ICP_exp.gather

#### GSE22886, array Expression profiles from a variety of resting and activated human immune cells -----
# gene info
GSE22886_gene <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/GSE22886","GPL96_gene_anno.txt")) %>%
  dplyr::select(ID,`Gene Symbol`,`ENTREZ_GENE_ID`)

# sample info
GSE22886_sample_info <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/GSE22886","sample_anno.txt"),col_names = F) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::select(V1,V2,V8) %>%
  dplyr::rename("title"="V1","sample"="V2","source_name"="V8") %>%
  .[-1,]

# load exp data
GSE22886_immunecell_ICP_exp <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/GSE22886","GSE22886-GPL96_series_matrix_RNA-expression-filter.txt")) %>%
  dplyr::rename("ID"="ID_REF") %>%
  dplyr::inner_join(GSE22886_gene,by="ID") %>%
  dplyr::filter(ENTREZ_GENE_ID %in% gene_list_exp_site$entrez_ID)

GSE22886_immunecell_ICP_exp$ENTREZ_GENE_ID %>% unique() %>% length() #57 genes 

# select transcript with biggest expression, and add sample info
GSE22886_immunecell_ICP_exp %>%
  dplyr::select(-ID) %>%
  tidyr::gather(-ENTREZ_GENE_ID,-`Gene Symbol`,key="sample",value="Exp") %>%
  dplyr::group_by(ENTREZ_GENE_ID,sample) %>%
  dplyr::mutate(Exp_max = max(Exp)) %>%
  dplyr::select(-Exp) %>%
  unique() %>%
  dplyr::rename("Exp" = "Exp_max") %>%
  dplyr::inner_join(GSE22886_sample_info,by="sample") %>%
  dplyr::ungroup()-> GSE22886_immunecell_ICP_exp.filter

#### GSE49910, An Expression Atlas of Human Primary Cells: Inference of Gene Function from Coexpression Networks -----
# gene info
GSE49910_gene <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/GSE49910","GPL570_gene_anno.txt")) %>%
  dplyr::select(ID,`Gene Symbol`,`ENTREZ_GENE_ID`)
# sample info
GSE49910_sample_info <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/GSE49910","sample_anno.txt"),col_names = F) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::select(V1,V2,V8) %>%
  dplyr::rename("title"="V1","sample"="V2","source_name"="V8") %>%
  .[-1,]

# load exp data
GSE49910_immunecell_ICP_exp <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/GSE49910","GSE49910_series_matrix_RNA-expression-filter.txt")) %>%
  dplyr::rename("ID"="ID_REF") %>%
  dplyr::inner_join(GSE22886_gene,by="ID") %>%
  dplyr::filter(ENTREZ_GENE_ID %in% gene_list_exp_site$entrez_ID)

GSE49910_immunecell_ICP_exp$ENTREZ_GENE_ID %>% unique() %>% length() #57 genes 

# select transcript with biggest expression, and add sample info
GSE49910_immunecell_ICP_exp %>%
  dplyr::select(-ID) %>%
  tidyr::gather(-ENTREZ_GENE_ID,-`Gene Symbol`,key="sample",value="Exp") %>%
  dplyr::group_by(ENTREZ_GENE_ID,sample) %>%
  dplyr::mutate(Exp_max = max(Exp)) %>%
  dplyr::select(-Exp) %>%
  unique() %>%
  dplyr::rename("Exp" = "Exp_max") %>%
  dplyr::inner_join(GSE49910_sample_info,by="sample") %>%
  dplyr::ungroup()-> GSE49910_immunecell_ICP_exp.filter

#### CCLE cancer cell lines expression data -----
# gene info
ICP_ENSG <- readr::read_tsv(file.path(basic_path,"data/ccle_database/Expression_Atlas_release_31","Homo_sapiens.GRCh38.96-filtered.gtf"),col_names = F) %>%
  dplyr::filter(X2 %in% gene_list_exp_site$symbol) %>%
  dplyr::rename("ENSG" = "X1", "symbol" = "X2")
gene_list_exp_site %>%
  dplyr::filter(!symbol %in% ICP_ENSG$symbol) 

readr::read_tsv(file.path(basic_path,"data/ccle_database/Expression_Atlas_release_31","Homo_sapiens.GRCh38.96-filtered.gtf"),col_names = F) %>%
  dplyr::filter(X1 %in% c("ENSG00000227357","ENSG00000107738", "ENSG00000278152","ENSG00000277885","ENSG00000275452","ENSG00000276258","ENSG00000275737","ENSG00000276139","ENSG00000274438","ENSG00000275253","ENSG00000275735","ENSG00000274518","ENSG00000276425","ENSG00000275583","ENSG00000278300","ENSG00000277725","ENSG00000276011","ENSG00000275407","ENSG00000276731","ENSG00000275546","ENSG00000275914","ENSG00000277251","ENSG00000278692","ENSG00000273661","ENSG00000278731","ENSG00000273578","ENSG00000274412","ENSG00000275960")) %>% 
  dplyr::rename("ENSG" = "X1", "symbol" = "X2") %>%
  rbind(ICP_ENSG)%>%
  dplyr::mutate(symbol= ifelse( symbol == "VSIR","C10orf54",symbol)) -> ICP_ENSG

# expression data
ICP_ccle_exp <- readr::read_tsv(file.path(basic_path,"data/ccle_database/Expression_Atlas_release_31","tpms.filtered.tsv")) %>%
  dplyr::filter(`Gene ID` %in% ICP_ENSG$ENSG )

# sample info
sample_info <- readr::read_tsv(file.path(basic_path,"data/ccle_database/Expression_Atlas_release_31","experiment-design")) %>%
  dplyr::mutate(sample = paste(`Sample Characteristic[disease]`,`Sample Characteristic[cell line]`,sep=", ")) %>%
  dplyr::select(sample,source_name = `Sample Characteristic[disease]`)

ICP_ccle_exp %>%
  dplyr::rename("ENSG" = "Gene ID", "symbol" = "Gene Name") %>%
  tidyr::gather(-ENSG,-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(sample_info,by="sample") -> ICP_ccle_exp_data

#### GTEx normal tissue expression data -----
GTEX_sym <- readr::read_rds(file.path(basic_path,"data/id_corresponding/id_correspond_between_NCBI_TCGA.rds.gz")) %>%
  dplyr::filter(GeneID %in% gene_list_exp_site$GeneID) %>%
  .$GTEX_sym %>% unique()

GTEx_expr <- readr::read_tsv(file.path(basic_path,"data/GSCALite/GTEx/expression","GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")) %>%
  dplyr::filter(Description %in% GTEX_sym)

Gtex.sample <- readr::read_tsv(file.path(basic_path,"data/GTEx","GTEx_v7_Annotations_SampleAttributesDS.txt")) %>%
  dplyr::select(SAMPID,SMTS) %>%
  dplyr::rename("sample"="SAMPID","source_name"="SMTS")

GTEx_expr %>%
  dplyr::rename("ENSG" = "Name", "symbol" = "Description") %>%
  tidyr::gather(-ENSG,-symbol,key="sample",value="Exp") %>%
  dplyr::inner_join(Gtex.sample,by="sample") -> ICP_GTEx_expr


# validation function -----------------------------------------------------
ICP_GTEx_expr -> .data
library(ggplot2)
library(ggpubr)

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
strip_color <- data.frame(Exp_site = c("Only_exp_on_Immune","Mainly_exp_on_Immune","Both_exp_on_Tumor_Immune","Mainly_exp_on_Tumor","Only_exp_on_Tumor" ),
                          site_cplor = c("blue", "green", "yellow", "pink","red"),
                          rank = c(5,4,3,2,1))
fn_plot_ICP_exp_in_dataset <- function(.data,ylab,title,filename){
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
    rotate() +
    # ylab(ylab) +
    xlab("Symbol") +
    # title(main = title) +
    scale_x_discrete(limits = .symbol_rank$symbol) +
    my_theme +
    theme(
      axis.text.y = element_text(size = 10,colour = .symbol_rank$site_cplor)
    )
  ggsave(file.path(res_path,"pattern_validation",paste(filename,"pdf",sep=".")),device = "pdf",width = 4, height = 10)
  ggsave(file.path(res_path,"pattern_validation",paste(filename,"png",sep=".")),device = "png", width = 4, height = 10)
  
  # output the rank
  nr <-  nrow(.symbol_rank)
  .symbol_rank %>%
    dplyr::ungroup() %>%
    dplyr::mutate(n = 1:nr) %>%
    dplyr::select(symbol,n)
}


# validation --------------------------------------------------------------
fn_plot_ICP_exp_in_dataset(.data=GSE22886_immunecell_ICP_exp.filter %>%
                             dplyr::rename("symbol"="Gene Symbol"),
                           ylab="log2(Array exp)",
                           title = "GSE22886, array, 20 immune cells",
                           filename = "3.1.GSE22886.20immunecells.overall") -> GSE22886.ICP_rank

fn_plot_ICP_exp_in_dataset(.data=GSE49910_immunecell_ICP_exp.filter %>%
                             dplyr::rename("symbol"="Gene Symbol"),
                           ylab="log2(Array exp)",
                           title = "GSE49910, array, 11 immune cells, 46 samples",
                           filename = "3.2.GSE49910.11immunecells.overall") -> GSE49910.ICP_rank

fn_plot_ICP_exp_in_dataset(.data=liver_Tcell_ICP_exp.gather,
                           ylab="log2(TPM)",
                           title = "cell.2017.05.035, single cell sequence, T cells from liver",
                           filename = "3.3.cell.2017.05.035.Tcells.overall") -> liverTcell.ICP_rank

fn_plot_ICP_exp_in_dataset(.data=ICP_ccle_exp_data %>%
                             dplyr::filter(!is.na(Exp)),
                           ylab="log2(TPM)",
                           title = "CCLE, RNAseq, tumor cell lines",
                           filename = "3.4.ccle.tumorcells.overall") -> ccle.ICP_rank

fn_plot_ICP_exp_in_dataset(.data=ICP_GTEx_expr,
                           ylab="log2(Array exp)",
                           title = "GTEx, array, 20 immune cells",
                           filename = "3.5.GTEx.NormalTissue.overall") -> GTEx.ICP_rank



save.image(file.path(res_path,"pattern_validation","ICP_exp.Rdata"))
