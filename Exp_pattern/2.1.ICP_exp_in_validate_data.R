############## verify the gene expression site in multiple data set ############
library(magrittr)

basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")
res_path <- file.path(immune_path,"result_20171025/ICP_exp_patthern-byratio")

#### gene list ----
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_list_exp_site <- readr::read_tsv(file.path(res_path,"pattern_info","ICP_exp_pattern_in_immune_tumor_cell-by-FC-pvalue.tsv")) %>%
  dplyr::select(entrez_ID,symbol,Exp_site,`log2FC(I/T)`) %>%
  dplyr::inner_join(gene_list,by="symbol")

#### Zhang ZM Single cell deep sequence of T cells (n=5063)------
liver_Tcell_ICP_exp <- readr::read_tsv(file.path(basic_path,"data/immune_cell_exp_data/liver_single_cell_TPM-ZhangZM","GSE98638_HCC.TCell.S5063.TPM.txt")) %>%
  dplyr::filter(geneID %in% gene_list_exp_site$entrez_ID)
dim(liver_Tcell_ICP_exp) # 68 rows


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
ccle_exp <- readr::read_tsv(file.path(basic_path,"data/ccle_database","CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct"))

#### GTEx normal tissue expression data -----
GTEX_sym <- readr::read_rds(file.path(basic_path,"data/id_corresponding/id_correspond_between_NCBI_TCGA.rds.gz")) %>%
  dplyr::filter(GeneID %in% gene_list_exp_site$GeneID) %>%
  .$GTEX_sym %>% unique()
GTEx_expr <- readr::read_tsv(file.path(basic_path,"data/GSCALite/GTEx/expression","GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct")) %>%
  dplyr::filter(Description %in% GTEX_sym)
