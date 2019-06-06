## Data source : FANTOM 5, Update of the FANTOM web resource: high resolution transcriptome of diverse cell types in mammals
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5210666/
## /project/huff/huff/data/FANTOM5/extra/matrix.hg38_fair_new_CAGE_peaks_phase1and2_tpm_ann.osc.txt
# FANTOM 5 data filter and classification ---------------------------------

library(magrittr)
# library(org.Hs.eg.db)
library(clusterProfiler)
# data path ---------------------------------------------------------------
basic_path <- "/home/huff/project"
fantom_path <- file.path(basic_path,"data/FANTOM5/extra")
immune_path <- file.path(basic_path,"immune_checkpoint")
gene_list_path <-file.path(immune_path,"checkpoint/20171021_checkpoint")


# load data ---------------------------------------------------------------

fantom_exp <- readr::read_tsv(file.path(fantom_path,"matrix.hg38_fair_new_CAGE_peaks_phase1and2_tpm_ann.osc.txt"))
fantom_sample_info <- readr::read_tsv(file.path(fantom_path,"HumanSamples2.0.classification.txt"))

# dealing with peaks such as p1@KIR2DL2,p1@KIR2DL3,p1@KIR2DS2 covered several genes ------
fantom_exp %>%
  dplyr::select(short_description,entrezgene_id,hgnc_id) %>%
  tidyr::separate(entrezgene_id,paste("entrez",1:10,sep=" "),"/")
# load gene list and correspond TCGA symbol with FANTOM symbol id ----
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)

ICP_HGNC_symbol <- readr::read_tsv(file.path(basic_path,"data/HGNC","All_approved_symbols.txt")) %>%
  dplyr::select(hgnc_id,entrez_id,symbol) %>%
  dplyr::filter(entrez_id %in% gene_list$GeneID)

xCell_cell_type <- readr::read_tsv(file.path(xcell_path,"cell_type_classification.13059_2017_1349_MOESM1_ESM.txt")) %>%
  dplyr::select(`Full name`,`Cell types`,Subgroup) 

xCell_Fantom_cellname_adjust <- readr::read_tsv(file.path(xcell_path,"xCell_Fantom5_cellname_uniform.txt")) %>%
  dplyr::rename("Cell types" = "xCell_Cell types") %>%
  dplyr::inner_join(xCell_cell_type,by="Cell types") %>%
  dplyr::mutate(Type = ifelse(Subgroup == "Lymphoid" | Subgroup == "Myeloid","Immune Cell","Stromal Cell")) %>%
  dplyr::select(Fantom5_cellname,Type) %>%
  unique()
# data filter -------------------------------------------------------------
##### filter ICP genes in expression data, and catch the key info of samples ----
# filter ICP genes
fantom_exp %>%
  dplyr::filter(hgnc_id %in% ICP_HGNC_symbol$hgnc_id) %>%
  dplyr::select(-`00Annotation`,-description,-association_with_transcript,-entrezgene_id,-uniprot_id) %>%
  tidyr::gather(-hgnc_id,-short_description,key="sample",value="TPM") -> ICP_fantom_exp

# catch the key info of samples 
fn_substr <- function(.x){
  .tmp <- as.character(t(stringr::str_split(.x,pattern = "\\.",n = 5)[[1]][2])) 
  .tmp_1 <- as.character(t(stringr::str_split(.x,pattern = "\\.",n = 5)[[1]][3]))
  if(substr(.tmp_1,1,4)!="CNhs"){
    .tmp_1 <-.tmp_1
  }else{
    .tmp_1 <-as.character(t(stringr::str_split(.x,pattern = "\\.",n = 5)[[1]][4]))
  }
  return(toupper(paste(.tmp,.tmp_1,sep=".")))
}
ICP_fantom_exp %>%
  dplyr::mutate(sample = purrr::map(sample,fn_substr)) %>%
  dplyr::mutate(sample = as.character(sample))-> ICP_fantom_exp.substr;ICP_fantom_exp.substr

ICP_fantom_exp.substr %>%
  dplyr::group_by(hgnc_id,sample) %>%
  dplyr::mutate(gene_tpm = sum(TPM)) %>%
  dplyr::select(hgnc_id,sample,gene_tpm) %>% 
  unique() %>%
  dplyr::ungroup() -> ICP_fantom_gene.exp.substr

##### filter cell line and primary cell samples ----
# filter samples we need and transcode the sample names by curl_escape() of curl R package
library(curl)
fantom_sample_info %>%
  dplyr::filter(`Characteristics [Category]` %in% c("cell lines","primary cells")) %>%
  dplyr::mutate(sample = toupper(paste(curl_escape(`Charateristics [description]`),`Source Name`,sep="."))) %>%
  dplyr::select(sample,`Characteristics[Tissue]`,`Characteristics [Cell type]`,`Characteristics [Category]`) -> fantom_sample_info_transcode


# merge expression data and sample info by sample keys --------------------

# cell lines = tumor expression ----
fantom_sample_info_transcode %>%
  dplyr::filter(`Characteristics [Category]` == "cell lines") %>%
  dplyr::inner_join(ICP_fantom_gene.exp.substr,by="sample") %>%
  dplyr::group_by(hgnc_id,`Characteristics[Tissue]`) %>%
  dplyr::mutate(gene_mean_exp = mean(gene_tpm)) %>%
  dplyr::select(hgnc_id,`Characteristics[Tissue]`,gene_mean_exp) %>%
  unique() %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(ICP_HGNC_symbol,by="hgnc_id")  %>%
  dplyr::mutate(Group = "Tumor Cell")-> ICP_fantom_gene.exp.cell_line
  
# primary cell = immune cell and stromal cell ----
fantom_sample_info_transcode %>%
  dplyr::filter(`Characteristics [Category]` == "primary cells") %>%
  dplyr::rename("Fantom5_cellname"="Characteristics [Cell type]") %>%
  dplyr::inner_join(xCell_Fantom_cellname_adjust,by="Fantom5_cellname") %>%
  # dplyr::filter(! is.na(`Characteristics[Tissue]`)) %>%
  dplyr::rename("Group" = "Type") %>%
  dplyr::inner_join(ICP_fantom_gene.exp.substr,by="sample") %>%
  dplyr::group_by(hgnc_id,Group) %>%
  dplyr::mutate(gene_mean_exp = mean(gene_tpm)) %>%
  dplyr::select(hgnc_id,Group,gene_mean_exp) %>%
  unique() %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(ICP_HGNC_symbol,by="hgnc_id") -> ICP_fantom_gene.exp.Immune_Stromal_cell


# data combination --------------------------------------------------------
ICP_fantom_gene.exp.Immune_Stromal_cell %>%
  dplyr::mutate(`Characteristics[Tissue]`="PrimaryCell") %>%
  rbind(ICP_fantom_gene.exp.cell_line) -> ICP_fantom.gene_exp.cell_line.Immune_cell.combine

ICP_fantom.gene_exp.cell_line.Immune_cell.combine %>%
  readr::write_rds(file.path(immune_path,"genelist_data","ICP_fantom.gene_exp.cell_line.Immune_cell.combine.rds.gz"),compress = "gz")

# all_cancer_TIL <- readr::read_rds("/project/huff/huff/data/TCGA/immune_infiltration/miao_TCAP_prediction_for_all_samples/All_TCGA_sample_TIL.rds.gz")
