############## ICP exp in inhibitory immune cell ###############
#################### FANTOM data ###############################
library(curl)
library(ggbeeswarm)

# data path ---------------------------------------------------------------
basic_path <- "/home/huff/project"
immune_path <- file.path(basic_path,"immune_checkpoint")
result_path <- file.path(immune_path,"result_20171025","ICP_exp_patthern")
fantom_path <- file.path(basic_path,"data/FANTOM5/extra")

# load data ---------------------------------------------------------------

ICP_fantom.gene_exp.Immune_cell.combine <-
  readr::read_rds(file.path(immune_path,"genelist_data","FANTOM5","ICP_fantom.gene_exp.cell_line.Immune_cell.raw.exp.rds.gz")) %>%
  dplyr::filter(Group == "Immune Cell")

fantom_sample_info <- readr::read_tsv(file.path(fantom_path,"HumanSamples2.0.classification.txt")) %>%
  dplyr::filter(`Characteristics [Category]` %in% c("primary cells")) %>%
  dplyr::mutate(sample = toupper(paste(curl_escape(`Charateristics [description]`),`Source Name`,sep="."))) %>%
  dplyr::select(sample,`Charateristics [description]`,`Characteristics[Tissue]`,`Characteristics [Cell type]`,`Characteristics [Category]`)  %>%
  dplyr::filter(sample %in% unique(ICP_fantom.gene_exp.Immune_cell.combine$sample))


# ICP in each cell types ----------------------------------------------
inhibitory_immune_cell<- c("monocyte","neutrophil")
ICP_fantom.gene_exp.Immune_cell.combine %>%
  dplyr::inner_join(fantom_sample_info,by="sample") %>%
  dplyr::mutate(inhibitory = ifelse(`Characteristics [Cell type]` %in% inhibitory_immune_cell, "yes", "no")) -> ICP_fantom.gene_exp.Immune_cell.combine.cell_type

ICP_fantom.gene_exp.Immune_cell.combine.cell_type %>%
  dplyr::filter(inhibitory == "yes") %>%
  dplyr::mutate(gene_tpm = log2(gene_tpm+1) ) %>%
  ggplot(aes(x=symbol,y=gene_tpm)) +
  geom_boxplot() +
  facet_wrap(~`Characteristics [Cell type]`,scales = "free") +
  rotate()
