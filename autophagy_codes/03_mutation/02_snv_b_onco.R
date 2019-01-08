library(magrittr)
out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
snv_path <- file.path(out_path, "m_2_snv")
tcga_path <- "/project/huff/huff/immune_checkpoint/data/TCGA_data"
expr_path <-c("/project/huff/huff/immune_checkpoint/result_20171025/expr_rds")

# load mc3
mc3_file <- file.path(tcga_path, 'syn_mutation_syn7824274_mc3_public.pass.maf')

suppressPackageStartupMessages(require(maftools))
gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)
# mc3_maf <- read.maf(mc3_file)
# readr::write_rds(mc3_maf, path = file.path(tcga_path, "syn_mutation_syn7824274_mc3_public.pass.maf.rds.gz"), compress = 'gz')
readr::read_rds(path = file.path(tcga_path, "syn_mutation_syn7824274_mc3_public.pass.maf.rds.gz")) -> mc3_maf

subsetMaf(mc3_maf, genes = gene_list %>% dplyr::pull(symbol), mafObj = T) -> gene_list_maf

plotmafSummary(gene_list_maf)

oncodrive(maf = gene_list_maf, AACol = "HGVSp", minMut = 5, pvalMethod = 'zscore') -> gene_list_maf_sig
plotOncodrive(res = gene_list_maf_sig, fdrCutOff = 0.5, useFraction = T)
pfamDomains(maf = gene_list_maf, top = 10)
oncoplot(maf = gene_list_maf, removeNonMutated = T)

# getFields(mc3_maf)
# getGeneSummary(mc3_maf)
pdf(file.path(snv_path, "genelist_oncoplot.pdf"), width=10, height=10)
oncoplot(maf = mc3_maf, genes = gene_list %>% dplyr::pull(symbol), removeNonMutated = TRUE)
dev.off()

# pdf(file.path(snv_path, "lys_oncoplot.pdf"), width=10, height=10)
# oncoplot(maf = mc3_maf, genes = gene_list %>% dplyr::filter(status == "l") %>% dplyr::pull(symbol), removeNonMutated = TRUE)
# dev.off()

# lollipop 
# fn_loli <- function(.x){
#   tryCatch({
#   path <- file.path(snv_path, 'lollipop')
#   pdf(file.path(path, paste(.x, "lollipop.pdf", sep = ".")), width=12, height=3)
#   lpop = lollipopPlot(maf = mc3_maf, gene = .x,  showMutationRate = F, defaultYaxis = FALSE)
#   dev.off()},
#   error = function(e){1},
#   warning = function(e){1})
#   
# }
# gene_list %>% 
#   dplyr::filter(pathway == "autophagesome formation-core") %>% 
#   dplyr::pull(symbol) %>% 
#   purrr::walk(.f = fn_loli)




save.image(file = file.path(snv_path, ".rda_02_snv_b_onco.rda"))
load(file = file.path(snv_path, ".rda_02_snv_b_onco.rda"))


