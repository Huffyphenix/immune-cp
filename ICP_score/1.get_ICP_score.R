library()

# load data ---------------------------------------------------------------

data_result_path <- "/project/huff/huff/immune_checkpoint/genelist_data"

## get from 1.data_prepare.R
load(file = file.path(data_result_path, ".rda_IMK_mutationburden_cancerSubtype_analysis.rda"))

# reduce the dataset to the top 5,000 most variable genes, measured by median absolute deviation(mad). 绝对中位差
# expression data
expr_mads=apply(gene_list_expr_with_mutation_load.matrix.combine,1,mad)
expr_top20_mad <- names(rev(sort(expr_mads))[1:30]) # 排序，获得rank值，倒置，取前20
 
# methylation data
methy_mads=apply(PanCan26_gene_list_methy_matrix.combine,1,mad)
methy_top20_mad <- names(rev(sort(methy_mads))[1:30]) # 排序，获得rank值，倒置，取前20

# CNV data
cnv_mads=apply(cnv_merge_snv_data.matrix.combine,1,mad)
cnv_top20_mad <- names(rev(sort(cnv_mads))[1:30]) # 排序，获得rank值，倒置，取前20

# genes as the top 30 variable in at least two data type
intersect(sort(expr_top20_mad[1:10]),sort(methy_top20_mad[1:10])) -> expr_methy_top
intersect(sort(expr_top20_mad[1:10]),sort(cnv_top20_mad[1:10])) -> expr_cnv_top
intersect(sort(cnv_top20_mad[1:10]),sort(methy_top20_mad[1:10])) -> methy_cnv_top

c(expr_methy_top,expr_cnv_top,methy_cnv_top) %>% unique() -> at_leat_two # 30 genes

# genes as the top 30 variable in both three data type
intersect(sort(expr_top20_mad),sort(methy_top20_mad)) %>% intersect(sort(cnv_top20_mad))-> expr_methy_cnv_top # "CD40"   "LGALS9" "SIRPA" 
