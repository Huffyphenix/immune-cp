library(magrittr)
imm_path<-c("/project/huff/huff/immune_checkpoint/data/10kImmunomes.All.2017-09-20/allDataForApp/Normalized Files")
cytof_pbmc<-read.csv(file.path(imm_path,"cytof_pbmc_normalized.csv",fsep=.Platform$file.sep))
gene_pbmc<-read.csv(file.path(imm_path,"gene_pbmc_normalized.csv",fsep=.Platform$file.sep))
out_path<-"/project/huff/huff/immune_checkpoint/data"
gene_list_path <- c("/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint")
gene_list <- read.table(file.path(gene_list_path, "all.entrez_id-gene_id"),header=T)
gene_type<-read.table(file.path(gene_list_path,"checkpoint.type"),header=T)
gene_list<-dplyr::left_join(gene_list,gene_type,by="symbol")

#absolute_immunity<-readr::read_rds("S:/study/生存分析/免疫检查点project/liucj_tcga_process_data/immunity_raw.rds.gz")

###change data form, get import data.
gene_pbmc %>%
  dplyr::select(-study_accession,-age,-race,-gender,-data_accession,-class,-sample_type) %>%
  tidyr::gather(symbol,expr,-subject_accession) %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(class=ifelse(expr>=mean(expr),"up","down")) ->gene_pbmc.1
cytof_pbmc %>%
  dplyr::select(-study_accession,-age,-race,-gender) %>%
  tidyr::gather(Immune_cell,count,-subject_accession) ->cytof_pbmc.1
gene_pbmc.1 %>%
  dplyr::inner_join(cytof_pbmc.1,by="subject_accession") ->gene_cellCount_pbmc_10kImmun.match
readr::write_rds(gene_cellCount_pbmc_10kImmun.match,
                 path = file.path(out_path,"10kImmun_gene_cellCount_pbmc.match.rds.gz"),compress = "gz")
