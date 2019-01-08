library(magrittr)

out_path <- "/project/huff/huff/immune_checkpoint/result_20171025"
stage_path <- file.path(out_path, "c_1_stage")
subtype_path <- file.path(out_path, "c_2_subtype")
survival_path <- file.path(out_path, "c_3_survival")
clinical_path<-file.path(out_path,"clinical_result")
#read in results
stage<-readr::read_csv(file.path(stage_path,"03_b_stage_gene_fdr0.05.csv"))
subtype<-readr::read_csv(file.path(subtype_path,"03_c_subtype_gene_fdr0.05.csv"))
survival<-readr::read_csv(file.path(survival_path,"c_3_survival_genelist_sig_pval.csv"))

#gene list
gene_list_path <- "/project/huff/huff/immune_checkpoint/checkpoint/20171021_checkpoint"
gene_list <- read.table(file.path(gene_list_path, "gene_list_type"),header=T)
gene_list$symbol<-as.character(gene_list$symbol)

#merge data
stage %>%
  dplyr::mutate(n=1) %>%
  dplyr::mutate(cancer=paste(cancer_types,"st",sep=".")) %>%
  dplyr::select(cancer,symbol,n) %>%
  tidyr::spread(key=cancer,value=n) ->stage.st

subtype %>%
  dplyr::mutate(n=1) %>%
  dplyr::mutate(cancer=paste(cancer_types,"sub",sep=".")) %>%
  dplyr::select(cancer,symbol,n) %>%
  tidyr::spread(key=cancer,value=n) ->stage.sub

survival %>%
  dplyr::mutate(cancer=paste(cancer_types,"sur",sep=".")) %>%
  dplyr::select(cancer,symbol,status) %>%
  tidyr::spread(key=cancer,value=status) ->stage.sur

stage.st %>%
  dplyr::full_join(stage.sub,by="symbol") %>%
  dplyr::full_join(stage.sur,by="symbol") ->clinical_result

clinical_result %>%
  #as.data.frame() %>%
  write.table(file=file.path(clinical_path,"clinical_result.txt"),sep="\t",quote = F)

clinical_result %>%
  dplyr::select(symbol) %>%
  dplyr::inner_join(gene_list,by="symbol") ->gene_list_format

gene_list_format %>%
  write.table(file=file.path(clinical_path,"gene_list_format.txt"),sep="\t",quote = F)

save.image(file = file.path(clinical_path,"clinical_merge.Rdata"))
