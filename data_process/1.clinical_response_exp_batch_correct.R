library(sva)
library(magrittr)
library(edgeR)

basic_path <- file.path("/home/huff/project")

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


# load expression count data ----------------------------------------------
exp_count <- readr::read_tsv(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp","all_count_expression_2.txt"))
exp_count.df <- as.data.frame(exp_count)[,-1]
rownames(exp_count.df) <- exp_count$gene_id

sample_info <- readr::read_tsv(file.path(basic_path,"immune_checkpoint/clinical_response_data","RNAseq-sample_info_complete.tsv")) %>%
  dplyr::filter(Library_strategy == "RNA-Seq")%>%
  dplyr::select(Run,blockade,Cancer.y,Biopsy_Time,Response,Gender,Bioproject) 

# expression normalization 
# CPM expression

exp_count.df.normalized <- calcNormFactors(DGEList(counts=as.matrix(exp_count.df)), method="upperquartile")
normalized_CPM_expr = cpm(exp_count.df.normalized$counts)

# remove genes with variance == 0 in each batch.
gene_no_variance <- c()
for (i in 1:length(unique(sample_info$batch))) {
  sample_info.df %>%
    dplyr::mutate(Run = rownames(sample_info.df)) %>%
    dplyr::filter(batch == i) %>%
    .$Run -> .Run
  exp_count.df[,.Run] -> .exp_count
  
  .mad <- apply(.exp_count,1,mad) 
  
  .mad %>%
    as.data.frame() %>%
    dplyr::mutate(gene = names(.mad)) %>%
    dplyr::filter(.==0) %>%
    .$gene -> .gene_no_variance
  
  c(.gene_no_variance,gene_no_variance) %>% unique()-> gene_no_variance
} 

exp_count.df.filter <- exp_count.df[setdiff(rownames(exp_count.df),gene_no_variance),]
normalized_CPM_expr.filter <- normalized_CPM_expr[setdiff(rownames(normalized_CPM_expr),gene_no_variance),]

# PCA analysis before batch correct -----------------------------------------------------------
library("Rtsne")

exp_count.df.filter %>%
  dplyr::mutate(gene_id = rownames(exp_count.df.filter)) %>%
  tidyr::gather(-gene_id,key="Run",value="count") %>%
  tidyr::spread(key="gene_id",value="count") %>%
  dplyr::inner_join(sample_info,by="Run")-> for_tsne

normalized_CPM_expr.filter %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(gene_id = rownames(normalized_CPM_expr.filter)) %>%
  tidyr::gather(-gene_id,key="Run",value="count") %>%
  tidyr::spread(key="gene_id",value="count") %>%
  dplyr::inner_join(sample_info,by="Run")-> for_tsne

# normalization
data_for_tSNE.matx.normalize <- normalize_input(for_tsne[,-c(1,(ncol(for_tsne)-6):ncol(for_tsne))] %>% as.matrix())

# do tSNE
tsne_res <- Rtsne(data_for_tSNE.matx.normalize,dims=2,pca=FALSE,theta=0.0)

# Show the objects in the 2D tsne representation
plot(tsne_res$Y,col=factor(for_tsne$Bioproject), asp=1)
tsne_res$Y %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  # dplyr::rename("tSNE 1"="V1","tSNE 2"="V2") %>%
  dplyr::mutate(Bioproject = for_tsne$Bioproject,
                blockage = for_tsne$blockade,
                Response = for_tsne$Response) %>%
  ggplot(aes(x=V1,y= V2)) +
  geom_jitter(aes(color=Bioproject,shape=blockage),size=1) +
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  ggpubr::color_palette("jco") +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp","Before_Batch_corrected_counts_exp.png"),device = "png",width = 4,height = 3)

ggsave(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp","Before_Batch_corrected_CPM_exp.png"),device = "png",width = 4,height = 3)

# 已知的batch effect可能有Bioproject 和 blockade


# known batch -------------------------------------------------------------
sample_info %>%
  dplyr::inner_join(tibble::tibble(
    Bioproject=c("PRJEB25780", "PRJNA312948", "PRJNA356761", "PRJNA476140","PRJNA82747"),
    batch = c(1,2,3,4,5))) -> sample_info
sample_info.df <- as.data.frame(sample_info)[,-1]
rownames(sample_info.df) <- sample_info$Run

tibble::tibble(sample=rownames(sample_info.df)) %>%
  dplyr::filter(sample %in% colnames(exp_count.df)) %>%
  .$sample -> common_sample
sample_info.df[common_sample,] -> sample_info.df

exp_count.df.filter[,common_sample] -> exp_count.df

tibble::tibble(sample=rownames(sample_info.df)) %>%
  dplyr::filter(sample %in% colnames(normalized_CPM_expr.filter)) %>%
  .$sample -> common_sample.cpm

normalized_CPM_expr.filter[,common_sample.cpm] -> exp_count.df.cpm

# 应用ComBat功能去除已知的batch effect ------------------------------------------

modcombat <-  model.matrix(~1, data = sample_info.df[,c("Bioproject", "batch")])
batch = sample_info.df$batch 
combat_edata = ComBat(dat=as.matrix(exp_count.df), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)


combat_edata.CPM = ComBat(dat=as.matrix(exp_count.df.cpm), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

# PCA analysis after batch correct -----------------------------------------------------------
library("Rtsne")

combat_edata.CPM %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(gene_id = rownames(exp_count.df.filter)) %>%
  tidyr::gather(-gene_id,key="Run",value="count") %>%
  tidyr::spread(key="gene_id",value="count") %>%
  dplyr::inner_join(sample_info,by="Run")-> for_tsne.after

combat_edata %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  dplyr::mutate(gene_id = rownames(exp_count.df.filter)) %>%
  tidyr::gather(-gene_id,key="Run",value="count") %>%
  tidyr::spread(key="gene_id",value="count") %>%
  dplyr::inner_join(sample_info,by="Run")-> for_tsne.after

# normalization
data_for_tSNE.matx.normalize.after <- normalize_input(for_tsne.after[,-c(1,(34692-6):34692)] %>% as.matrix())

# do tSNE
tsne_res.after <- Rtsne(data_for_tSNE.matx.normalize.after,dims=2,pca=FALSE,theta=0.0)

# Show the objects in the 2D tsne representation
plot(tsne_res.after$Y,col=factor(tsne_res.after$Bioproject), asp=1)
tsne_res.after$Y %>%
  as.data.frame() %>%
  dplyr::as.tbl() %>%
  # dplyr::rename("tSNE 1"="V1","tSNE 2"="V2") %>%
  dplyr::mutate(Bioproject = for_tsne.after$Bioproject,
                blockage = for_tsne.after$blockade,
                Response = for_tsne.after$Response) %>%
  ggplot(aes(x=V1,y= V2)) +
  geom_jitter(aes(color=Bioproject,shape=blockage),size=1) +
  xlab("tSNE 1") +
  ylab("tSNE 2") +
  ggpubr::color_palette("jco") +
  my_theme +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp","After_Batch_corrected_counts_exp.png"),device = "png",width = 4,height = 3)

ggsave(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp","After_Batch_corrected_CPM_exp.png"),device = "png",width = 4,height = 3)

combat_edata %>%
  readr::write_rds(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp","Batch_corrected_counts_exp.rds.gz"))

combat_edata.CPM  %>%
  readr::write_rds(file.path("/home/huff/project/immune_checkpoint/clinical_response_data/mRNA_exp","Batch_corrected_CPM_exp.rds.gz"))
