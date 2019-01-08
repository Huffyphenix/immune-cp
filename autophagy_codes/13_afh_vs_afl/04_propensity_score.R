
# library -----------------------------------------------------------------

library(ggplot2)
library(magrittr)
library(methods)

# path --------------------------------------------------------------------

tcga_path <- "/home/cliu18/liucj/projects/6.autophagy/TCGA"
afhl_path <- "/home/cliu18/liucj/projects/6.autophagy/09_afh_vs_afl"
afhl_class <- file.path(afhl_path, "01_af_h_l_classification")
props_path <- file.path(afhl_path, "04_props")
setwd(props_path)
scripts.dir <- props_path
# load data ---------------------------------------------------------------
library(dummies)
source("/home/cliu18/liucj/github/RstudioWithGit/autophagy_codes/13_afh_vs_afl/yy_code/cal.r")

clinical <- readr::read_rds(file.path(tcga_path, "pancan34_clinical.rds.gz"))
p62_sample_classification <- readr::read_rds(path = file.path(afhl_class, '.rds_01_p62_sample_classification.rds.gz')) %>% 
  dplyr::select(-barcode) %>% 
  tidyr::nest(-cancer_types)

# ----
# expr <- readr::read_rds(file.path(tcga_path, "pancan33_expr.rds.gz")) 
# expr %>% 
#   dplyr::mutate(expr = purrr::map(
#     .x = expr,
#     .f = function(.x) {
#       .x %>% tidyr::unite("gene", symbol, entrez_id, sep = "|") %>% as.data.frame() -> mRNAseq
#       
#       tibble::tibble(names = colnames(mRNAseq)[-1]) %>% 
#         dplyr::filter(stringr::str_sub(names, 14, 15) == "01") %>% 
#         dplyr::mutate(sample = stringr::str_sub(names, 1, 12)) %>% 
#         dplyr::distinct(sample, .keep_all = T) %>% 
#         tibble::deframe() -> uni_names
#       if (length(uni_names) == 0) return(NULL)
#       
#       mRNAseq <- mRNAseq[, c("gene", names(uni_names))]
#       names(mRNAseq) <- c("gene", uni_names)
#       mRNAseq.pri <- data.frame(t(mRNAseq[,2:ncol(mRNAseq)]))
#       colnames(mRNAseq.pri) <- as.vector(mRNAseq[,1])
#       mRNAseq.pri <- rm.zero.col(mRNAseq.pri)
#     }
#   )) %>% 
#   dplyr::filter(purrr::map_lgl(expr, Negate(is.null))) -> expr_dum
# 
# readr::write_rds(expr_dum, file.path(props_path, ".rds_01_expr_dum.rds.gz"), compress = "gz")

# ----

# propensity score --------------------------------------------------------

expr_dum <- readr::read_rds(file.path(props_path, ".rds_01_expr_dum.rds.gz"))
sum.mRNAAll <- data.frame()
analysis <- "analysis"
clinical %>% 
  dplyr::inner_join(expr_dum, by = "cancer_types") %>% 
  dplyr::inner_join(p62_sample_classification, by = "cancer_types") %>% 
  dplyr::mutate(
    clinical = purrr::map(
      .x = clinical,
      .f = function(.x){
        # not take pathologic stage into account.
        .x %>% 
          dplyr::select(barcode, age_at_initial_pathologic_diagnosis,gender,race) %>% 
          dplyr::filter(!race %in% c("ASIAN", NA_character_)) %>% 
          dplyr::mutate(race = as.factor(race)) %>% 
          dplyr::rename(sample = barcode)
      }
    )
  ) %>%
  dplyr::mutate(
    merge = purrr::map2(
      .x = clinical,
      .y = data,
      .f = function(.x, .y) {
        .y %>% 
          dplyr::mutate(myclusters = ifelse(rppa >= median(rppa), 0, 1)) %>% 
          dplyr::select(sample, Z = myclusters) %>% 
          dplyr::inner_join(.x, by = "sample") 
      }
    )
  ) %>% 
  dplyr::mutate(x = purrr::pmap(
    .l = list(
      .x = merge, 
      .y = expr,
      .z = cancer_types
    ),
    .f = function(.x, .y, .z){
      data <- .x %>% as.data.frame()
      rownames(data) <- data[,1]
      data <- data[,-1]
      mRNAseq.pri <- .y
      cancer <- .z
      
      dummy.feature <- setdiff(colnames(.x),c("Z","age_at_initial_pathologic_diagnosis"))
      data.dum <- dummy.data.frame(data, names = dummy.feature)
      dummy.list <- attr(data.dum,"dummies")
      rm.col <- c()
      for (i in 1:length(dummy.list))
      {
        rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
      }
      data.dum <- data.dum[,-rm.col]
      data.dum$X0 <- rep(1, nrow(data.dum))
      exclude.col <- match(c("Z","X0"), colnames(data.dum))
      colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
      form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
      
      .path <- file.path(props_path, cancer)
      if (!file.exists(.path)) dir.create(.path)
      
      mRNAseq.result <- weight.test(
        data.dum, form, mRNAseq.pri, is.continuous=TRUE,
        weight=ifelse(analysis=="myclusters","MW","ATT"), mirror.plot=FALSE, cancer, data.type = "mRNAseq", 
        outdir = .path, perm=FALSE)
      
      sum.mRNA <- summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE)
      summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE, cutoff=0.05)
      
      write.summary(sum.mRNA, cancer, analysis, "mRNA")
      write.result(mRNAseq.result, cancer,analysis, "mRNA")
      save(mRNAseq.result, file = paste(cancer,"_result.RData",sep=""))
      
      
      
      if(length(which(mRNAseq.result$fdr < 0.05)) > 0){
        
        sum.mRNA <- data.frame(sum.mRNA)
        sum.mRNA$class <- rep(cancer,times=nrow(sum.mRNA))
        if(nrow(sum.mRNAAll) == 0){
          sum.mRNAAll <- sum.mRNA
        }else{
          sum.mRNAAll <- rbind(sum.mRNAAll,sum.mRNA)
        }
        library(doMC)
        library(foreach)
        registerDoMC(40)
        perm.result <- foreach (seed =1:100, .combine='rbind') %dopar%
        {
          ## mRNAseq, only for KIRC        
          perm.mRNAseq.result <- weight.test(data.dum, form, mRNAseq.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, "mRNA", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=TRUE, seed=seed)
          perm.sum.mRNA <- summarize.fdr(mRNAseq.pri, perm.mRNAseq.result)
          
          write(c(seed, perm.sum.mRNA$n.sig), file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_count_",seed,".txt", sep=""), ncolumns=7)
          save(seed,perm.mRNAseq.result, file=paste(scripts.dir, "/",cancer,"_",analysis,"/perm_sig_result_",seed,".RData", sep=""))
          
        }
        
        
        cutoff <- 0.05
        seedV <- 1:100
        perm.cal(cancer, analysis, "mRNAseq", mRNAseq.pri, cutoff=cutoff, seedV=seedV)
        
        if(FALSE)
        {
          mRNA.ttest <- myttest(data.dum, mRNAseq.pri, cancer,"mRNA")
          sum.mRNA <- summarize.fdr(mRNAseq.pri, mRNA.ttest)
          save(mRNAseq.ttest, file=paste(cancer,"_mRNA_ttest.RData", sep=""))
        }
      }
    }
  ))
write.table(sum.mRNAAll,file="mRNAseq.genes.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)


# save --------------------------------------------------------------------

save.image(file.path(props_path, 'rda_01_propensity_score_mrna.rda'))
