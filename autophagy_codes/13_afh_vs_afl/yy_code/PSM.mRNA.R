# check data
rm(list=ls())
setwd("/extraspace/yye1/analysis/Hypoxia/PSM")
folder <- "mRNA"
if (!file.exists(folder)) { dir.create(folder) }
scripts.dir <- paste("/extraspace/yye1/analysis/Hypoxia/PSM/",folder,sep="")
setwd(scripts.dir)
analysis="myclusters" #Oxygen_Content
sum.mRNAAll <- data.frame()
####must exist stum  clinical, stratification, mRNA files.
mRNA.files.names <- list.files(path="/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",pattern="_20160513")
mRNA.files.Abs <- gsub("_mRNA_each_exp_20160513","",mRNA.files.names)
clinical.files.names <- list.files(path="/extraspace/TCGA/TCGA_clinical/",pattern=".txt")
clinical.files.Abs <- gsub("_clinical_clean.txt","",clinical.files.names)
stratification.files.names <- list.files(path="/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/",pattern=".txt")
stratification.files.Abs <- gsub(".Hypoxia.stratification.txt","",stratification.files.names)
SampleCount <- read.delim("/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/Hypoxia and normoxia sample size large than 30.txt",header=T)
stratification.files.Abs <- intersect(SampleCount$type,intersect(clinical.files.Abs,mRNA.files.Abs))
#stratification.files.Abs <- stratification.files.Abs[which(stratification.files.Abs != "OV" & stratification.files.Abs !="UCEC")]
for(cancer in stratification.files.Abs){
  clinical.file <- paste("/extraspace/TCGA/TCGA_clinical/",cancer,"_clinical_clean.txt", sep="")
  clinical.raw <- read.table(clinical.file, header=TRUE, check.names=FALSE, sep="\t",quote = "", comment.char="")
  print(paste("Total clinical samples for", cancer, ":", nrow(clinical.raw)))
  ###add hypoxic and normoxic stratification data
  stratification <- read.delim(paste("/extraspace/yye1/analysis/Hypoxia/Stratification.tumor/",cancer,".Hypoxia.stratification.txt",sep=""),header=T)
  stratification$barcode <- gsub("\\.","\\-",substr(stratification$SampleID,1,12))
  stratification <- stratification[which(stratification$myclusters %in% c("hypoxic","normoxic")),c("myclusters","barcode")]
  clinical.raw <- merge(stratification,clinical.raw,by="barcode")
  if(nrow(clinical.raw) >=50){
    end.col <- which(colnames(clinical.raw)=="os_days")-1
    #str(clinical.raw[,1:end.col])
    #summary(clinical.raw[,1:end.col])
    # apply(clinical.raw[1:end.col], 2, function(x) length(which(is.na(x))))
    
    # process data
    data <- clinical.raw[,1:end.col]
    ## rm race_ASIAN and race_NA
    consider_factors <- c("barcode","myclusters","age_at_initial_pathologic_diagnosis","gender","race","pathologic_stage")
    consider_factors <- intersect(consider_factors,colnames(data))
    if(length(table(data$pathologic_stage))==0){
      consider_factors <- consider_factors[which(consider_factors != "pathologic_stage")]
      data <- data[,consider_factors]
    }else{
      data <- data[,consider_factors]
    }
    
    if(length(which(is.na(data$race))) >0){
      data <- data[-which(data$race=="ASIAN" | is.na(data$race)),]
    }
    if(length(which(is.na(data$pathologic_stage))) > 0){
      data <- data[-which(is.na(data$pathologic_stage)),]
    }
    
    data$race <- as.vector(data$race)
    #data$race[which(is.na(data$race))] <- "UNKNOWN"
    data$race <- factor(data$race)
    data <- unique(data)
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    analysis <- "myclusters"
    if(analysis=="myclusters"){
      # convert hypoxic and normoxic to numeric 1,0 to suppress the warning message in lm
      data$myclusters <- ifelse(data$myclusters=="hypoxic",1,0)
      colnames(data)[which(colnames(data)=="myclusters")] <- "Z"
    }
    
    # convert to dummy
    library(dummies)
    dummy.feature <- setdiff(colnames(data),c("Z","age_at_initial_pathologic_diagnosis"))#,"pathologic_stage"))
    data.dum <- dummy.data.frame(data, names=dummy.feature)
    dummy.list <- attr(data.dum,"dummies")
    rm.col <- c()
    for (i in 1:length(dummy.list))
    {
      rm.col <- c(rm.col, dummy.list[[i]][length(dummy.list[[i]])])
    }
    data.dum <- data.dum[,-rm.col]
    data.dum$X0 <- rep(1, nrow(data.dum))
    #form <- as.formula("Z~.") # should exclude X0
    exclude.col <- match(c("Z","X0"), colnames(data.dum))
    colnames(data.dum) <- gsub(" ", ".", colnames(data.dum))
    form <- as.formula(paste("Z~",paste(colnames(data.dum)[-exclude.col],collapse="+"),sep=""))
    # perform calculation
    source("~/code/cal.R")
    library(doMC)
    library(foreach)
    registerDoMC(15)
    # mRNA.exp
    mRNAseq <- read.delim(paste("/extraspace/TCGA/TCGA_exp_DataPortal/mRNA_exp/",cancer,"_mRNA_each_exp_20160513",sep=""),header=T)
    mRNAseq <- mRNAseq[,colnames(mRNAseq)[-which(colnames(mRNAseq) %in% colnames(mRNAseq)[duplicated(substr(colnames(mRNAseq),1,12))] | colnames(mRNAseq)=="X")]]
    mRNAseq[,2:ncol(mRNAseq)] <- mRNAseq[,2:ncol(mRNAseq)]
    mRNAseq.pri <- data.frame(t(mRNAseq[,2:ncol(mRNAseq)]))
    colnames(mRNAseq.pri) <- as.vector(mRNAseq[,1])
    rownames(mRNAseq.pri) <- substr(rownames(mRNAseq.pri),1,12)
    mRNAseq.pri <- rm.zero.col(mRNAseq.pri)
    folder <- paste(cancer,"_",analysis,sep="")
    if (!file.exists(folder)) { dir.create(folder) }
    
    mRNAseq.result <- weight.test(data.dum, form, mRNAseq.pri, is.continuous=TRUE,weight=ifelse(analysis=="myclusters","MW","ATT"),mirror.plot=FALSE, cancer, data.type= "mRNAseq", outdir=paste(scripts.dir, "/",cancer,"_",analysis,sep=""),perm=FALSE)
    sum.mRNA <- summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE)
    summarize.fdr(mRNAseq.pri, mRNAseq.result, print=TRUE, cutoff=0.05)
    write.summary(sum.mRNA, cancer, analysis,"mRNA")
    write.result(mRNAseq.result, cancer, analysis,"mRNA")
    save(mRNAseq.result, file=paste(cancer,"_", analysis,"_result.RData",sep=""))
    if(length(which(mRNAseq.result$fdr < 0.05)) > 0){
      sum.mRNA <- data.frame(sum.mRNA)
      sum.mRNA$class <- rep(cancer,times=nrow(sum.mRNA))
      if(nrow(sum.mRNAAll) == 0){
        sum.mRNAAll <- sum.mRNA
      }else{
        sum.mRNAAll <- rbind(sum.mRNAAll,sum.mRNA)
      }
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
  
}

write.table(sum.mRNAAll,file="mRNAseq.genes.across.cancer.typesAll.txt",quote = F,sep="\t",row.names = F)

