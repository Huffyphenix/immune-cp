#!/usr/bin/Rscript
#输入表达量文件，行名基因，列名样本名
args<-commandArgs(TRUE)
Row_data_noiseq_RNA<-read.table(args[1],header=T,sep='\t',row.names=1)
for(i in 1:nrow(Row_data_noiseq_RNA)){
rownames(Row_data_noiseq_RNA)[i]<-strsplit(rownames(Row_data_noiseq_RNA)[i],'[|]')[[1]][2]}
#按照样本名分类样本，即需要做差异的两类样本
#Tumor_sample<-grep("0.",substr(colnames(Row_data_noiseq_RNA),14,15))
Control_sample<-colnames(Row_data_noiseq_RNA)[grep("1.",substr(colnames(Row_data_noiseq_RNA),14,15))]
if(args[9]=="A"){
Tumor_sample<-grep("0.",substr(colnames(Row_data_noiseq_RNA),14,15))}else if(args[9]=="B"){
	Tumor_sample<-c();j=1
	for(i in Control_sample){
		patient<-substr(i,9,12)
		if(length(grep(patient,colnames(Row_data_noiseq_RNA)))>=2){
		grep(paste(patient,".0.",sep=""),colnames(Row_data_noiseq_RNA),value=T)->Tumor_sample[j]
		j=j+1
	}else{
		Control_sample[-grep(patient,Control_sample)]->Control_sample
	}
}}
Case<-Row_data_noiseq_RNA[,Tumor_sample]
Control<-Row_data_noiseq_RNA[,Control_sample]

col_case<-ncol(Case)
col_con<-ncol(Control)
sample_info<-matrix(c(col_case,col_con),nrow=1)
rownames(sample_info)=args[2]
colnames(sample_info)=c("num_case","num_con")

Tmp_case<-as.data.frame(apply(t(Case),1,function(x) as.integer(x)))
Tmp_con<-as.data.frame(apply(t(Control),1,function(x) as.integer(x)))


test_data_noiseq_RNA<-data.frame(Tmp_con,Tmp_case)
row.names(test_data_noiseq_RNA)<-row.names(Row_data_noiseq_RNA)
x<-c(paste("con_",1:col_con,sep=""))
y<-c(paste("case_",1:col_case,sep=""))
library(NOISeq)

myfactors <- data.frame(Tissue = c(rep("con",col_con),rep("case",col_case)),TissueRun = c(x,y))
mydata <- readData(data = test_data_noiseq_RNA,factors = myfactors)

#设置cpm，为1时，筛除表达量全小于1的基因不对其进行差异表达
#mynoiseq <- noiseqbio(mydata, k = 0.5, norm = "n", factor = "Tissue",lc = 0, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345,filter = 1,cv.cutoff = 100, cpm = as.numeric(as.character(args[6])))
CPM=as.numeric(as.character(args[6]))
mynoiseq <- noiseqbio(mydata, k = 0.5, norm = "n", factor = "Tissue",lc = 0, r = 20, adj = 1.5, plot = FALSE, a0per = 0.9, random.seed = 12345,filter = 1,cv.cutoff = 100, cpm =CPM)


#提取差异表达基因
noiseq_result_rna<-mynoiseq@results[[1]] 
noiseq_degene_rna<-subset(noiseq_result_rna,noiseq_result_rna$prob>=1-as.numeric(as.character(args[4]))) #set FDR, 0.99=1-FDR
noiseq_up_rna<-subset(noiseq_degene_rna,noiseq_degene_rna$log2FC>=args[5]) #set fold change
noiseq_down_rna<-subset(noiseq_degene_rna,noiseq_degene_rna$log2FC<=(-as.numeric(as.character(args[5])))) #set fold change

#output results
result_path<-paste(args[3],"/",args[2],"/",args[2],"-",sep="")
sample_info_path<-paste(result_path,"sample_info",sep="")
write.table(sample_info,sample_info_path,quote = F,sep='\t',row.names=T)
all_result<-paste(result_path,"cpm_",args[6],sep="")
write.table(noiseq_result_rna,all_result,quote = F,sep='\t',row.names=T)
FDR<-1-as.numeric(as.character(args[4]))
up_res<-paste(result_path,"cpm_",args[6],"-FC_",args[5],"-FDR_",FDR,"-up",sep="")
write.table(noiseq_up_rna,up_res,quote = F,sep='\t',row.names=T)

down_res<-paste(result_path,"cpm_",args[6],"-FC_",args[5],"-FDR_",FDR,"-down",sep="")
write.table(noiseq_down_rna,down_res,quote = F,sep='\t',row.names=T)

all_res<-paste(result_path,"cpm_",args[6],"-FC_",args[5],"-FDR_",FDR,"-all",sep="")
noiseq_de_gene<-rbind(noiseq_up_rna,noiseq_down_rna)
write.table(noiseq_de_gene,all_res,quote = F,sep='\t',row.names=T)

#extract interest gene list DE result
gene_list<-as.character(as.vector(t(read.table(args[7],sep='\t'))))
gene_res<-paste(result_path,"cpm_",args[6],"-FC_",args[5],"-FDR_",FDR,"-",args[8],sep="")
write.table(noiseq_de_gene[intersect(rownames(noiseq_de_gene),gene_list),],gene_res,quote = F,sep='\t',row.names=T)

gene_res_all<-paste(result_path,"cpm_",args[6],"-",args[8],sep="")
write.table(noiseq_result_rna[intersect(rownames(noiseq_result_rna),gene_list),],gene_res_all,quote = F,sep='\t',row.names=T)
