library(dplyr)

#########################
#read in De genes
#########################
DE_path<-"/project/huff/huff/immune_checkpoint/1.DE/B/result_count"
DE_data<-read.table(file.path(DE_path,"samples_over_5/symbol-gene-FC-trend-cancer.txt"),sep="\t")
colnames(DE_data)=c("symbol","gene_id","log2FC","trend","cancer")
#get cancer type.
DE_data %>%
  dplyr::select(cancer) %>%
  unique()%>% t() %>% as.vector()->cancer_type

################################
#do kruskal.test for each gene De in cancer.
################################
#checpoint<-read.table("/project/huff/huff/immune_checkpoint/checkpoint/all.entrez_id-gene_id",sep="\t",header=T)
symbol<-c()
for(i in 1:nrow(DE_data)){
        c(symbol,paste(DE_data[i,"symbol"],DE_data[i,"gene_id"],sep="|"))->symbol}
unique(symbol)->symbol
TCGA_path<-"/project/huff/huff/TCGA_survival/data/"
result<-matrix(NA,nrow = length(symbol),ncol = length(cancer_type))
colnames(result)=cancer_type
rownames(result)=unique(DE_data[,"gene_id"])
#FC result
FCresult<-matrix(NA,nrow = length(symbol),ncol = length(cancer_type))
colnames(FCresult)=cancer_type
rownames(FCresult)=unique(DE_data[,"gene_id"])

for(cancer in cancer_type){
#  DE_data[DE_data[,"cancer"]==cancer,]->symbol_tmp
#  symbol<-c()
#  for(i in 1:nrow(symbol_tmp)){
#	c(symbol,paste(symbol_tmp[i,"symbol"],symbol_tmp[i,"gene_id"],sep="|"))->symbol}
    exp_tmp<-t(read.table(paste(TCGA_path,cancer,
                              "/gdac.broadinstitute.org_",
                              cancer,".mRNAseq_Preprocess.Level_3.2016012800.0.0/",
                              cancer,".uncv2.mRNAseq_RSEM_all.txt",sep=""),
                              header = T,sep = "\t",row.names = 1))
    exp_tmp<-exp_tmp[grep("0.",substr(rownames(exp_tmp),14,15)),]
    exp_tmp<-data.frame(sample=substr(rownames(exp_tmp),9,12),exp_tmp[,symbol])
    for(i in 2:ncol(exp_tmp)){
        colnames(exp_tmp)[i]<-strsplit(colnames(exp_tmp)[i],"\\.")[[1]][length(strsplit(colnames(exp_tmp)[i],"\\.")[[1]])]}

    clin_tmp<-t(read.table(paste(TCGA_path,cancer,"/gdac.broadinstitute.org_",
                               cancer,".Clinical_Pick_Tier1.Level_4.2016012800.0.0/All_CDEs.txt",
                               sep=""),sep="\t",header = T,row.names = 1,fill=T,quote=""))
    clin_need<-c("patient_id","pathologic_stage")
    clin_tmp<-clin_tmp[,clin_need]
    for(i in 1:nrow(clin_tmp)){
	if(!is.na(clin_tmp[i,"pathologic_stage"])){
	if(clin_tmp[i,"pathologic_stage"]=="stage i" | clin_tmp[i,"pathologic_stage"]=="stage ia" | clin_tmp[i,"pathologic_stage"]=="stage ib" | clin_tmp[i,"pathologic_stage"]=="stage ic"){
	clin_tmp[i,"pathologic_stage"]="StageI"}
	else if(clin_tmp[i,"pathologic_stage"]=="stage ii" | clin_tmp[i,"pathologic_stage"]=="stage iia" | clin_tmp[i,"pathologic_stage"]=="stage iib"| clin_tmp[i,"pathologic_stage"]=="stage iic"){
	clin_tmp[i,"pathologic_stage"]="StageII"}
	else if(clin_tmp[i,"pathologic_stage"]=="stage iii" | clin_tmp[i,"pathologic_stage"]=="stage iiia" |clin_tmp[i,"pathologic_stage"]=="stage iiib"|clin_tmp[i,"pathologic_stage"]=="stage iiic"){
	clin_tmp[i,"pathologic_stage"]="StageIII"}
	else if(clin_tmp[i,"pathologic_stage"]=="stage iv"|clin_tmp[i,"pathologic_stage"]=="stage iva"|clin_tmp[i,"pathologic_stage"]=="stage ivb"|clin_tmp[i,"pathologic_stage"]=="stage ivc"){
	clin_tmp[i,"pathologic_stage"]="StageIV"}
      }
   }
    toupper(clin_tmp[,"patient_id"])->clin_tmp[,"patient_id"]
    colnames(clin_tmp)<-c("sample","stage")
    grep("TRUE",is.na(clin_tmp[,"stage"]))->del
    if(length(del)!=dim(clin_tmp)[1]){
       if(length(del)!=0){clin_tmp[-del,]->clin_tmp}else{clin_tmp->clin_tmp}
       merge(clin_tmp,exp_tmp,by="sample")->data_tmp 
       for(gene in colnames(data_tmp)[3:ncol(data_tmp)]){
         -round(log10(kruskal.test(data_tmp[,gene] ~ data_tmp[,"stage"])$p.value),3)->p
	 round(kruskal.test(data_tmp[,gene] ~ data_tmp[,"stage"])$p.value,4)->pvalue
	 res<-c()
	 for(i in levels(data_tmp$stage)){
		 res<-c(res,median(data_tmp[grep("TRUE",data_tmp[,"stage"]==i),gene]))}
	 if(min(res)!=0){FCresult[gene,cancer]<-max(res)/min(res)}else{FCresult[gene,cancer]<-max(res)/0.5}
	 if(p>=1 && FCresult[gene,cancer]>=1.5){
	    result[gene,cancer]<-p
	    unique(DE_data[grep("TRUE",DE_data[,"gene_id"]==gene),"symbol"])->gene_s
	    jpegname<-paste("plot/",cancer,"_",gene_s,".tif",sep="")
            tiff(jpegname,width = 480, height = 480)
	    boxplot(log2(data_tmp[,gene]) ~ data_tmp[,"stage"],main=paste(cancer,"_",gene_s,sep=""),col = rainbow(length(levels(data_tmp[,"stage"]))),xlab = "Stage",ylab = "log2(Gene Expression)")
	    FC<-round(FCresult[gene,cancer],3)
	    text(1,max(log2(data_tmp[,gene]))-1,c(paste("P.value=",pvalue,"\n","Maximal FC=",FC)))
	    dev.off()
	}
     }}else{print(paste("WORNING:",cancer,"have no stage data to use!"))}
}
result[is.na(result)]<-0
write.table(result,"log10pva-DEcheckpoint.719",sep="\t",quote=F)
write.table(FCresult,"FC-DEcheckpoint.719",sep="\t",quote=F)
