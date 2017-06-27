#!/bin/sh
OPTIND=1
SCRIPTS_HOME=/project/huff/huff/immune_checkpoint/script/1.DE
BASE_PATH=$PWD
CPM=1
FDR=0.01
log2FC=1
SAMPLE_PAIR="B"
usage() {
    cat <<EOF
Usage:bash $0 -e <exp_matrix> -t <cancer type> -g <gene list> -l <label of gene list> -c [CPM] -p [FDR] -f [log2 FC] -d [workdir] -s <DE_type>
        -e <string> Expression matrix of the genes, colnames are gene id, rownames are sample id. Seperated by tab.
	-t <cancer type> cancer type/abbreviation of cancer.
        -l <string> label of your gene list, some label you want to show at the end of your result filename.
	-c [numeric] CPM of noiseq, whose sum of expression < CPM*num(samples) will be filtered.
        -p [numeric] Cutoff pvalue of FDR of noiseq, range from 0 to 1. Default is 0.01.
        -g <list> gene list of your interest. Devided by \\n.
	-f [numeric] Log2 FC of noiseq.
	-d [string] working direcotry. Default is $PWD.
	-s <DE type> A: All disease samples compare with normal samples.
		     B: Paired samples to do differentiate expression.
		     Default is B type.
EOF
    exit 1
}

while getopts :e:t:g:l:c:p:f:d:s: ARGS; do
    case $ARGS in
                e) EXP_PATH=$OPTARG ;;
		t) CANCER=$OPTARG ;;
                g) GENE_LIST=$OPTARG ;;
		l) LABEL=$OPTARG ;;
		c) CPM=$OPTARG ;;
                p) FDR=$OPTARG ;;
                f) log2FC=$OPTARG ;;
                d) BASE_PATH=$OPTARG ;;
		s) SAMPLE_PAIR=$OPTARG ;;
        *) usage ;;
    esac
done

if [ "$EXP_PATH" == "" -o "$GENE_LIST" == "" -o "$CANCER" == "" ];
then
    usage
fi

mkdir $BASE_PATH/$CANCER
#echo $EXP_PATH $CANCER $BASE_PATH $FDR $log2FC $CPM $GENE_LIST $LABEL

Rscript $SCRIPTS_HOME/noiseq.R $EXP_PATH $CANCER $BASE_PATH $FDR $log2FC $CPM $GENE_LIST $LABEL $SAMPLE_PAIR

