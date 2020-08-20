#!/bin/bash

## Run this script as
##
## $bash create_go_annotations.sh
##


source ../local/global.sh
JOB_ROOT=${PROJECT_ROOT}/data_analysis
GODATADIR=${JOB_ROOT}/data/go

mkdir -p ${GODATADIR}
cd ${GODATADIR}



date '+%Y-%m-%d'
## 2017-04-19

wget ftp://ftp.arabidopsis.org/Ontologies/Gene_Ontology/ATH_GO_GOSLIM.txt.gz
gunzip ATH_GO_GOSLIM.txt.gz
wget http://geneontology.org/ontology/go.obo
python ${SCRIPTDIR}/create_go_annotations.py --obo go.obo \
            --goslim ATH_GO_GOSLIM.txt \
            --ann ${GENOMEDIR}/chirsuta/chi_m25.txt \
            --go2carhr go2carhr.tsv 

date '+%Y-%m-%d'
## 2017-04-19



