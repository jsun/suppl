#!/bin/bash

## Run this script as
## 
## $bash calc_stats.sh camara
## $bash calc_stats.sh crivularis
##



## Bash arguments
## Set up 'camara' or 'crivularis'.
##------------------------------------------------------------------------------
if [ $# -ne 1 ]; then
    echo 'Set the species name as camara or crivularis'
    exit 1
fi
species=$1



source ../local/global.sh
JOB_ROOT=${PROJECT_ROOT}/assembly
SPDIR=${JOB_ROOT}/${species}


echo "Number of uniquely mapped reads:"
for i in 1 2 3 4 5 6 7 8 9 10; do
    log_file=${SPDIR}/tmp${i}/Log.final.out
    cat ${log_file} | awk 'BEGIN{FS="\t"}
    # $1 ~ /Number of input reads/ {print $2}
    $1 ~ /Uniquely mapped reads number/{print $2}' 
done



cd ${SPDIR}/tmp10
echo "Number of assembled genes:"
python ${SCRIPTDIR}/htseq-count2.py -m divided -f bam -r pos -s reverse -t gene -i ID Aligned.out.sorted.bam ${GENOMELEAFDIR}/${species}/genome.modified.gff > mapped.counts.txt
awk 'BEGIN{FS="\t"}{if($1 ~ /CARHR/ && $2 > 0){num_of_genes++}}END{print num_of_genes}' mapped.counts.txt


java -jar ${BIN}/picard.jar CollectAlignmentSummaryMetrics R=genome.fa I=Aligned.out.sorted.bam O=Aligned.out.sorted.bam.summary.Xls
java -jar ${BIN}/picard.jar QualityScoreDistribution I=Aligned.out.sorted.bam O=Aligned.out.sorted.bam.qualscores.dist.txt CHART=Aligned.out.sorted.bam.qualscores.dist.pdf 







