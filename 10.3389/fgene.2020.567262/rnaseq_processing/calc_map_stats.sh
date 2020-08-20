#!/bin/bash

## Run this script as
## 
## $bash calc_stats.sh camara
## $bash calc_stats.sh crivularis
##



## Bash arguments
## Set up 'camara' or 'crivularis'.
##------------------------------------------------------------------------------



source ../local/global.sh
JOB_ROOT=${PROJECT_ROOT}/rnaseq_processing
BAM_DIR=${JOB_ROOT}/results/bam
RC_DIR=${JOB_ROOT}/results/readclassify
COUNT_DIR=${JOB_ROOT}/results/counts
TMP_DIR=${JOB_ROOT}/tmp

TAGS=(orig common other)
GENOMES=(camara crivularis)



cd ${TMP_DIR}/star
echo "Number of uniquely mapped reads on A-genome:"
for seq_lib in ${RNASEQ_LIBS[@]}; do
    map_log_file=${seq_lib}__on__camara/Log.final.out
    cat ${map_log_file} | awk 'BEGIN{FS="\t"}
    $1 ~ /Number of input reads/ {i = $2}
    $1 ~ /Uniquely mapped reads number/{o = $2}
    END{print i, o}' 
done

echo "Number of uniquely mapped reads on R-genome:"
for seq_lib in ${RNASEQ_LIBS[@]}; do
    map_log_file=${seq_lib}__on__crivularis/Log.final.out
    cat ${map_log_file} | awk 'BEGIN{FS="\t"}
    $1 ~ /Number of input reads/ {i = $2}
    $1 ~ /Uniquely mapped reads number/{o = $2}
    END{print i, o}' 
done
cd -




cd ${RC_DIR}
echo "Number of read pairs after read classificaton with HomeoRoq."

for seq_lib in ${RNASEQ_LIBS[@]}; do
    for ref_genome in ${GENOMES[@]}; do
        for tag in ${TAGS[@]}; do
            genome_fasta=${GENOMELEAFDIR}/${ref_genome}/genome.fa
            bam_file=${seq_lib}__on__${ref_genome}.filtered.rc_${tag}.sorted.bam
            java -jar ${BIN}/picard.jar CollectAlignmentSummaryMetrics R=${genome_fasta} I=${bam_file} O=${bam_file}.xls
        done
    done
done



for ref_genome in ${GENOMES[@]}; do
    for tag in ${TAGS[@]}; do
        echo ${ref_genome}" - "${tag}
        for seq_lib in ${RNASEQ_LIBS[@]}; do
            bam_stats=${seq_lib}__on__${ref_genome}.filtered.rc_${tag}.sorted.bam.xls
            cat ${bam_stats} | awk 'BEGIN{FS="\t"} $1 ~ /^FIRST_OF_PAIR/ {print $9}'
        done
    done
done













