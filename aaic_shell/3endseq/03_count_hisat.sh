#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=4:00:00
#PBS -N WGENOME_CNT



PROJECT_DIR=/home/jqsun/research/wheat_genome
cd ${PROJECT_DIR}

BIN=/home/jqsun/local/bin
UTILS=/home/jqsun/local/utilfunc

DATA_DIR=${PROJECT_DIR}/data
BAM_DIR=${DATA_DIR}/bam

GENOME_DIR=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS

EX_GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext.gff3
EX_COUNTS_DIR=${DATA_DIR}/counts

GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3
COUNTS_DIR=${DATA_DIR}/counts_iwgsc



mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
for bam_filepath in `ls *.bam`; do
    bam_filename=`basename ${bam_filepath} .bam`
     ${BIN}/featureCounts -t gene -g ID -a ${EX_GTF_FILEPATH} -s 1 \
                         -o ${EX_COUNTS_DIR}/${bam_filename}.counts.gene.tsv ${bam_filepath}
    
    ${BIN}/featureCounts -t gene -g ID -a ${GTF_FILEPATH} -s 1 \
                         -o ${COUNTS_DIR}/${bam_filename}.counts.gene.tsv ${bam_filepath}

done
    
cd ${COUNTS_DIR}

python ${UTILS}/merge_counts.py . counts.gene.tsv all.counts.gene.tsv
python ${UTILS}/merge_counts.py . counts.gene.tsv.summary all.counts.gene.tsv.summary
gzip all.counts.gene.tsv
gzip all.counts.gene.tsv.summary



