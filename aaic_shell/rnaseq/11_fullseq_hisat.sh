#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -N WCT_FUL_HISAT_CS_COUNT1k

nCPU=16
PROJECT_DIR=/home/jqsun/research/wheat_cold_treatment

BIN=/home/jqsun/local/bin
UTILS=/home/jqsun/local/utilfunc

BATCH_DIR=${PROJECT_DIR}/data/fullseq_hisat
BAM_DIR=${BATCH_DIR}/bam
FASTQ_DIR=${PROJECT_DIR}/data/clean_fullseq

GENOME_DIR=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2


COUNTS_DIR=${BATCH_DIR}/counts_iwgsc
GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3

COUNTS_DIR=${BATCH_DIR}/counts_1.0k
GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.gff3

pHISAT=0
pCount=1



cd ${PROJECT_DIR}



if [ ${pHISAT} -eq 1 ]; then
    
    mkdir -p ${BAM_DIR}
    
    cd ${FASTQ_DIR}
    for fastq_filepath in `ls *_R1.cleaned.fastq.gz | grep _CS_ `; do
        echo ${fastq_filepath}
        
        fastq_filename=`basename ${fastq_filepath} _R1.cleaned.fastq.gz`
        
        # mapping
        time ${BIN}/hisat2 -p ${nCPU} -x ${GENOME_INDEX} \
               -1 ${fastq_filename}_R1.cleaned.fastq.gz -2 ${fastq_filename}_R2.cleaned.fastq.gz \
               -S ${BAM_DIR}/${fastq_filename}.sam
        
        # sorting
        cd ${BAM_DIR}
        ${BIN}/samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}.sam
        ${BIN}/samtools index -c ${fastq_filename}.bam
        
        # remove sam
        rm ${fastq_filename}.sam
        
        cd -
    done

fi




if [ ${pCount} -eq 1 ]; then
    mkdir -p ${COUNTS_DIR}
    cd ${BAM_DIR}

    for bam_filepath in `ls *.bam`; do
        bam_filename=`basename ${bam_filepath} .bam`

        ${BIN}/featureCounts -t gene -g ID -a ${GTF_FILEPATH} \
               -o ${COUNTS_DIR}/${bam_filename}.counts.gene.tsv ${bam_filepath}
        ${BIN}/featureCounts -t exon -g Parent -a ${GTF_FILEPATH} \
               -o ${COUNTS_DIR}/${bam_filename}.counts.transcript.tsv ${bam_filepath}
    done
 

    cd ${COUNTS_DIR}
    python ${UTILS}/merge_counts.py . counts.gene.tsv all.counts.gene.tsv
    python ${UTILS}/merge_counts.py . counts.transcript.tsv all.counts.transcript.tsv
    python ${UTILS}/merge_counts.py . counts.gene.tsv.summary all.counts.gene.tsv.summary
    python ${UTILS}/merge_counts.py . counts.transcript.tsv.summary all.counts.transcript.tsv.summary
    gzip all.counts.gene.tsv
    gzip all.counts.transcript.tsv
    gzip all.counts.gene.tsv.summary
    gzip all.counts.transcript.tsv.summary

fi




