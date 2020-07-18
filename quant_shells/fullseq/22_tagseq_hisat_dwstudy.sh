#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l memsz_job=64gb
#PBS -l elapstim_req=24:00:00
#PBS -N WCT_TAG_DWSTUDY_0.80



nCPU=16
PROJECT_DIR=/home/jqsun/research/wheat_cold_treatment

BATCH_DIR=${PROJECT_DIR}/data/clean_tagseq_0.80



BIN=/home/jqsun/local/bin
UTILS=/home/jqsun/local/utilfunc


BAM_DIR=${BATCH_DIR}/bam
FASTQ_DIR=${BATCH_DIR}/fastq
GENOME_DIR=/home/jqsun/research/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2
EX_LEN=1.0
COUNTS_DIR=${BATCH_DIR}/counts_${EX_LEN}k
GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_${EX_LEN}k.gff3



pHISAT=1
pCount=1


cd ${PROJECT_DIR}

if [ ${pHISAT} -eq 1 ]; then
    mkdir -p ${BAM_DIR}
    cd ${FASTQ_DIR}
    for fastq_filepath in `ls *.cleaned.fastq.gz `; do
        echo ${fastq_filepath}

        fastq_filename=`basename ${fastq_filepath} .cleaned.fastq.gz`

        # mapping
        time ${BIN}/hisat2 -p ${nCPU} --no-spliced-alignment -x ${GENOME_INDEX} -U ${fastq_filepath} -S ${BAM_DIR}/${fastq_filename}.sam

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
        
        ${BIN}/featureCounts -t gene -g ID -a ${GTF_FILEPATH} -s 1 \
               -o ${COUNTS_DIR}/${bam_filename}.counts.gene.tsv ${bam_filepath}
    done

    cd ${COUNTS_DIR}
    python ${UTILS}/merge_counts.py . counts.gene.tsv all.counts.gene.tsv
    python ${UTILS}/merge_counts.py . counts.gene.tsv.summary all.counts.gene.tsv.summary
    gzip all.counts.gene.tsv
    gzip all.counts.gene.tsv.summary

fi



