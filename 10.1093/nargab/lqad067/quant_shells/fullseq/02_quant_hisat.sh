#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l elapstim_req=72:00:00
#PBS -b 1
#PBS -N HB_FUL_QUANT_HISAT


nCPU=16


PROJECT_DIR=~/projects/HuoberBrezel
DATA_DIR=${PROJECT_DIR}/data/fullseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bam
COUNTS_DIR=${DATA_DIR}/counts
GENOME_DIR=~/projects/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2
GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706


cd ${PROJECT_DIR}


mkdir -p ${BAM_DIR}
cd ${CLEAN_FASTQ_DIR}

for fastq_filepath in `ls *_R1.cleaned.fastq.gz `; do
    echo ${fastq_filepath}
    fastq_filename=`basename ${fastq_filepath} _R1.cleaned.fastq.gz`
        
    time hisat2 -p ${nCPU} -x ${GENOME_INDEX} \
                       -1 ${fastq_filename}_R1.cleaned.fastq.gz \
                       -2 ${fastq_filename}_R2.cleaned.fastq.gz \
                       -S ${BAM_DIR}/${fastq_filename}.sam
        
    cd ${BAM_DIR}
    samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}.sam
    samtools index -c ${fastq_filename}.bam
    #rm ${fastq_filename}.sam
    cd -
done




mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
gff_versions=("iwgsc" "ext_0.5k" "ext_1.0k" "ext_1.5k" "ext_2.0k" "ext_2.5k" "ext_3.0k" "ext_3.5k" "ext_4.0k")
for gff in ${gff_versions[@]}; do
    featureCounts -p --countReadPairs -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH}.${gff}.gff3 \
                         -o ${COUNTS_DIR}/counts.gene.${gff}.tsv *.bam
done
cd ${COUNTS_DIR}
gzip *.tsv
gzip *.tsv.summary




