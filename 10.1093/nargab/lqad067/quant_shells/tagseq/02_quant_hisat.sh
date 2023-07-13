#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=16
#PBS -l memsz_job=64gb
#PBS -l elapstim_req=72:00:00
#PBS -N HB_TAG_QUANT_COUNT


nCPU=16
PROJECT_DIR=~/projects/HuoberBrezel
DATA_DIR=${PROJECT_DIR}/data/tagseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bam
COUNTS_DIR=${DATA_DIR}/counts
GENOME_DIR=~/projects/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX=${GENOME_DIR}/index/dna_hisat2
GTF_FILEPATH=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706



cd ${PROJECT_DIR}


# mapping reads with HISAT2
mkdir -p ${BAM_DIR}
cd ${CLEAN_FASTQ_DIR}
for fastq_filepath in `ls *.clean.fastq.gz`; do
    echo ${fastq_filepath}
    fastq_filename=`basename ${fastq_filepath} .clean.fastq.gz`

    echo 'start mapping ...'
    echo `date`
    time hisat2 -p ${nCPU} --no-spliced-alignment -x ${GENOME_INDEX} \
                       -U ${fastq_filepath} -S ${BAM_DIR}/${fastq_filename}.sam
    echo `date`
    echo 'finished.'

    cd ${BAM_DIR}
    samtools sort -@ ${nCPU} -O bam -o ${fastq_filename}.bam ${fastq_filename}.sam
    samtools index -c ${fastq_filename}.bam
        
    #rm ${fastq_filename}.sam
    cd -
done



# count mapped reads
mkdir -p ${COUNTS_DIR}
cd ${BAM_DIR}
gff_versions=("iwgsc" "ext_0.5k" "ext_1.0k" "ext_1.5k" "ext_2.0k" "ext_2.5k" "ext_3.0k" "ext_3.5k" "ext_4.0k")
for gff in ${gff_versions[@]}; do
    featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH}.${gff}.gff3 -s 1 \
                         -o ${COUNTS_DIR}/tcs.counts.gene.${gff}.tsv TaeRS2728_*.bam
    featureCounts -T ${nCPU} -t gene -g ID -a ${GTF_FILEPATH}.${gff}.gff3 -s 1 \
                         -o ${COUNTS_DIR}/cs.counts.gene.${gff}.tsv 20181109*.bam
done
cd ${COUNTS_DIR}
gzip *.tsv
gzip *.tsv.summary



