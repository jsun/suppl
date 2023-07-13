#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
# #PBS -l memsz_job=128gb
#PBS -l elapstim_req=72:00:00
#PBS -N HB_TAG_QUANT_EAGLENGI

nCPU=16

PROJECT_DIR=~/projects/HuoberBrezel
UTILS=~/local/utilfunc
EAGLE_SCRIPT_DIR=~/local/src/eagle/scripts

DATA_DIR=${PROJECT_DIR}/data/tagseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bameaglerc
COUNTS_DIR=${DATA_DIR}/countseaglrc
EAGLE_DIR=${DATA_DIR}/eaglerc


GENOME_DIR=~/projects/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX_PREFIX=${GENOME_DIR}/index/dna_
GENOME_INDEX_SUFFIX=_star
GTF_A=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.chrA.gff3
GTF_B=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.chrB.gff3
GTF_D=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.chrD.gff3



pSTAR=1
pEAGLERC=1
pCountHomeolog=1


echo "shell: 05_quant_eaglercngi.sh"
echo "pSTAR: ${pSTAR}"
echo "pEAGLERC: ${pEAGLERC}"
echo "pCount: ${pCountHomeolog}"
echo "nCPU: ${nCPU}"
echo "----- start -----"
echo `date`


cd ${PROJECT_DIR}


## mapping with STAR
if [ ${pSTAR} -eq 1 ]; then
    mkdir -p ${BAM_DIR}
    
    cd ${BAM_DIR}
    for fastq_filepath in `ls ${CLEAN_FASTQ_DIR}/*.clean.fastq.gz`; do
        echo ${fastq_filepath}

        fastq_filename=`basename ${fastq_filepath} .clean.fastq.gz`
        
        # mapping
        chr_names=(chrA chrB chrD)
        for chr_name in ${chr_names[@]}; do
            starsam=${BAM_DIR}/${fastq_filename}__on__${chr_name}
            mkdir -p ${starsam}
            cd ${starsam}
            
            echo "map ${fastq_filename} on ${chr_name} chromosome."
            echo `date`
            time STAR \
                --genomeDir ${GENOME_INDEX_PREFIX}${chr_name}${GENOME_INDEX_SUFFIX} \
                --readFilesCommand zcat \
                --readFilesIn ${fastq_filepath} \
                --runThreadN ${nCPU} --genomeLoad NoSharedMemory \
                --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout
            echo `date`
            echo "finished."
            
            echo "start samtools for sorting BAM."
            echo `date`
            samtools view -@ ${nCPU} -Shb Aligned.out.sam > Aligned.out.bam
            samtools sort -@ ${nCPU} -o Aligned.out.refsort.bam Aligned.out.bam
            samtools index -c Aligned.out.refsort.bam
            echo `date`
            echo "fnished."
            
            rm Aligned.out.sam
            cd -
        done
    done
fi





## EAGLERC (ngi) processes 
if [ ${pEAGLERC} -eq 1 ]; then
    cd ${BAM_DIR}
    for seqlib in `ls | awk 'BEGIN{FS="__"}{print $1}' | sort | uniq` ; do
        echo ${seqlib}
        mkdir -p ${EAGLE_DIR}/${seqlib}
        
        bam_on_chrA=${seqlib}__on__chrA/Aligned.out.refsort.bam
        bam_on_chrB=${seqlib}__on__chrB/Aligned.out.refsort.bam
        bam_on_chrD=${seqlib}__on__chrD/Aligned.out.refsort.bam
        
        time eagle-rc --ngi --listonly --splice --isc \
                    --ref1=${GENOME_DIR}/chrA.fa --ref2=${GENOME_DIR}/chrB.fa \
                    --bam1=${bam_on_chrA} --bam2=${bam_on_chrB} > ${EAGLE_DIR}/${seqlib}/eaglerc.out.A.vs.B.list
        time eagle-rc --ngi --listonly --splice --isc \
                    --ref1=${GENOME_DIR}/chrA.fa --ref2=${GENOME_DIR}/chrD.fa \
                    --bam1=${bam_on_chrA} --bam2=${bam_on_chrD} > ${EAGLE_DIR}/${seqlib}/eaglerc.out.A.vs.D.list
        time eagle-rc --ngi --listonly --splice --isc \
                    --ref1=${GENOME_DIR}/chrB.fa --ref2=${GENOME_DIR}/chrD.fa \
                    --bam1=${bam_on_chrB} --bam2=${bam_on_chrD} > ${EAGLE_DIR}/${seqlib}/eaglerc.out.B.vs.D.list
        
        python ${EAGLE_SCRIPT_DIR}/ref3_ngi_consensus.py -u -d \
                    -o ${EAGLE_DIR}/${seqlib}/${seqlib}.ref \
                    -AB ${EAGLE_DIR}/${seqlib}/eaglerc.out.A.vs.B.list \
                    -AD ${EAGLE_DIR}/${seqlib}/eaglerc.out.A.vs.D.list \
                    -BD ${EAGLE_DIR}/${seqlib}/eaglerc.out.B.vs.D.list
    done
fi




# EAGLE-RC & featureCounts
cd ${PROJECT_DIR}
pyenv local eagle
mkdir -p ${COUNTS_DIR}

if [ ${pCountHomeolog} -eq 1 ]; then
    mkdir -p ${EAGLE_DIR}
    cd ${EAGLE_DIR}
    
    for seqlib in `ls`; do
        cd ${seqlib}
        
        bam_on__A=${BAM_DIR}/${seqlib}__on__chrA
        bam_on__B=${BAM_DIR}/${seqlib}__on__chrB
        bam_on__D=${BAM_DIR}/${seqlib}__on__chrD
        
        # homeolog reads
        time eagle-rc --refonly --readlist -a ${bam_on__A}/Aligned.out.refsort.bam -o ${seqlib}.chrA.homeolog ${seqlib}.ref.chrA.list
        time eagle-rc --refonly --readlist -a ${bam_on__B}/Aligned.out.refsort.bam -o ${seqlib}.chrB.homeolog ${seqlib}.ref.chrB.list
        time eagle-rc --refonly --readlist -a ${bam_on__D}/Aligned.out.refsort.bam -o ${seqlib}.chrD.homeolog ${seqlib}.ref.chrD.list
        
        featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.counts.homeolog.txt ${seqlib}.chrA.homeolog.ref.bam
        featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.counts.homeolog.txt ${seqlib}.chrB.homeolog.ref.bam
        featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.counts.homeolog.txt ${seqlib}.chrD.homeolog.ref.bam
        mv ${seqlib}.chr*.counts.homeolog.* ${COUNTS_DIR}
        
        # subgenome specific reads
        echo "" > dummy.txt
        time eagle-rc --refonly --readlist -a ${bam_on__A}/Aligned.out.refsort.bam \
                                  -u ${bam_on__B}/Aligned.out.refsort.bam,${bam_on__D}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrA.specific dummy.txt
        time eagle-rc --refonly --readlist -a ${bam_on__B}/Aligned.out.refsort.bam \
                                  -u ${bam_on__A}/Aligned.out.refsort.bam,${bam_on__D}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrB.specific dummy.txt
        time eagle-rc --refonly --readlist -a ${bam_on__D}/Aligned.out.refsort.bam \
                                  -u ${bam_on__A}/Aligned.out.refsort.bam,${bam_on__B}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrD.specific dummy.txt
        
        featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.counts.specific.txt ${seqlib}.chrA.specific.ref.bam
        featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.counts.specific.txt ${seqlib}.chrB.specific.ref.bam
        featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.counts.specific.txt ${seqlib}.chrD.specific.ref.bam
        mv ${seqlib}.chr*.counts.specific.* ${COUNTS_DIR}
        
        cd -
    done
    
fi




echo `date`
echo "----- finish -----"
