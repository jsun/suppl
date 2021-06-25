#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_hb_tag_eaglercngi_count
#$ -mods l_hard h_rt 720:00:00


##### ##### ##### ##### ##### #####
# RUN 04_quant_eaglerc.sh FIRST!  #
##### ##### ##### ##### ##### #####

pEAGLE=0
pCount=1


nCPU=16

PROJECT_DIR=~/projects/HuoberBrezel

BIN=~/local/bin
EAGLE_SCRIPT_DIR=~/local/src/eagle/scripts
PY2=~/.pyenv/versions/eagle/bin/python

DATA_DIR=${PROJECT_DIR}/data/tagseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bameaglercngi
COUNTS_DIR=${DATA_DIR}/countseaglercngi
EAGLE_DIR=${DATA_DIR}/eaglercngi

GENOME_DIR=~/projects/db/genome/IWGSC_RefSeq_v2.1_CS
GENOME_INDEX_PREFIX=${GENOME_DIR}/index/dna_
GENOME_INDEX_SUFFIX=_star

GTF_A=${GENOME_DIR}/seqdata/dna.ext_1.0k.chrA.gff3
GTF_B=${GENOME_DIR}/seqdata/dna.ext_1.0k.chrB.gff3
GTF_D=${GENOME_DIR}/seqdata/dna.ext_1.0k.chrD.gff3



echo "shell: 05_tagseq.eaglerc.sh"
echo "pEAGLE: ${pEAGLE}"
echo "pCount: ${pCount}"
echo "nCPU: ${nCPU}"
echo "----- start -----"
echo `date`


cd ${PROJECT_DIR}


if [ ${pEAGLE} -eq 1 ]; then

    cd ${DATA_DIR}
    cp -r bameaglerc ${BAM_DIR}
    rm ${BAM_DIR}/*/eaglerc*
    
    cd ${BAM_DIR}
    for seqlib in `ls | awk 'BEGIN{FS="__"}{print $1}' | sort | uniq` ; do
        echo ${seqlib}
        mkdir -p ${EAGLE_DIR}/${seqlib}
        
        bam_on_chrA=${seqlib}__on__chrA/Aligned.out.refsort.bam
        bam_on_chrB=${seqlib}__on__chrB/Aligned.out.refsort.bam
        bam_on_chrD=${seqlib}__on__chrD/Aligned.out.refsort.bam
        
        time ${BIN}/eagle-rc --ngi --listonly --splice --isc \
                    --ref1=${GENOME_DIR}/seqdata/dna_chrA.fa --ref2=${GENOME_DIR}/seqdata/dna_chrB.fa \
                    --bam1=${bam_on_chrA} --bam2=${bam_on_chrB} > ${EAGLE_DIR}/${seqlib}/eaglerc.out.A.vs.B.list
        time ${BIN}/eagle-rc --ngi --listonly --splice --isc \
                    --ref1=${GENOME_DIR}/seqdata/dna_chrA.fa --ref2=${GENOME_DIR}/seqdata/dna_chrD.fa \
                    --bam1=${bam_on_chrA} --bam2=${bam_on_chrD} > ${EAGLE_DIR}/${seqlib}/eaglerc.out.A.vs.D.list
        time ${BIN}/eagle-rc --ngi --listonly --splice --isc \
                    --ref1=${GENOME_DIR}/seqdata/dna_chrB.fa --ref2=${GENOME_DIR}/seqdata/dna_chrD.fa \
                    --bam1=${bam_on_chrB} --bam2=${bam_on_chrD} > ${EAGLE_DIR}/${seqlib}/eaglerc.out.B.vs.D.list
        
        ${PY2} ${EAGLE_SCRIPT_DIR}/ref3_ngi_consensus.py -u -d \
                    -o ${EAGLE_DIR}/${seqlib}/${seqlib}.ref \
                    -AB ${EAGLE_DIR}/${seqlib}/eaglerc.out.A.vs.B.list \
                    -AD ${EAGLE_DIR}/${seqlib}/eaglerc.out.A.vs.D.list \
                    -BD ${EAGLE_DIR}/${seqlib}/eaglerc.out.B.vs.D.list
    done
fi




# EAGLE-RC & featureCounts
cd ${PROJECT_DIR}
mkdir -p ${COUNTS_DIR}

if [ ${pCount} -eq 1 ]; then
    mkdir -p ${EAGLE_DIR}
    cd ${EAGLE_DIR}
    
    for seqlib in `ls`; do
        echo ${seqlib}
        cd ${seqlib}
        
        bam_on__A=${BAM_DIR}/${seqlib}__on__chrA
        bam_on__B=${BAM_DIR}/${seqlib}__on__chrB
        bam_on__D=${BAM_DIR}/${seqlib}__on__chrD
        
        # homeolog reads
        echo "count homeolog reads"
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__A}/Aligned.out.refsort.bam -o ${seqlib}.chrA.homeolog ${seqlib}.ref.chrA.list
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__B}/Aligned.out.refsort.bam -o ${seqlib}.chrB.homeolog ${seqlib}.ref.chrB.list
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__D}/Aligned.out.refsort.bam -o ${seqlib}.chrD.homeolog ${seqlib}.ref.chrD.list
        
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.counts.homeolog.txt ${seqlib}.chrA.homeolog.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.counts.homeolog.txt ${seqlib}.chrB.homeolog.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.counts.homeolog.txt ${seqlib}.chrD.homeolog.ref.bam
        mv ${seqlib}.chr*.counts.homeolog.* ${COUNTS_DIR}
        
        # subgenome specific reads
        echo "count subgenome specific reads"
        echo "" > dummy.txt
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__A}/Aligned.out.refsort.bam \
                                  -u ${bam_on__B}/Aligned.out.refsort.bam,${bam_on__D}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrA.specific dummy.txt
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__B}/Aligned.out.refsort.bam \
                                  -u ${bam_on__A}/Aligned.out.refsort.bam,${bam_on__D}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrB.specific dummy.txt
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__D}/Aligned.out.refsort.bam \
                                  -u ${bam_on__A}/Aligned.out.refsort.bam,${bam_on__B}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrD.specific dummy.txt
        
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.counts.specific.txt ${seqlib}.chrA.specific.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.counts.specific.txt ${seqlib}.chrB.specific.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.counts.specific.txt ${seqlib}.chrD.specific.ref.bam
        mv ${seqlib}.chr*.counts.specific.* ${COUNTS_DIR}
        
        cd -
    done
    
fi




echo `date`
echo "----- finish -----"
