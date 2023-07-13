#!/bin/bash
#$ -S /bin/bash
#$ -jc hostos_g1
#$ -cwd
#$ -N qslog_hb_tag_eaglerc_counts
#$ -mods l_hard h_rt 720:00:00

pSTAR=0
pEAGLE=0
pCount=1



nCPU=16

PROJECT_DIR=~/projects/HuoberBrezel

BIN=~/local/bin
EAGLE_SCRIPT_DIR=~/local/src/eagle/scripts
PY2=~/.pyenv/versions/eagle/bin/python


DATA_DIR=${PROJECT_DIR}/data/tagseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bameaglerc
COUNTS_DIR=${DATA_DIR}/countseaglerc
EAGLE_DIR=${DATA_DIR}/eaglerc

GENOME_DIR=~/projects/db/genome/IWGSC_RefSeq_v2.1_CS
GENOME_INDEX_PREFIX=${GENOME_DIR}/index/dna_
GENOME_INDEX_SUFFIX=_star

GTF_A=${GENOME_DIR}/seqdata/dna.ext_1.0k.chrA.gff3
GTF_B=${GENOME_DIR}/seqdata/dna.ext_1.0k.chrB.gff3
GTF_D=${GENOME_DIR}/seqdata/dna.ext_1.0k.chrD.gff3



echo "shell: 04_quant_eaglerc.sh"
echo "pSTAR: ${pSTAR}"
echo "pEAGLE: ${pEAGLE}"
echo "pCount: ${pCount}"
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
            time ${BIN}/STAR \
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
            ${BIN}/samtools view -@ ${nCPU} -Shb Aligned.out.sam > Aligned.out.bam
            ${BIN}/samtools sort -@ ${nCPU} -o Aligned.out.refsort.bam Aligned.out.bam
            ${BIN}/samtools index -c Aligned.out.refsort.bam
            echo `date`
            echo "fnished."
            
            rm Aligned.out.sam
            cd -
        done
    done
fi





## EAGLE processes (multi-CPUs)
if [ ${pEAGLE} -eq 1 ]; then
    cd ${BAM_DIR}
    for seqlib in `ls | awk 'BEGIN{FS="__"}{print $1}' | sort | uniq` ; do
        echo ${seqlib}
        chr_names=(A B D)
        for chr_name in ${chr_names[@]}; do
            echo ${chr_name}
            cd ${seqlib}__on__chr${chr_name}
            for vcf_file in `ls ${GENOME_DIR}/seqdata/${chr_name}.vs.*.gtf.vcf`; do
                vcf_basename=`basename ${vcf_file} .gtf.vcf`
                echo "start eagle"
                echo `date`
                time ${BIN}/eagle -t ${nCPU} -a Aligned.out.refsort.bam \
                     -r ${GENOME_DIR}/seqdata/dna_chr${chr_name}.fa \
                     -v ${vcf_file} \
                     --splice --rc 1> eagle.out.${vcf_basename}.txt 2> eagle.out.${vcf_basename}.readinfo.txt
                echo `date`
                echo "finished eagle"
                echo "start eaglerc"
                echo `date`
                time ${BIN}/eagle-rc --listonly -a Aligned.out.refsort.bam \
                                -o eaglerc.${vcf_basename} \
                                -v eagle.out.${vcf_basename}.txt eagle.out.${vcf_basename}.readinfo.txt > eaglerc.out.${vcf_basename}.list
                echo `date`
                echo "finished eaglerc"
            done
            cd -
        done
    done
fi








# EAGLE-RC & featureCounts
cd ${PROJECT_DIR}

mkdir -p ${COUNTS_DIR}

if [ ${pCount} -eq 1 ]; then
    mkdir -p ${EAGLE_DIR}
    cd ${EAGLE_DIR}
    
    for seqlib in `ls  ${BAM_DIR} | awk 'BEGIN{FS="__"}{print $1}' | sort | uniq`; do
        bam_on__A=${BAM_DIR}/${seqlib}__on__chrA
        bam_on__B=${BAM_DIR}/${seqlib}__on__chrB
        bam_on__D=${BAM_DIR}/${seqlib}__on__chrD

        mkdir -p ${seqlib}
        cd ${seqlib}
        
        # homeolog reads
        echo "couting homeolog reads"
        ${PY2} ${EAGLE_SCRIPT_DIR}/ref3_consensus.py -d -u -o ${seqlib}.ref        \
           -A ${bam_on__A}/eaglerc.out.A.vs.B.list ${bam_on__A}/eaglerc.out.A.vs.D.list \
           -B ${bam_on__B}/eaglerc.out.B.vs.A.list ${bam_on__B}/eaglerc.out.B.vs.D.list \
           -D ${bam_on__D}/eaglerc.out.D.vs.A.list ${bam_on__D}/eaglerc.out.D.vs.B.list

        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__A}/Aligned.out.refsort.bam -o ${seqlib}.chrA.homeolog ${seqlib}.ref.chrA.list
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__B}/Aligned.out.refsort.bam -o ${seqlib}.chrB.homeolog ${seqlib}.ref.chrB.list
        time ${BIN}/eagle-rc --refonly --readlist -a ${bam_on__D}/Aligned.out.refsort.bam -o ${seqlib}.chrD.homeolog ${seqlib}.ref.chrD.list
        
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.counts.homeolog.txt ${seqlib}.chrA.homeolog.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.counts.homeolog.txt ${seqlib}.chrB.homeolog.ref.bam
        ${BIN}/featureCounts -T 1 -s 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.counts.homeolog.txt ${seqlib}.chrD.homeolog.ref.bam
        mv ${seqlib}.chr*.counts.homeolog.* ${COUNTS_DIR}
        
        
        # subgenome specific reads
        echo "couting subgenome specific reads"
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
