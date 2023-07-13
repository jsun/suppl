#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l memsz_job=128gb
#PBS -l elapstim_req=72:00:00
#PBS -N HB_FUL_QUANT_EAGLE_4

nCPU=16

PROJECT_DIR=~/projects/HuoberBrezel

EAGLE_SCRIPT_DIR=~/local/src/eagle/scripts

DATA_DIR=${PROJECT_DIR}/data/fullseq
CLEAN_FASTQ_DIR=${DATA_DIR}/clean_fastq
BAM_DIR=${DATA_DIR}/bameaglerc
COUNTS_DIR=${DATA_DIR}/countseaglrc
EAGLE_DIR=${DATA_DIR}/eaglerc

GENOME_DIR=~/projecs/data/genome/IWGSC_RefSeq_v1.1_CS
GENOME_INDEX_PREFIX=${GENOME_DIR}/index/dna_
GENOME_INDEX_SUFFIX=_star
GTF_A=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.chrA.gff3
GTF_B=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.chrB.gff3
GTF_D=${GENOME_DIR}/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.ext_1.0k.chrD.gff3


pSTAR=0
pEAGLE=0
pEAGLERC=0
pCountHomeolog=1


echo "shell: 02_fullseq.eaglerc.sh"
echo "LIBPTN: ${LIBPTN}"
echo "pSTAR: ${pSTAR}"
echo "pEAGLE: ${pEAGLE}"
echo "pEAGLERC: ${pEAGLERC}"
echo "pCountHomeolog: ${pCountHomeolog}"
echo "pCountUnique: ${pCountUnique}"
echo "nCPU: ${nCPU}"
echo "----- start -----"
echo `date`


cd ${PROJECT_DIR}


## mapping with STAR
if [ ${pSTAR} -eq 1 ]; then

    mkdir -p ${BAM_DIR}
    
    cd ${BAM_DIR}
    for fastq_filepath in `ls ${CLEAN_FASTQ_DIR}/*_R1.cleaned.fastq.gz`; do
        echo ${fastq_filepath}

        fastq_filename=`basename ${fastq_filepath} _R1.cleaned.fastq.gz`
        
        # mapping
        chr_names=(chrA chrB chrD)
        for chr_name in ${chr_names[@]}; do
            starsam=${BAM_DIR}/${fastq_filename}__on__${chr_name}
            mkdir -p ${starsam}
            cd ${starsam}
            
            echo "map ${fastq_filename} on ${chr_name} chromosome."
            time STAR \
                --genomeDir ${GENOME_INDEX_PREFIX}${chr_name}${GENOME_INDEX_SUFFIX} \
                --readFilesCommand zcat \
                --readFilesIn ${CLEAN_FASTQ_DIR}/${fastq_filename}_R1.cleaned.fastq.gz ${CLEAN_FASTQ_DIR}/${fastq_filename}_R2.cleaned.fastq.gz \
                --runThreadN ${nCPU} --genomeLoad NoSharedMemory \
                --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                --outSJfilterCountUniqueMin 3 2 2 2 --outMultimapperOrder Random --outFilterType BySJout
            echo "finished."
            
            echo "start bamsort."
            samtools view -@ ${nCPU} -Shb Aligned.out.sam > Aligned.out.bam
            samtools sort -@ ${nCPU} -o Aligned.out.refsort.bam Aligned.out.bam
            samtools index -c Aligned.out.refsort.bam
            echo "fnished bamsort."
            
            rm Aligned.out.sam
            cd -
        done
        
    done
fi





## EAGLE processes (multi-CPUs)
if [ ${pEAGLE} -eq 1 ]; then
    cd ${BAM_DIR}
    for seqlib in `ls | awk 'BEGIN{FS="__on__"}{print $1}' | sort | uniq`; do
        echo ${seqlib}
        chr_names=(A B D)
        for chr_name in ${chr_names[@]}; do
            cd ${seqlib}__on__chr${chr_name}
            for vcf_file in `ls ${GENOME_DIR}/${chr_name}.vs.*.gtf.vcf`; do
                vcf_basename=`basename ${vcf_file} .gtf.vcf`
                echo "start eagle"
                time eagle -t ${nCPU} -a Aligned.out.refsort.bam \
                     -r ${GENOME_DIR}/chr${chr_name}.fa \
                     -v ${vcf_file} \
                     --splice --rc 1> eagle.out.${vcf_basename}.txt 2> eagle.out.${vcf_basename}.readinfo.txt
                echo "finished eagle"
            done
            cd -
        done
    done
fi





## EAGLE-RC processes (single-CPUs)
if [ ${pEAGLERC} -eq 1 ]; then
    cd ${BAM_DIR}
    for seqlib in `ls | awk 'BEGIN{FS="__on__"}{print $1}' | sort | uniq`; do
        echo ${seqlib}
        chr_names=(A B D)
        for chr_name in ${chr_names[@]}; do
            cd ${seqlib}__on__chr${chr_name}
            for vcf_file in `ls ${GENOME_DIR}/${chr_name}.vs.*.gtf.vcf`; do
                vcf_basename=`basename ${vcf_file} .gtf.vcf`
                echo "start eaglerc"
                time eagle-rc --listonly -a Aligned.out.refsort.bam \
                                -o eaglerc.${vcf_basename} \
                                -v eagle.out.${vcf_basename}.txt eagle.out.${vcf_basename}.readinfo.txt > eaglerc.out.${vcf_basename}.list
                echo "finished eaglerc"
            done
            cd ../
        done
    done
fi






# EAGLE-RC & featureCounts
cd ${PROJECT_DIR}
pyenv local eagle
mkdir -p ${COUNTS_DIR}

if [ ${pCountHomeolog} -eq 1 ]; then
    mkdir -p ${EAGLE_DIR}
    cd ${EAGLE_DIR}
    
    for seqlib in `ls  ${BAM_DIR} | awk 'BEGIN{FS="__on__"}{print $1}' | sort | uniq`; do
        bam_on__A=${BAM_DIR}/${seqlib}__on__chrA
        bam_on__B=${BAM_DIR}/${seqlib}__on__chrB
        bam_on__D=${BAM_DIR}/${seqlib}__on__chrD

        mkdir -p ${seqlib}
        cd ${seqlib}
        
        # homeolog reads
        python ${EAGLE_SCRIPT_DIR}/ref3_consensus.py --pe -d -u -o ${seqlib}.ref        \
           -A ${bam_on__A}/eaglerc.out.A.vs.B.list ${bam_on__A}/eaglerc.out.A.vs.D.list \
           -B ${bam_on__B}/eaglerc.out.B.vs.A.list ${bam_on__B}/eaglerc.out.B.vs.D.list \
           -D ${bam_on__D}/eaglerc.out.D.vs.A.list ${bam_on__D}/eaglerc.out.D.vs.B.list
        
        time eagle-rc --refonly --paired --readlist -a ${bam_on__A}/Aligned.out.refsort.bam -o ${seqlib}.chrA.homeolog ${seqlib}.ref.chrA.list
        time eagle-rc --refonly --paired --readlist -a ${bam_on__B}/Aligned.out.refsort.bam -o ${seqlib}.chrB.homeolog ${seqlib}.ref.chrB.list
        time eagle-rc --refonly --paired --readlist -a ${bam_on__D}/Aligned.out.refsort.bam -o ${seqlib}.chrD.homeolog ${seqlib}.ref.chrD.list
        
        featureCounts -p -T 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.counts.homeolog.txt ${seqlib}.chrA.homeolog.ref.bam
        featureCounts -p -T 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.counts.homeolog.txt ${seqlib}.chrB.homeolog.ref.bam
        featureCounts -p -T 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.counts.homeolog.txt ${seqlib}.chrD.homeolog.ref.bam
        mv ${seqlib}.chr*.counts.homeolog.* ${COUNTS_DIR}
        
        
        # subgenome specific reads
        echo "" > dummy.txt
        time eagle-rc --refonly --paired --readlist -a ${bam_on__A}/Aligned.out.refsort.bam \
                                  -u ${bam_on__B}/Aligned.out.refsort.bam,${bam_on__D}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrA.specific dummy.txt
        time eagle-rc --refonly --paired --readlist -a ${bam_on__B}/Aligned.out.refsort.bam \
                                  -u ${bam_on__A}/Aligned.out.refsort.bam,${bam_on__D}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrB.specific dummy.txt
        time eagle-rc --refonly --paired --readlist -a ${bam_on__D}/Aligned.out.refsort.bam \
                                  -u ${bam_on__A}/Aligned.out.refsort.bam,${bam_on__B}/Aligned.out.refsort.bam \
                                  -o ${seqlib}.chrD.specific dummy.txt
        
        featureCounts -p -T 1 -t gene -g ID -a ${GTF_A} -o ${seqlib}.chrA.counts.specific.txt ${seqlib}.chrA.specific.ref.bam
        featureCounts -p -T 1 -t gene -g ID -a ${GTF_B} -o ${seqlib}.chrB.counts.specific.txt ${seqlib}.chrB.specific.ref.bam
        featureCounts -p -T 1 -t gene -g ID -a ${GTF_D} -o ${seqlib}.chrD.counts.specific.txt ${seqlib}.chrD.specific.ref.bam
        mv ${seqlib}.chr*.counts.specific.* ${COUNTS_DIR}
        
        cd -
    done
fi



echo `date`
echo "----- finish -----"
