#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=24:00:00
#PBS -b 1
#PBS -m be
#PBS -N calcpval

#d=0418_KT2_IR1
#d=0418_KT5_IR1
#d=0418_KT5_KT2
#d=0516_KT2_IR1
#d=0516_KT5_KT2
#d=0516_KT5_IR1
#d=0502_KT5_IR1
#d=0502_KT2_IR1
d=0502_KT5_KT2


cd /home/jqsun/homeoroq_tmp/${d}
rscript=/home/jqsun/homeoroq_tmp/calc.pval.R

for i in {1..10}
do
    echo $i
    R --vanilla --slave --args . $i < $rscript
done




