#!/bin/bash
#PBS --group=g-mlbi
#PBS -q cq
#PBS -l gpunum_job=0
#PBS -l cpunum_job=1
#PBS -l elapstim_req=24:00:00
#PBS -N calcpval

# d=IR1_0418vs0502
# d=IR1_0502vs0516
# d=KT2_0418vs0502
# d=KT2_0502vs0516
# d=KT5_0418vs0502
# d=KT5_0502vs0516
d=IR1_0516vs0418
# d=KT2_0516vs0418
# d=KT5_0516vs0418



cd /home/jqsun/homeoroq_site/${d}
rscript=/home/jqsun/homeoroq_site/calc.pval.R

for i in {1..10}
do
    echo $i
    R --vanilla --slave --args . $i < $rscript
done




