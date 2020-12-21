#


rscript=~/research/ecolokam/bin/ratioTest/calc.pval.mean.R

for d in `ls | grep -v sh | grep -v log | grep -v calc.pval.R`
do
    echo $d
    cd $d
    R --vanilla --slave --args . < $rscript
    cd -
done




