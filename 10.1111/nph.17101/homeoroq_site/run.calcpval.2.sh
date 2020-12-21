#


cd /home/jqsun/homeoroq_site/${d}
rscript=/home/jqsun/homeoroq_site/calc.pval.mean.R

for d in `ls | grep _ | grep  vs`
do
    echo $d
    cd $d
    R --vanilla --slave --args . < $rscript
    cd -
done




