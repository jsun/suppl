# Data preparation

```
# R
source('scripts/format_data.R')




#pp <- function() {
#x <- read.table('test.result.percente.txt', sep = '\t', header = TRUE)
#x$prefecture <- factor(x$prefecture, levels = rev(pref.name.en))
#x$predicted <- x$predicted

#g <- ggplot(x, aes(x = month, y = prefecture, fill =  predicted, label = round(predicted, digits=1))) +
#    geom_tile() + geom_text() + theme_bw() + scale_fill_gradient(low = 'white', high = 'red')
#print(g)
#}


```




# Modeling


```
cd nn_model


python train.py --model L1 --output ./weights/test_modelarch.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 10 --batch-size 1024 --train-mode train


python train.py --model L1 --output ./weights/test_modelarch.pth \
                --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                --epochs 10 --batch-size 1024 --train-mode cv


# grid search
qsub train.sh

```


```
cd nn_model


python trainL2.py --model L2 --output ./weights/test_modelarch.pth \
                  --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                  --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                  --epochs 10 --batch-size 1024 --train-mode train


python trainL2.py --model L2 --output ./weights/test_modelarch.pth \
                  --train-dataset ../formatted_data/kyuuri_honpo_percent.train.tsv \
                  --valid-dataset ../formatted_data/kyuuri_honpo_percent.valid.tsv \
                  --epochs 10 --batch-size 1024 --train-mode cv


# grid search
qsub train.sh

```







