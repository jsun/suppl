library(tidyverse)
library(ggsci)

set.seed(2020)
for (tag in c('genus',  'species')) {

x <- read_tsv(paste0(tag, '.log'), col_names = TRUE)
x <- x[x$model != 'mobilenet',]
x <- x[x$model != 'densenet',]
x$data <- factor(x$data, levels = c('W1', 'W2', 'F', 'W1F', 'W2F'))
x$model <- str_replace(x$model, 'vgg19', 'VGG19')
x$model <- str_replace(x$model, 'vgg', 'VGG16')
x$model <- str_replace(x$model, 'resnet152', 'ResNet152')
x$model <- str_replace(x$model, 'resnet', 'ResNet18')
x$model <- factor(x$model, levels = c('VGG16', 'VGG19', 'ResNet18', 'ResNet152'))
x_mu <- x %>%
        group_by(data, model) %>%
        summarise(mean = mean(time))
x_sd <- x %>%
        group_by(data, model) %>%
        summarise(sd = sd(time))

x_mu$sd <- x_sd$sd

g <- ggplot() +
    geom_jitter(aes(x = data, y = time), alpha=0.5, data = x) +
    facet_grid(~ model) +
    ylab('training time [h]')

png(paste0(tag, '.png'), 2000, 800, res = 300)
print(g)
dev.off()

}

