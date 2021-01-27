library(tidyverse)
library(ggsci)

LV_MODEL <- c('MobileNetV2', 'VGGNet16', 'VGGNet19', 'ResNet18', 'ResNet152', 'DenseNet121')
LV_DATA <- c('W1', 'W2', 'F', 'W1F', 'W2F')

norm_model_name <- function(x) {
    if (x == 'vgg') {x = 'VGGNet16'}
    if (x == 'vgg19') {x = 'VGGNet19'}
    if (x == 'resnet') {x = 'ResNet18'}
    if (x == 'resnet152') {x = 'ResNet152'}
    if (x == 'densenet') {x = 'DenseNet121'}
    if (x == 'mobilenet') {x = 'MobileNetV2'}
    x
}

summarise_train_stats <- function(dpath) {
    trainstats <- NULL
    
    for (fpath in sort(list.files(dpath, pattern = '.train_hisotry.tsv'))) {
        fname <- str_split(fpath, '__', simplify = TRUE)[1, ]
        dataset_name <- fname[1]
        model_name   <- norm_model_name(fname[2])
        run_id       <- str_split(fname[3], '\\.', simplify = TRUE)[, 1]
        
        d <- read_tsv(paste0(dpath, '/', fpath), col_names = TRUE)
        trainstats <- bind_rows(trainstats,
                d %>%
                select(c('train_acc', 'val_acc')) %>%
                mutate(epoch = 1:length(train_acc)) %>%
                pivot_longer(c('train_acc', 'val_acc'), names_to = 'mode', values_to = 'acc') %>%
                mutate(mode = str_split(mode, '_', simplify = TRUE)[, 1],
                       dataset = dataset_name, model = model_name, run = run_id))
    }
    
    trainstats$dataset <- str_replace(trainstats$dataset, 'g', '')
    trainstats$dataset <- factor(trainstats$dataset, levels = LV_DATA)
    trainstats$model   <- factor(trainstats$model, levels = LV_MODEL)
    fig_acc <- ggplot(trainstats, aes(x = epoch, y = acc, color = mode, group = interaction(mode, run))) +
                    geom_line() + facet_grid(dataset ~ model) +
                    scale_color_npg()
    finalstats <- trainstats %>% filter(epoch == max(epoch)) %>%
                    group_by(mode, dataset, model) %>%
                    summarise(mean = mean(acc), sd = sd(acc))

    list(stats = list(train = trainstats, final = finalstats), fig = list(acc_history = fig_acc))
}



summarise_valid_stats <- function(dpath, model_type) {
    
    .sumvalid <- function(fpath, full_return = FALSE) {
        # validation condition
        fname <- str_split(fpath, '__', simplify = TRUE)[1, ]
        dataset_name <- fname[1]
        model_name   <- norm_model_name(fname[2])
        run_id       <- str_split(fname[3], '\\.', simplify = TRUE)[, 1]
        
        # validation results
        d <- read_tsv(paste0(dpath, '/', fpath), col_names = TRUE)
        true_labels <- str_split(d$X1, '/', simplify = TRUE)[, 9]
        probs <- as.matrix(d %>% select(-X1))
        pred_labels <- apply(probs, 1, function(i) {
            names(i)[i == max(i)][1]
        })
        
        # confusion matrix
        # vertical classes are true label, whereas horizental classes are predicted label
        if (full_return) {
            confmatrix <- matrix(0, ncol = length(colnames(probs)),
                                    nrow = length(colnames(probs)))
            colnames(confmatrix) <- rownames(confmatrix) <- colnames(probs)
            for (i in 1:length(true_labels)) {
                confmatrix[pred_labels[i], true_labels[i]] <- confmatrix[pred_labels[i], true_labels[i]] + 1
            }
            #keep_ij <- (colSums(confmatrix) > 0) | (rowSums(confmatrix) > 0)
            #confmatrix <- confmatrix[keep_ij, keep_ij]
            df <- data.frame(true = rownames(confmatrix), confmatrix) %>%
                        pivot_longer(-true, names_to = 'predict', values_to = 'counts')
            df$true <- factor(df$true, levels = colnames(probs))
            df$predict <- factor(df$predict, levels = colnames(probs))
            df$counts[df$counts == 0] <- NA
            hmap <- ggplot(df, aes(x = predict, y = true, fill = counts)) +
                geom_tile() + coord_fixed() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  strip.background = element_rect(fill = "white", colour = "white")) +
                scale_fill_continuous(na.value = 'white')
        }

        # calculate top-N accuracy
        acc <- NULL
        for (top_n in c(1, 3, 5)) {
            pred_labels <- apply(probs, 1, function(i) {
                paste(names(i)[rank(-i, ties.method = 'min') <= top_n], collapse = '__')
            })
            crr <- apply(cbind(true_labels, pred_labels), 1, function(i) {
                length(grep(i[1], i[2]))
            })
            acc <- c(acc, sum(crr) / length(crr))
        }
        
        stat_df <- data.frame(dataset = dataset_name, model = model_name, run = run_id,
                              topacc = c(1, 3, 5), acc = acc)
        
        if (full_return) {
            list(stat = stat_df, confmat = confmatrix, conffig = hmap)
        } else {
            stat_df
        }
    }
    
    .figout <- function(pattern) {
        validstats <- NULL
        for (fpath in sort(list.files(dpath, pattern = paste0('.', pattern, '.tsv')))) {
            print(fpath)
            validstats <- rbind(validstats, .sumvalid(fpath))
        }
        validstats_df <- validstats %>%
            group_by(dataset, model, topacc) %>%
            summarise(mean = mean(acc), sd = sd(acc))
        validstats_df$model <- factor(validstats_df$model, levels = LV_MODEL)
        validstats_df$dataset <- str_replace(validstats_df$dataset, 'g', '')
        validstats_df$dataset <- factor(validstats_df$dataset, levels = LV_DATA)
        validstats_df$topacc <- str_replace_all(validstats_df$topacc, '1', 'Top-1')
        validstats_df$topacc <- str_replace_all(validstats_df$topacc, '3', 'Top-3')
        validstats_df$topacc <- str_replace_all(validstats_df$topacc, '5', 'Top-5')
        
        validstats_fig <- ggplot(validstats_df) +
            geom_bar(aes(x = dataset, y = mean), stat = 'identity', fill = '#808180') +
            geom_errorbar(aes(x = dataset, ymax = mean + sd, ymin = mean - sd, width = 0.4)) +
            facet_grid(topacc ~ model) + ylim(0, 1.0) +
            ylab('accuracy') + xlab('dataset')
        
        list(data = validstats_df, fig = validstats_fig)
    }
    
    stats_image   <- .figout(model_type)
    if (dpath == 'weights_species') {
        vgg19_confmat <- .sumvalid(paste0('W2F__vgg19__1.T_valid.', model_type, '.tsv'), TRUE)
    } else {
        vgg19_confmat <- .sumvalid(paste0('W2Fg__vgg19__1.T_valid.', model_type, '.tsv'), TRUE)
    }
    
    list(accstats = stats_image$data, accfig = stats_image$fig, confmat = vgg19_confmat$confmat)
}





if (FALSE) {

# species identification
spmodel_train_stats <- summarise_train_stats('weights_species')
spmodel_valid_stats_image <- summarise_valid_stats('weights_species', 'image')
spmodel_valid_stats_mesh10  <- summarise_valid_stats('weights_species', 'd10')
spmodel_valid_stats_mesh20  <- summarise_valid_stats('weights_species', 'd20')
spmodel_valid_stats_mesh50  <- summarise_valid_stats('weights_species', 'd50')
spmodel_valid_stats_mesh100 <- summarise_valid_stats('weights_species', 'd100')

spmodel_valid_stats_image$accstats %>% ungroup() %>% filter(topacc == 'Top-1') %>% filter(mean == max(mean))
spmodel_valid_stats_mesh10$accstats %>% ungroup() %>% filter(topacc == 'Top-1') %>% filter(mean == max(mean))
spmodel_valid_stats_mesh20$accstats %>% ungroup() %>% filter(topacc == 'Top-1') %>% filter(mean == max(mean))
spmodel_valid_stats_mesh50$accstats %>% ungroup() %>% filter(topacc == 'Top-1') %>% filter(mean == max(mean))
spmodel_valid_stats_mesh100$accstats %>% ungroup() %>% filter(topacc == 'Top-1') %>% filter(mean == max(mean))

spmodel_valid_stats_image$accstats %>% ungroup() %>% filter(topacc == 'Top-3') %>% filter(mean == max(mean))
spmodel_valid_stats_mesh10$accstats %>% ungroup() %>% filter(topacc == 'Top-3') %>% filter(mean == max(mean))
spmodel_valid_stats_mesh20$accstats %>% ungroup() %>% filter(topacc == 'Top-3') %>% filter(mean == max(mean))
spmodel_valid_stats_mesh50$accstats %>% ungroup() %>% filter(topacc == 'Top-3') %>% filter(mean == max(mean))
spmodel_valid_stats_mesh100$accstats %>% ungroup() %>% filter(topacc == 'Top-3') %>% filter(mean == max(mean))


valid_sum_table <- data.frame(dataset = spmodel_valid_stats_image$accstats$dataset,
                              model   = spmodel_valid_stats_image$accstats$model,
                              topacc  = spmodel_valid_stats_image$accstats$topacc,
                              image_mean = spmodel_valid_stats_image$accstats$mean,
                              image_sd   = spmodel_valid_stats_image$accstats$sd,
                              binded_mean = spmodel_valid_stats_mesh50$accstats$mean,
                              binded_sd   = spmodel_valid_stats_mesh50$accstats$sd) %>%
                    arrange(topacc, dataset, model)
write_tsv(valid_sum_table,
          path = 'eval_results/validacc.tsv', col_names = TRUE)
write.table(spmodel_valid_stats_image$confmat, 
            file = 'eval_results/validconfmatrix_W2F_vgg19_1_image.tsv', sep = '\t')
write.table(spmodel_valid_stats_mesh50$confmat, 
            file = 'eval_results/validconfmatrix_W2F_vgg19_1_mesh.tsv', sep = '\t')


png('eval_results/validacc_barplot_image.png', 4400, 1800, res = 460)
spmodel_valid_stats_image$accfig
dev.off()

png('eval_results/validacc_barplot_mesh.png', 4400, 1800, res = 460)
spmodel_valid_stats_mesh50$accfig
dev.off()


png('eval_results/trainstats_species_history.png', 3200, 1800, res = 220)
print(spmodel_train_stats$fig$acc_history)
dev.off()
png('eval_results/trainstats_species_acc.png', 3200, 1800, res = 220)
print(spmodel_train_stats$fig$acc)
dev.off()



# genus identification
gnmodel_train_stats <- summarise_train_stats('weights_genus')
gnmodel_valid_stats_image <- summarise_valid_stats('weights_genus', 'image')
gnmodel_valid_stats_mesh50 <- summarise_valid_stats('weights_genus', 'd50')

valid_sum_table <- data.frame(dataset = gnmodel_valid_stats_image$accstats$dataset,
                              model   = gnmodel_valid_stats_image$accstats$model,
                              topacc  = gnmodel_valid_stats_image$accstats$topacc,
                              image_mean = gnmodel_valid_stats_image$accstats$mean,
                              image_sd   = gnmodel_valid_stats_image$accstats$sd,
                              binded_mean = gnmodel_valid_stats_mesh50$accstats$mean,
                              binded_sd   = gnmodel_valid_stats_mesh50$accstats$sd) %>%
                    arrange(topacc, dataset, model)

write_tsv(valid_sum_table,
          path = 'eval_results/validaccg.tsv', col_names = TRUE)
write.table(gnmodel_valid_stats_image$confmat, 
            file = 'eval_results/validconfmatrixg_W2F_vgg19_1_image.tsv', sep = '\t')
write.table(gnmodel_valid_stats_mesh50$confmat, 
            file = 'eval_results/validconfmatrixg_W2F_vgg19_1_mesh.tsv', sep = '\t')


png('eval_results/validaccg_barplot_image.png', 4400, 1800, res = 460)
gnmodel_valid_stats_image$accfig
dev.off()

png('eval_results/validaccg_barplot_mesh.png', 4400, 1800, res = 460)
gnmodel_valid_stats_mesh50$accfig
dev.off()


png('eval_results/trainstats_genus_history.png', 3200, 1800, res = 220)
print(gnmodel_train_stats$fig$acc_history)
dev.off()
png('eval_results/trainstats_genus_acc.png', 3200, 1800, res = 220)
print(gnmodel_train_stats$fig$acc)
dev.off()


}


