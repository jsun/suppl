library(tidyverse)
library(ggsci)


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
    trainstats$dataset <- factor(trainstats$dataset, levels = c('W1', 'W2', 'F', 'W1F', 'W2F'))
    fig_acc <- ggplot(trainstats, aes(x = epoch, y = acc, color = mode, group = interaction(mode, run))) +
                    geom_line() + facet_grid(dataset ~ model) +
                    scale_color_npg()
    finalstats <- trainstats %>% filter(epoch == max(epoch)) %>%
                    group_by(mode, dataset, model) %>%
                    summarise(mean = mean(acc), sd = sd(acc))

    list(stats = list(train = trainstats, final = finalstats), fig = list(acc_history = fig_acc))
}



summarise_valid_stats <- function(dpath, model_type, valid_tag) {
    
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
    
    .figout <- function(valid_tag, model_type) {
        pattern <- paste0('.', valid_tag, '.', model_type, '.tsv')
        validstats <- NULL
        for (fpath in sort(list.files(dpath, pattern = pattern))) {
            print(fpath)
            validstats <- rbind(validstats, .sumvalid(fpath))
        }
        validstats_df <- validstats %>%
            dplyr::filter(model %in% c('VGGNet16', 'VGGNet19', 'ResNet18', 'ResNet152')) %>%
            dplyr::filter(topacc != 5)  %>%
            group_by(dataset, model, topacc) %>%
            summarise(mean = mean(acc), sd = sd(acc))
        validstats_df$model <- factor(validstats_df$model, levels = c('VGGNet16', 'VGGNet19', 'ResNet18', 'ResNet152'))
        validstats_df$dataset <- str_replace(validstats_df$dataset, 'g', '')
        validstats_df$dataset <- factor(validstats_df$dataset, levels = c('W1', 'W2', 'F', 'W1F', 'W2F'))
        validstats_df$topacc <- str_replace_all(validstats_df$topacc, '1', 'Top-1')
        validstats_df$topacc <- str_replace_all(validstats_df$topacc, '3', 'Top-3')
        
        validstats_fig <- ggplot(validstats_df) +
            geom_bar(aes(x = dataset, y = mean), stat = 'identity', fill = '#808180') +
            geom_errorbar(aes(x = dataset, ymax = mean + sd, ymin = mean - sd, width = 0.4)) +
            facet_grid(topacc ~ model) + ylim(0, 1.0) +
            ylab('accuracy') + xlab('dataset')
        
        list(data = validstats_df, fig = validstats_fig)
    }
    
    # calculate micro averages of accuracy
    stats_image   <- .figout(valid_tag, model_type)
    
    # calculate macro averages of accuracy and plots figures
    macroave <- microave <- NULL
    pattern <- paste0('.', valid_tag, '.', model_type, '.tsv')
    for (fpath in sort(list.files(dpath, pattern = pattern))) {
        print(fpath)
        xconfmat <- .sumvalid(fpath, TRUE)$confmat
        keys <- unlist(strsplit(fpath, '__'))
        keys[2] <- norm_model_name(keys[2])
        i <- as.integer(unlist(strsplit(keys[3], '\\.'))[1])
        
        microave <- rbind(microave, data.frame(dataset = keys[1], model = keys[2], topacc = 'Top-1'),
                                    acc = sum(diag(xconfmat)) / sum(xconfmat))
        if (!(paste(keys[c(1, 2)], collapse = '__') %in% names(macroave))) {
            macroave[[paste(keys[c(1, 2)], collapse = '__')]] <- matrix(NA, ncol = 10, nrow = nrow(xconfmat))
            rownames(macroave[[paste(keys[c(1, 2)], collapse = '__')]]) <- rownames(xconfmat)
        }
        macroave[[paste(keys[c(1, 2)], collapse = '__')]][, i] <- diag(xconfmat) / rowSums(xconfmat)
    }
    if (dpath == 'weights_species') {
        xconfmat <- macroave[['W2F__ResNet152']]
    } else {
        xconfmat <- macroave[['W2Fg__ResNet152']]
    }
    xconfdf <- data.frame(predicted = rownames(xconfmat), xconfmat) %>%
         pivot_longer(- predicted, names_to = 'try', values_to = 'precision')
    xconfdf$predicted <- factor(xconfdf$predicted, levels  = rownames(xconfmat))
    #xconfdfsum <- xconfdf %>% group_by(predicted) %>% summarise(precision = mean(precision[!is.nan(precision)]))
    #figconf <- ggplot() +
    #    geom_bar(aes(x = predicted, y = precision), data = xconfdfsum, stat = 'identity', fill = '#999999') +
    #    geom_jitter(aes(x = predicted, y = precision), data = xconfdf, width = 0.4, size = 0.5) +
    #    theme(axis.text = element_text(angle = 45, vjust=1, hjust=1))

    # calculate confusion matrix for ResNet152#1
    if (dpath == 'weights_species') {
        net_confmat <- .sumvalid(paste0('W2F__resnet152__1.', valid_tag, '.', model_type, '.tsv'), TRUE)
    } else {
        net_confmat <- .sumvalid(paste0('W2Fg__resnet152__1.', valid_tag, '.', model_type, '.tsv'), TRUE)
    }
    
    list(accstats = stats_image$data, accfig = stats_image$fig, conffig = xconfdf, confmat = net_confmat$confmat)
}


plot_confmat <- function(confmat, remove_zeros = FALSE) {
    if (remove_zeros) {
        allzero_col <- (colSums(confmat) == 0)
        allzero_row <- (rowSums(confmat) == 0)
        allzero_colrow <- (allzero_col & allzero_row)
        confmat_mini <- confmat[!allzero_colrow, !allzero_colrow]
    } else {
        confmat_mini <- confmat
    }
    confdf <- data.frame(predicted_label = rownames(confmat_mini), confmat_mini)
    rownames(confdf) <- NULL
    confdf <- confdf %>% pivot_longer(- predicted_label, names_to = 'true_label', values_to = 'value')
    confdf$predicted_label <- factor(confdf$predicted_label, levels = rownames(confmat_mini))
    confdf$true_label <- factor(confdf$true_label, levels = rownames(confmat_mini))
    confdf$value[confdf$value == 0] <- NA
    f <- ggplot(confdf, aes(x = predicted_label, y = true_label, fill = value)) +
        geom_tile() + coord_fixed() +
        theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_fill_gradientn(limits = c(1, 40), colours = c('#20854E', '#BC3C29'), na.value = 'white')
    
    f
}


plot_acc_boxplot <- function(confmat1, confmat2) {
    acc1 <- diag(confmat1) / rowSums(confmat1)
    acc2 <- diag(confmat2) / rowSums(confmat2)
    
}



# change TRUE to perform summarization
if (TRUE) {

# species identification
spmodel_train_stats <- summarise_train_stats('weights_species')
spmodel_valid_stats_image <- summarise_valid_stats('weights_species', 'image', 'T_valid')
spmodel_valid_stats_d50   <- summarise_valid_stats('weights_species', 'd50', 'T_valid')

valid_sum_table <- data.frame(dataset = spmodel_valid_stats_image$accstats$dataset,
                              model   = spmodel_valid_stats_image$accstats$model,
                              topacc  = spmodel_valid_stats_image$accstats$topacc,
                              image_mean = spmodel_valid_stats_image$accstats$mean,
                              image_sd   = spmodel_valid_stats_image$accstats$sd,
                              binded_mean = spmodel_valid_stats_d50$accstats$mean,
                              binded_sd   = spmodel_valid_stats_d50$accstats$sd)

valid_sum_table$dataset <- factor(valid_sum_table$dataset, levels = c('W1', 'W2', 'F', 'W1F', 'W2F'))
valid_sum_table$model <- factor(valid_sum_table$model, levels = c('VGGNet16', 'VGGNet19', 'ResNet18', 'ResNet152'))
valid_sum_table <- valid_sum_table %>%
        arrange(dataset, model)

write_tsv(valid_sum_table[valid_sum_table$topacc == 'Top-1', -3],
          path = 'eval_results/validacc_top1.tsv', col_names = TRUE)
write_tsv(valid_sum_table[valid_sum_table$topacc == 'Top-3', -3],
          path = 'eval_results/validacc_top3.tsv', col_names = TRUE)
write.table(spmodel_valid_stats_image$confmat, 
            file = 'eval_results/validconfmatrix_W2F_resnet152_1_image.tsv', sep = '\t')
write.table(spmodel_valid_stats_d50$confmat, 
            file = 'eval_results/validconfmatrix_W2F_resnet152_1_mesh.tsv', sep = '\t')


png('eval_results/validacc_barplot_image.png', 3800, 1800, res = 460)
print(spmodel_valid_stats_image$accfig)
dev.off()

png('eval_results/validacc_barplot_mesh.png', 3800, 1800, res = 460)
print(spmodel_valid_stats_d50$accfig)
dev.off()


png('eval_results/trainstats_species_history.png', 2200, 1800, res = 220)
print(spmodel_train_stats$fig$acc_history)
dev.off()
png('eval_results/trainstats_species_acc.png', 2200, 1800, res = 220)
print(spmodel_train_stats$fig$acc)
dev.off()


cf1 <- spmodel_valid_stats_image$confmat
cf2 <- spmodel_valid_stats_d50$confmat
keep1 <- (colSums(cf1) > 0) | (rowSums(cf1) > 0)
keep2 <- (colSums(cf2) > 0) | (rowSums(cf2) > 0)
keep <- (keep1 | keep2)

png('eval_results/spmodel_confmat_image.png', 2900, 2900, res = 220)
print(plot_confmat(cf1[keep, keep]))
dev.off()
png('eval_results/spmodel_confmat_d50.png', 2900, 2900, res = 220)
print(plot_confmat(cf2[keep, keep]))
dev.off()




confdf <- rbind(data.frame(spmodel_valid_stats_image$conffig, model = 'image-based'),
                data.frame(spmodel_valid_stats_d50$conffig, model = 'combined'))
confdfsum <- confdf %>% group_by(model, predicted) %>% summarise(precision = mean(precision[!is.nan(precision)]))
conffig <- ggplot() +
        geom_bar(aes(x = predicted, y = precision), data = confdfsum, stat = 'identity', fill = '#999999') +
        geom_jitter(aes(x = predicted, y = precision), data = confdf, width = 0.4, size = 0.5) +
        theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
        facet_grid(model ~ .)

png('eval_results/macro_precision_barplot.png', 7000, 1500, res=220)
print(conffig)
dev.off()




# genus identification
gnmodel_train_stats <- summarise_train_stats('weights_genus')
gnmodel_valid_stats_image <- summarise_valid_stats('weights_genus', 'image', 'Tg_valid')
gnmodel_valid_stats_d50   <- summarise_valid_stats('weights_genus', 'd50', 'Tg_valid')

valid_sum_table <- data.frame(dataset = gnmodel_valid_stats_image$accstats$dataset,
                              model   = gnmodel_valid_stats_image$accstats$model,
                              topacc  = gnmodel_valid_stats_image$accstats$topacc,
                              image_mean = gnmodel_valid_stats_image$accstats$mean,
                              image_sd   = gnmodel_valid_stats_image$accstats$sd,
                              binded_mean = gnmodel_valid_stats_d50$accstats$mean,
                              binded_sd   = gnmodel_valid_stats_d50$accstats$sd)

valid_sum_table$dataset <- factor(valid_sum_table$dataset, levels = c('W1', 'W2', 'F', 'W1F', 'W2F'))
valid_sum_table$model <- factor(valid_sum_table$model, levels = c('VGGNet16', 'VGGNet19', 'ResNet18', 'ResNet152'))
valid_sum_table <- valid_sum_table %>%
        arrange(dataset, model)

write_tsv(valid_sum_table[valid_sum_table$topacc == 'Top-1', -3],
          path = 'eval_results/validaccg_top1.tsv', col_names = TRUE)
write_tsv(valid_sum_table[valid_sum_table$topacc == 'Top-3', -3],
          path = 'eval_results/validaccg_top3.tsv', col_names = TRUE)
write.table(gnmodel_valid_stats_image$confmat, 
            file = 'eval_results/validconfmatrixg_W2F_resnet152_1_image.tsv', sep = '\t')
write.table(gnmodel_valid_stats_d50$confmat, 
            file = 'eval_results/validconfmatrixg_W2F_resnet152_1_mesh.tsv', sep = '\t')


png('eval_results/validaccg_barplot_image.png', 3800, 1800, res = 460)
print(gnmodel_valid_stats_image$accfig)
dev.off()

png('eval_results/trainstats_genus_history.png', 2200, 1800, res = 220)
print(gnmodel_train_stats$fig$acc_history)
dev.off()
png('eval_results/trainstats_genus_acc.png', 2200, 1800, res = 220)
print(gnmodel_train_stats$fig$acc)
dev.off()


cf1 <- gnmodel_valid_stats_image$confmat
cf2 <- gnmodel_valid_stats_d50$confmat
keep1 <- (colSums(cf1) > 0) | (rowSums(cf1) > 0)
keep2 <- (colSums(cf2) > 0) | (rowSums(cf2) > 0)
keep <- (keep1 | keep2)

png('eval_results/gnmodel_confmat_image.png', 2900, 2900, res = 220)
print(plot_confmat(cf1[keep, keep]))
dev.off()
png('eval_results/gnmodel_confmat_d50.png', 2900, 2900, res = 220)
print(plot_confmat(cf2[keep, keep]))
dev.off()







confdf <- rbind(data.frame(gnmodel_valid_stats_image$conffig, model = 'image-based'),
                data.frame(gnmodel_valid_stats_d50$conffig, model = 'combined'))
confdfsum <- confdf %>% group_by(model, predicted) %>% summarise(precision = mean(precision[!is.nan(precision)]))
conffig <- ggplot() +
        geom_bar(aes(x = predicted, y = precision), data = confdfsum, stat = 'identity', fill = '#999999') +
        geom_jitter(aes(x = predicted, y = precision), data = confdf, width = 0.4, size = 0.5) +
        theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
        facet_grid(model ~ .)

png('eval_results/macro_precision_barplot_genus.png', 2500, 1500, res=220)
print(conffig)
dev.off()






}







