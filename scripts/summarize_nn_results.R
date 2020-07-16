library(tidyverse)



nn_result_prefix <- 'model/weights'
nn_result_permu <- 'model/permutest'

cv_L1 <- read_tsv(paste0(nn_result_prefix, '/L1.pth_cv_stats.tsv'), col_names = TRUE)

L1p1 <- cv_L1 %>%
    filter(activate_fun == 'relu') %>%
    select(- valid_mse) %>%
    group_by(n_hidden_0, dropout) %>%
    summarize(mse = mean(train_mse)) %>%
    ungroup() %>%
    ggplot(aes(y = n_hidden_0, x = dropout)) +
        geom_tile(aes(fill = mse)) +
        geom_text(aes(label = round(mse, 2))) +
        scale_fill_gradient(low = 'white', high = 'steelblue')

png('~/Desktop/L1p1.png', 1300, 1300, res=250)
print(L1p1)
dev.off()


   

L1p2 <- cv_L1 %>%
    filter(activate_fun == 'relu') %>%
    select(- train_mse) %>%
    group_by(n_hidden_0, dropout) %>%
    summarize(mse = mean(valid_mse)) %>%
    ungroup() %>%
    ggplot(aes(y = n_hidden_0, x = dropout)) +
        geom_tile(aes(fill = mse)) +
        geom_text(aes(label = round(mse, 2))) +
        scale_fill_gradient(low = 'white', high = 'steelblue')

png('~/Desktop/L1p2.png', 1300, 1300, res=250)
print(L1p2)
dev.off()







cv_L2 <- read_tsv(paste0(nn_result_prefix, '/L2.pth_cv_stats.tsv'), col_names = TRUE)


L2p1 <- cv_L2 %>%
    filter(activate_fun == 'relu') %>%
    select(- valid_mse) %>%
    group_by(n_hidden_0, n_hidden_1, dropout) %>%
    summarize(mse = mean(train_mse)) %>%
    ungroup() %>%
    ggplot(aes(x = n_hidden_0, y = n_hidden_1)) +
        geom_tile(aes(fill = mse)) +
        geom_text(aes(label = round(mse, 2))) +
        scale_fill_gradient(low = 'white', high = 'steelblue') +
        facet_grid(. ~ dropout) 
png('~/Desktop/L2p1.png', 2400, 1300, res=250)
print(L2p1)
dev.off()

L2p2 <- cv_L2 %>%
    filter(activate_fun == 'relu') %>%
    select(- train_mse) %>%
    group_by(n_hidden_0, n_hidden_1, dropout) %>%
    summarize(mse = mean(valid_mse)) %>% #min, mean
    ungroup() %>%
    ggplot(aes(x = n_hidden_0, y = n_hidden_1)) +
        geom_tile(aes(fill = mse)) +
        geom_text(aes(label = round(mse, 2))) +
        scale_fill_gradient(low = 'white', high = 'steelblue') +
        facet_grid(. ~ dropout) 
png('~/Desktop/L2p2.png', 2400, 1300, res=250)
print(L2p2)
dev.off()






cv_L3 <- read_tsv(paste0(nn_result_prefix, '/L3.pth_cv_stats.tsv'), col_names = TRUE)

L3p2 <- cv_L3 %>%
    mutate(n_hidden = str_c(n_hidden_0, n_hidden_1, n_hidden_2, sep='_')) %>%
    filter(activate_fun == 'relu') %>%
    select(- train_mse) %>%
    group_by(n_hidden_0, n_hidden_1, n_hidden_2, dropout) %>%
    summarize(mse = mean(valid_mse)) %>% #min, mean
    ungroup() %>%
    ggplot(aes(x = n_hidden_0, y = n_hidden_1)) +
        geom_tile(aes(fill = mse)) +
        geom_text(aes(label = round(mse, 2))) +
        scale_fill_gradient(low = 'white', high = 'steelblue') +
        facet_wrap(. ~ n_hidden_2, ncol = 3) 


png('~/Desktop/L3p2.png', 4600, 2600, res=250)
print(L3p2)
dev.off()





calc_permutation_stats <- function(dpath) {
    .mse <- function(x, y) sum((x - y) ^ 2) / length(x)
    
    mse <- NULL
    for (fpath in list.files(dpath, pattern = '.tsv', full.names = TRUE)) {
        x <- read.table(fpath, header = TRUE, sep = '\t')
        mse <- c(mse, .mse(x[, 1], x[, 2]))
    }
    
    fig <- ggplot(data.frame(mse = mse), aes(x = mse)) +
                geom_histogram(binwidth = 0.2)
    fig
}


fig1 <- calc_permutation_stats(paste0(nn_result_permu, '/Type1'))
fig2 <- calc_permutation_stats(paste0(nn_result_permu, '/Type2'))


png('~/Desktop/permu1.png', 800, 600, res=250)
print(fig1)
dev.off()
png('~/Desktop/permu2.png', 800, 600, res=250)
print(fig2)
dev.off()



df <- read.table('formatted_data/kyuuri_honpo_percent.test.tsv', header = TRUE, sep = '\t')
df$incidence <- read.table('model/valid_outputs/testresult.tsv', header = TRUE, sep = '\t')$predicted

h1 <- df %>% select(prefecture, month, incidence) %>%
        ggplot(aes(x = month, y = prefecture, fill = incidence)) +
        geom_tile() + 
        scale_fill_gradientn('value', colours = brewer.pal(9, 'YlOrRd'), na.value = 'white')

png('~/Desktop/heatmap.png', 1500, 1600, res=250)
print(h1)
dev.off()









