library(tidyverse)



nn_result_prefix <- 'nn_model/weights'

cv_L1 <- read_tsv(paste0(nn_result_prefix, '/cv_modelarch.pth_cv_stats.tsv'), col_names = TRUE)

cv_L1 %>%
    filter(activate_fun == 'relu') %>%
    select(- valid_mse) %>%
    group_by(n_hidden, dropout) %>%
    summarize(mse = mean(train_mse)) %>%
    ungroup() %>%
    pivot_wider(names_from = dropout, values_from = mse)

cv_L1 %>%
    filter(activate_fun == 'relu') %>%
    select(- train_mse) %>%
    group_by(n_hidden, dropout) %>%
    summarize(mse = mean(valid_mse)) %>%
    ungroup() %>%
    pivot_wider(names_from = dropout, values_from = mse)


cv_L1 %>%
    filter(activate_fun == 'sigmoid') %>%
    select(- valid_mse) %>%
    group_by(n_hidden, dropout) %>%
    summarize(mse = mean(train_mse)) %>%
    ungroup() %>%
    pivot_wider(names_from = dropout, values_from = mse)


cv_L1 %>%
    filter(activate_fun == 'sigmoid') %>%
    select(- train_mse) %>%
    group_by(n_hidden, dropout) %>%
    summarize(mse = mean(valid_mse)) %>%
    ungroup() %>%
    pivot_wider(names_from = dropout, values_from = mse)






cv_L2 <- read_tsv(paste0(nn_result_prefix, '/cv_modelarch_L2.pth_cv_stats.tsv'), col_names = TRUE)


p1 <- cv_L2 %>%
    filter(activate_fun == 'relu') %>%
    select(- valid_mse) %>%
    group_by(n_hidden_1, n_hidden_2, dropout) %>%
    summarize(mse = mean(train_mse)) %>%
    ungroup() %>%
    ggplot(aes(x = n_hidden_1, y = n_hidden_2)) +
        geom_tile(aes(fill = mse)) +
        geom_text(aes(label = round(mse, 2))) +
        scale_fill_gradient(low = 'white', high = 'steelblue') +
        coord_fixed() + facet_grid(. ~ dropout) 
png('~/Desktop/p1.png', 2400, 1300, res=250)
print(p1)
dev.off()

p2 <- cv_L2 %>%
    filter(activate_fun == 'relu') %>%
    select(- train_mse) %>%
    group_by(n_hidden_1, n_hidden_2, dropout) %>%
    summarize(mse = mean(valid_mse)) %>%
    ungroup() %>%
    ggplot(aes(x = n_hidden_1, y = n_hidden_2)) +
        geom_tile(aes(fill = mse)) +
        geom_text(aes(label = round(mse, 2))) +
        scale_fill_gradient(low = 'white', high = 'steelblue') +
        coord_fixed() + facet_grid(. ~ dropout) 
png('~/Desktop/p2.png', 2400, 1300, res=250)
print(p2)
dev.off()









