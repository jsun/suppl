library(tidyverse)
library(maps)


world <- map_data('world')

png('~/Desktop/japan_map.png', 1800, 1800, res = 220)

world %>% 
    filter(region == 'Japan') %>% 
    ggplot(aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = 'white', colour = 'black', size = 1) +
    coord_fixed() + theme_void()

dev.off()




x <- read_tsv('data/data_summary/dataset_T.exif.tsv', col_names = F)


japan <- world %>% filter(region == 'Japan') 

png('~/Desktop/japan_map.png', 1800, 1800, res = 220)
ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill = 'white', colour = 'black', size = 1, data = japan) +
    coord_fixed() + theme_void() +
    geom_point(aes(x = X4, y = X3), size = 5, color = '#C7100030', data = x)
dev.off()







