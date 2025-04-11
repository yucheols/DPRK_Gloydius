#####  revise maps

# clean working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(terra)
library(sf)
library(dplyr)
library(tidyterra)
library(ggplot2)
library(ggpubr)

#### load NK polygon
nk <- st_read('E:/Asia_shp/North Korea/PRK_adm0.shp')

####  continous predictions
# load rasters
cont <- rast(list.files(path = 'ENMs/preds/cont_nk/', pattern = '.tif$', full.names = T))
names(cont) = c('G. brevicauda', 'G. intermedius', 'G. ussuriensis')
cont <- c(cont[[3]], cont[[1]], cont[[2]])
print(cont)

# plot
cont_plot <- ggplot() +
  geom_spatraster(data = cont) +
  geom_sf(data = nk, fill = NA, color = 'black', linewidth = 1) + 
  facet_grid(~ lyr) +
  scale_x_continuous(breaks = seq(125, 131, by = 2)) + 
  scale_fill_grass_c(name = 'Suitability',
                     palette = 'inferno',
                     direction = -1,
                     na.value = NA,
                     breaks = c(0.1, 0.9),
                     labels = c('Low', 'High')) +
  theme_bw() +
  theme(strip.text = element_text(size = 18, face = 'italic'), 
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = 'bold', margin = margin(b = 10)),
        legend.text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm'))



####  binary predictions
# load rasters
bin <- rast(list.files(path = 'ENMs/preds/bin/', pattern = '.tif$', full.names = T))
names(bin) = c('G. brevicauda', 'G. intermedius', 'G. ussuriensis')
bin <- c(bin[[3]], bin[[1]], bin[[2]])
print(bin)

# set colors
# Convert raster to long-format dataframe for ggplot
bin_df <- as.data.frame(bin, xy = T, cells = T) %>%
  pivot_longer(cols = starts_with('G.'), names_to = 'species', values_to = 'presence')

# ensure presence is a factor
bin_df <- bin_df %>%
  mutate(presence = as.character(presence),
         presence = ifelse(presence == '1', paste0('1_', species), '0'),
         presence = factor(presence, levels = c('0', '1_G. ussuriensis', '1_G. brevicauda', '1_G. intermedius')))

# rearrange plotting order
bin_df$species = factor(bin_df$species, levels = c('G. ussuriensis', 'G. brevicauda', 'G. intermedius'))

# manually set color
colors <- c('0' = '#dddddd', 
            '1_G. ussuriensis' = '#A9D18E', 
            '1_G. brevicauda' = '#8FAADC', 
            '1_G. intermedius' = '#FFD966')

# plot
bin_plot <- ggplot() +
  geom_tile(data = bin_df, aes(x = x, y = y, fill = presence)) +
  geom_sf(data = nk, fill = NA, color = 'black', linewidth = 1) +
  facet_grid(~ species) +
  scale_x_continuous(breaks = seq(125, 131, by = 2)) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(0.5, 0, 0.5, 0), 'cm'))

####  combine continuous & binary
ggarrange(cont_plot, bin_plot, nrow = 2, ncol = 1, heights = c(1, 1), align = 'hv')
ggsave('ENMs/plots/gloydius_ENMs.png', dpi = 600, width = 300, height = 400, units = 'mm')
