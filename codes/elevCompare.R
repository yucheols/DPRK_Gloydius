###### compare elevation, slope, forest cover between species
library(raster)
library(dplyr)
library(ggplot2)

### look at binaries
plot(u.bin)
plot(b.bin)
plot(i.bin)

### make copies
u.bin2 <- u.bin
b.bin2 <- b.bin
i.bin2 <- i.bin

### make 0 cells into NA
u.bin2[u.bin2 < 1] <- NA
b.bin2[b.bin2 < 1] <- NA
i.bin2[i.bin2 < 1] <- NA

plot(u.bin2)
plot(b.bin2)
plot(i.bin2)

### elev
elev <- raster('envs_all/elevation.bil')

# ussuriensis
u.elev <- raster::crop(elev, u.bin2)
u.elev <- raster::mask(u.elev, u.bin2)
plot(u.elev)

# brevicauda
b.elev <- raster::crop(elev, b.bin2)
b.elev <- raster::mask(b.elev, b.bin2)
plot(b.elev)

# intermedius
i.elev <- raster::crop(elev, i.bin2)
i.elev <- raster::mask(i.elev, i.bin2)
plot(i.elev)

### to dataframe
u.elev.df <- as.data.frame(u.elev) %>% na.omit()
b.elev.df <- as.data.frame(b.elev) %>% na.omit()
i.elev.df <- as.data.frame(i.elev) %>% na.omit()

# look at boxplot
boxplot(u.elev.df)
boxplot(b.elev.df)
boxplot(i.elev.df)

### format for ggplot2
u.elev.df$species = 'G.ussuriensis'
b.elev.df$species = 'G.brevicauda'
i.elev.df$species = 'G.intermedius'

elev.df <- rbind(u.elev.df, b.elev.df, i.elev.df)
elev.df$species = factor(elev.df$species, levels = c('G.ussuriensis', 'G.brevicauda', 'G.intermedius'))

# plot
elev.df %>%
  ggplot(aes(x = species, y = elevation)) +
  geom_boxplot(outlier.shape = NA, fill = c('#A9D18E', '#8FAADC', '#FFD966'), size = 1.2) +
  xlab('Species') + ylab('Elevation (m)') +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.text.x = element_text(face = 'italic'),
        axis.title = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

# export
ggsave('elev_distribution.png', width = 25, height = 22, dpi = 1000, units = 'cm')
  

