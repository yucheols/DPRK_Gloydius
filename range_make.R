#### make "rough" ranges for each sp for distribution visualization
library(ConR)

## load occs
ussuriensis <- read.csv('occs_thin/ussuriensis_thin_20km.csv')
brevicauda <- read.csv('occs_thin/brevicauda_thin_1km.csv')
intermedius <- read.csv('occs_thin/intermedius_thin_1km.csv')

head(ussuriensis[, c(2,3,4)])
head(brevicauda[, c(2,3,4)])
head(intermedius[, c(2,3,4)])

## prep data
u <- ussuriensis[, c(4,3,2)]
b <- brevicauda[, c(4,3,2)]
i <- intermedius[, c(4,3,2)]

## import mask
mask <- rgdal::readOGR('E:/Asia_shp/Asia/Asia.shp')

## make range poly
# ussuriensis
u.range <- EOO.computing(XY = u, exclude.area = T, country_map = mask, export_shp = T, write_shp = T, write_results = T, 
                         method.range = 'convex.hull', file.name = 'G.ussuriensis', show_progress = T)

# brevicauda
b.range <- EOO.computing(XY = b, exclude.area = T, country_map = mask, export_shp = T, write_shp = T, write_results = T, 
                         method.range = 'convex.hull', file.name = 'G.brevicauda', show_progress = T)

# intermedius
i.range <- EOO.computing(XY = i, exclude.area = T, country_map = mask, export_shp = T, write_shp = T, write_results = T, 
                         method.range = 'convex.hull', file.name = 'G.intermedius', show_progress = T)
