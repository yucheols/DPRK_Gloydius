##################  Gloydius ENMs == use SDMtune
# set random seed
set.seed(555)

# encoding error prevention
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

# load libraries
library(SDMtune)
library(terra)
library(dplyr)
library(raster)
library(ggplot2)
library(extrafont)

##### part 1 ::: environmental data #####
envs <- raster::stack(list.files(path = 'envs_all', pattern = '.bil$', full.names = T))
plot(envs[[1]])

##### part 2 ::: occurrence data #####
### NES data
nes <- read.csv('occs/4th NES_2014_2018_herps_raw.csv')
head(nes)

# filter per species
unique(nes$species)

ussuriensis <- nes %>% dplyr::filter(species == 'Gloydius ussuriensis (Emelianov, 1929)') %>% dplyr::select(2,4,5)
brevicauda <- nes %>% dplyr::filter(species == 'Gloydius brevicaudus (Stejneger, 1907)') %>% dplyr::select(2,4,5)
intermedius <- nes %>% dplyr::filter(species == 'Gloydius saxatilis (Emelianov, 1937)') %>% dplyr::select(2,4,5)

# fix species names
ussuriensis$species = 'Gloydius_ussuriensis'
brevicauda$species = 'Gloydius_brevicauda'
intermedius$species = 'Gloydius_intermedius'

### GBIF / iNat data
ussuri_other <- read.csv('occs/Gloydius_ussuriensis.csv')
brevi_other <- read.csv('occs/Gloydius_brevicauda.csv')

inter_other <- read.csv('occs/G.intermedius_iNat_20231127.csv') %>% dplyr::select(36,24,23)
inter_other$scientific_name = unique(intermedius$species)

colnames(ussuri_other) = colnames(ussuriensis)
colnames(brevi_other) = colnames(ussuri_other)
colnames(inter_other) = colnames(brevi_other)

### bind occs
ussuriensis <- rbind(ussuriensis, ussuri_other)
brevicauda <- rbind(brevicauda, brevi_other)
intermedius <- rbind(intermedius, inter_other)

### thin occs ::: ussuriensis == 20km // brevicauda == 1km // intermedius == 1km
ussuriensis <- humboldt::humboldt.occ.rarefy(in.pts = ussuriensis, colxy = 2:3, rarefy.dist = 20, rarefy.units = 'km', run.silent.rar = T)
brevicauda <- humboldt::humboldt.occ.rarefy(in.pts = brevicauda, colxy = 2:3, rarefy.dist = 1, rarefy.units = 'km', run.silent.rar = T)
intermedius <- humboldt::humboldt.occ.rarefy(in.pts = intermedius, colxy = 2:3, rarefy.dist = 1, rarefy.units = 'km', run.silent.rar = T)


##### part 3 ::: background point sampling #####
## import target group points
targ.pts <- list.files(path = 'backgrounds/target_group', pattern = '.csv', full.names = T) %>%
  lapply(read.csv) %>%
  plyr::rbind.fill() %>%
  dplyr::select(4,5,6) %>%
  na.omit()

targ.pts <- thinData(coords = targ.pts[, c(2,3)], env = terra::rast(envs), x = colnames(targ.pts[2]), y = colnames(targ.pts[3]), 
                     verbose = T, progress = T)

## make density raster
targ.ras <- rasterize(targ.pts, envs, 1)
plot(targ.ras)

targ.pres <- which(values(targ.ras) == 1) 
targ.pres.locs <- coordinates(targ.ras)[targ.pres, ]

targ.dens <- MASS::kde2d(targ.pres.locs[,1], targ.pres.locs[,2], n = c(nrow(targ.ras), ncol(targ.ras)),
                         lims = c(extent(envs)[1], extent(envs)[2], extent(envs)[3], extent(envs)[4]))

targ.dens.ras <- raster(targ.dens, envs)
targ.dens.ras2 <- resample(targ.dens.ras, envs)
dens.ras <- mask(targ.dens.ras2, envs[[1]])
plot(dens.ras)

## sample bg
bg <- xyFromCell(dens.ras,
                 sample(which(!is.na(values(subset(envs, 1)))), 10000,
                        prob = values(dens.ras)[!is.na(values(subset(envs, 1)))])) %>% as.data.frame() 

points(bg)

## export bg
write.csv(bg, 'backgrounds/extracted_bg.csv')

##### part 4 ::: select environmental variables #####
ntbox::run_ntbox()

## selected == Pearson |r| > 0.7 removed == bio1 bio2 bio3 bio8 bio12 bio13 bio15 forest slope 
envs <- raster::stack(subset(envs, c('bio1', 'bio2', 'bio3', 'bio8', 'bio12', 'bio13', 'bio15', 'forest', 'slope')))

##### part 5 ::: model tuning #####
## prep swd
u.swd <- prepareSWD(species = 'G.ussuriensis', env = terra::rast(envs), p = ussuriensis[, c(2,3)], a = bg, verbose = T)
b.swd <- prepareSWD(species = 'G.brevicauda', env = terra::rast(envs), p = brevicauda[, c(2,3)], a = bg, verbose = T)
i.swd <- prepareSWD(species = 'G.intermedius', env = terra::rast(envs), p = intermedius[, c(2,3)], a = bg, verbose = T)

## NA-cleaned occs
u.cleaned <- data.frame(species = u.swd@species, long = u.swd@coords$X, lat = u.swd@coords$Y, pa = u.swd@pa)
u.cleaned <- u.cleaned %>% dplyr::filter(pa == 1)

b.cleaned <- data.frame(species = b.swd@species, long = b.swd@coords$X, lat = b.swd@coords$Y, pa = b.swd@pa)
b.cleaned <- b.cleaned %>% dplyr::filter(pa == 1)

i.cleaned <- data.frame(species = i.swd@species, long = i.swd@coords$X, lat = i.swd@coords$Y, pa = i.swd@pa)
i.cleaned <- i.cleaned %>% dplyr::filter(pa == 1)


### NA-cleaned occs
write.csv(u.cleaned, 'occs_thin/ussuriensis_thin_20km.csv')
write.csv(b.cleaned, 'occs_thin/brevicauda_thin_1km.csv')
write.csv(i.cleaned, 'occs_thin/intermedius_thin_1km.csv')


## prep randomkfold
u.fold <- randomFolds(data = u.swd, k = 10, only_presence = T, seed = 555)
b.fold <- randomFolds(data = b.swd, k = 10, only_presence = T, seed = 555)
i.fold <- randomFolds(data = i.swd, k = 10, only_presence = T, seed = 555)

## train default random 10-fold CV models
u.mod <- SDMtune::train(method = 'Maxent', data = u.swd, folds = u.fold, progress = T, iter = 5000)
b.mod <- SDMtune::train(method = 'Maxent', data = b.swd, folds = b.fold, progress = T, iter = 5000)
i.mod <- SDMtune::train(method = 'Maxent', data = i.swd, folds = i.fold, progress = T, iter = 5000)

##### tune ENMs
## make tuner function
enm.tune <- function(sp.names, models, hypers, metric, save, interactive, progress, test = NULL, env = NULL){
  out.models <- list()
  out.metrics <- list()
  
  print(paste('============================ begin ENM tuning ============================'))
  
  # init
  for (i in 1:length(models)) {
    print(paste('======== initiate ENM tuning for', sp.names[[i]], '========'))
    run <- gridSearch(model = models[[i]], hypers = hypers, metric = metric, test = test, env = env,
                      save_models = save, interactive = interactive, progress = progress)
    
    out.models[[i]] <- run@models
    out.models[[i]]$species <- sp.names[[i]]
    
    out.metrics[[i]] <- run@results
    out.metrics[[i]]$species <- sp.names[[i]]
  }
  
  # save
  out.models <- out.models
  out.metrics <- out.metrics
  
  # done
  print(paste('== ENM tuning completed for', sp.names, '!! =='))
  return(list(out.models, out.metrics))
} 

## run tuning
tune <- enm.tune(sp.names = c('G.ussuriensis', 'G.brevicauda', 'G.intermedius'),
                 models = list(u.mod, b.mod, i.mod),
                 hypers = list(fc = c('l', 'lq', 'h', 'lqh', 'lqhp', 'lqhpt'),
                               reg = seq(0.5,5, by = 0.5)),
                 metric = 'auc', save = T, interactive = F, progress = T)

print(tune)

## save models
#dir.create('model_tuning')
#saveRDS(object = tune, 'model_tuning/model_tuning.rds')

## read in saved models
#tune <- readr::read_rds('model_tuning/model_tuning.rds')

##### select optimal models
# G.ussuriensis == lqhp 1.5 == 17th model
(u.opt <- tune[[2]][[1]] %>%  dplyr::filter(test_AUC == max(test_AUC)) %>% dplyr::filter(diff_AUC == min(diff_AUC)))
u.opt.mod <- tune[[1]][[1]][[17]]

# G.brevicauda == lqhp 0.5 == 5th model
(b.opt <- tune[[2]][[2]] %>% dplyr::filter(test_AUC ==max(test_AUC)) %>% dplyr::filter(diff_AUC == min(diff_AUC)))
b.opt.mod <- tune[[1]][[2]][[5]]

# G.intermedius == lqhpt 1.0 == 12th model
(i.opt <- tune[[2]][[3]] %>% dplyr::filter(test_AUC == max(test_AUC)) %>% dplyr::filter(diff_AUC == min(diff_AUC)))
i.opt.mod <- tune[[1]][[3]][[12]]

## save individual models
#saveRDS(u.opt.mod, 'model_tuning/ussuriensis_mx.rds')
#saveRDS(b.opt.mod, 'model_tuning/brevicauda_mx.rds')
#saveRDS(i.opt.mod, 'model_tuning/intermedius_mx.rds')

## read in models
u.opt.mod <- readr::read_rds('model_tuning/ussuriensis_mx.rds')
b.opt.mod <- readr::read_rds('model_tuning/brevicauda_mx.rds')
i.opt.mod <- readr::read_rds('model_tuning/intermedius_mx.rds')

##### part 6 ::: model prediction  #####
# make function
enm.predict <- function(species, models, envs, type, clamp, progress) {
  out.preds <- raster::stack()
  
  print('============================ begin ENM predictions ============================')
  for (i in 1:length(models)) {
    print(paste('====== ENM prediction for', species[[i]], '======'))
    preds <- SDMtune::predict(object = models[[i]], data = terra::rast(envs), type = type, clamp = clamp, progress = progress) %>% 
      raster()
    names(preds) <- species[[i]]
    out.preds <- raster::stack(out.preds, preds)
  }
  out.preds <- out.preds
  print('========================== ENM predictions completed ==========================')
  return(out.preds)
}

# run
mod.pred <- enm.predict(species = c('G.ussuriensis', 'G.brevicauda', 'G.intermedius'),
                        models = list(u.opt.mod, b.opt.mod, i.opt.mod), 
                        envs = envs, type = 'cloglog', clamp = T, progress = T)
print(mod.pred)

# export
for (i in 1:nlayers(mod.pred)) {
  layer <- mod.pred[[i]]
  file <- paste0('preds/cont/', names(mod.pred)[[i]], '_cont.tif')
  writeRaster(layer, filename = file, overwrite = T)
}


#### visualize predictions
# import NK poly
nk <- rgdal::readOGR('E:/Asia_shp/North Korea/PRK_adm0.shp')
plot(nk)

#### vis each model
# G.ussuriensis
u.pred <- raster::crop(mod.pred[[1]], extent(nk))
u.pred <- raster::mask(u.pred, nk)
plot(u.pred)

# G.brevicauda
b.pred <- raster::crop(mod.pred[[2]], extent(nk))
b.pred <- raster::mask(b.pred, nk)
plot(b.pred)

# G.intermedius
i.pred <- raster::crop(mod.pred[[3]], extent(nk))
i.pred <- raster::mask(i.pred, nk)
plot(i.pred)

### export NK masked preds
nk.preds <- raster::stack(u.pred, b.pred, i.pred)

for (i in 1:nlayers(nk.preds)) {
  layer <- nk.preds[[i]]
  file <- paste0('preds/cont/', names(nk.preds)[[i]], '_nk_cont.tif')
  writeRaster(layer, filename = file, overwrite = T)
}


#####  part 7 ::: variable importance and model eval  #####
## premutation importance
u.varimp <- maxentVarImp(u.opt.mod)
b.varimp <- maxentVarImp(b.opt.mod)
i.varimp <- maxentVarImp(i.opt.mod)

## export
write.csv(u.varimp, 'variable_importance/ussuriensis.csv')
write.csv(b.varimp, 'variable_importance/brevicauda.csv')
write.csv(i.varimp, 'variable_importance/intermedius.csv')

## AUC
auc(model = u.opt.mod, test = T)
auc(model = b.opt.mod, test = T)
auc(model = i.opt.mod, test = T)

## TSS
tss(model = u.opt.mod, test = T)
tss(model = b.opt.mod, test = T)
tss(model = i.opt.mod, test = T)


#####  part 8 ::: binary maps  #####
## calculate thresholds for each species == can use the "maxentTh" function to retrieve more thresholds 
u.thresh <- thresholds(model = combineCV(u.opt.mod), type = 'cloglog')
b.thresh <- thresholds(model = combineCV(b.opt.mod), type = 'cloglog')
i.thresh <- thresholds(model = combineCV(i.opt.mod), type = 'cloglog')

print(u.thresh)
print(b.thresh)
print(i.thresh)

## make binary == use Max Train Sens Spec cloglog threshold
u.bin <- ecospat::ecospat.binary.model(Pred = u.pred, Threshold = u.thresh[3,2])
b.bin <- ecospat::ecospat.binary.model(Pred = b.pred, Threshold = b.thresh[3,2])
i.bin <- ecospat::ecospat.binary.model(Pred = i.pred, Threshold = i.thresh[3,2])

plot(u.bin)
plot(b.bin)
plot(i.bin)

## export binary
bins <- raster::stack(u.bin, b.bin, i.bin)
names(bins)

for (i in 1:nlayers(bins)) {
  r <- bins[[i]]
  layer <- paste0('preds/bin/', names(bins)[[i]], '.tif')
  writeRaster(r, filename = layer, overwrite = T)
}

## calculate ranges
bins <- raster::stack(u.bin, b.bin, i.bin)

for (i in 1:nlayers(bins)) {
  r <- bins[[i]]
  
  r[r < 1] <- NA
  r <- raster::area(r, na.rm = T, weights = F)
  r <- r[!is.na(r)]
  r_ext <- length(r)*median(r)
  print(paste0('The caluclated are for ', names(bins)[[i]], 'is ', round(r_ext, digits = 0), ' km2'))
}

#####  part 9 ::: response curves  #####
## function to pull out response curve data
respDataPull <- function(model, var, type, only_presence, marginal, species_name) {
  
  plotdata.list <- list()
  
  for (i in 1:length(var)) {
    plotdata <- plotResponse(model = model, var = var[[i]], type = type, only_presence = only_presence, marginal = marginal)
    plotdata <- ggplot2::ggplot_build(plotdata)$data
    plotdata <- plotdata[[1]]
    
    plotdata <- plotdata[, c(1:4)]
    plotdata$species <- species_name
    plotdata$var <- var[[i]]
    
    plotdata.list[[i]] <- plotdata
  }
  plotdata.df <- dplyr::bind_rows(plotdata.list) 
  return(plotdata.df)
}

## data per species
# ussuriensis
u.resp.data <- respDataPull(model = u.opt.mod, var = names(envs), type = 'cloglog', 
                            only_presence = T, marginal = F, species_name = 'G.ussuriensis')

# brevicauda
b.resp.data <- respDataPull(model = b.opt.mod, var = names(envs), type = 'cloglog',
                            only_presence = T, marginal = F, species_name = 'G.brevicauda')

# intermedius
i.resp.data <- respDataPull(model = i.opt.mod, var = names(envs), type = 'cloglog',
                            only_presence = T, marginal = F, species_name = 'G.intermedius')

## bind all, recode, and reorder variables
# bind
all.resp.data <- rbind(u.resp.data, b.resp.data, i.resp.data)
colnames(all.resp.data) = c('x', 'y', 'ymin', 'ymax', 'Species', 'var')

# recode
all.resp.data$var = dplyr::recode_factor(all.resp.data$var,
                                         'bio1' = 'Bio 1',
                                         'bio2' = 'Bio 2',
                                         'bio3' = 'Bio 3',
                                         'bio8' = 'Bio 8',
                                         'bio12' = 'Bio 12',
                                         'bio13' = 'Bio 13',
                                         'bio15' = 'Bio 15',
                                         'forest' = 'Forest cover',
                                         'slope' = 'Slope')

# reorder species & var
all.resp.data$Species <- factor(all.resp.data$Species, levels = c('G.ussuriensis', 'G.brevicauda', 'G.intermedius'))
all.resp.data$var <- factor(all.resp.data$var, 
                            levels = c('Bio 1', 'Bio 2', 'Bio 3', 'Bio 8', 'Bio 12', 'Bio 13', 'Bio 15','Forest cover', 'Slope'))

## set font
windowsFonts(a = windowsFont(family = 'Arial'))

## plot
all.resp.data %>%
  ggplot(aes(x = x, y = y, group = Species, color = Species)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ var, scales = 'free') +
  scale_color_manual(values = c('#A9D18E', '#8FAADC', '#FFD966')) +
  xlab('Value') + ylab('Suitability') +
  theme_bw() +
  theme(text = element_text(family = 'a'),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, face = 'italic'),
        legend.position = 'top')

## save
ggsave('response.png', width = 30, height = 22, dpi = 1000, units = 'cm')


