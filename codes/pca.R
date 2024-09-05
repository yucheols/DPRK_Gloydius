### load package
library(dplyr)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(caret)
library(ggfortify)

# clear working environment
rm(list = ls(all.names = T))

# turn off scientific notation
options(scipen = 999)

### import and handle data
list.files()

# full data
morpho <- readxl::read_xlsx('morpho/morphometrics_copy.xlsx') %>% as.data.frame()
head(morpho)

# check if null value is present
colSums(is.na(morpho))

# trim down to numerical
data <- morpho[, c(6,7,8,9,11,12,18,19,20,21,22)]
head(data)


### ln transform the data
data$lnTOL <- log(data$TOL)
data$lnHL <- log(data$HL)
data$lnHW <- log(data$HW)
data$lnHH <- log(data$HH)
data$lnED <- log(data$ED)
data$lnIOD <- log(data$IOD)
data$lnHLTOL <- log(data$`HL/TOL`)
data$lnHWTOL <- log(data$`HW/TOL`)
data$lnHHTOL <- log(data$`HH/TOL`)
data$lnEDTOL <- log(data$`ED/TOL`)
data$lnIODTOL <- log(data$`IOD/TOL`)

data_log <- data[12:22]
head(data_log)

### comput correlation matrix
corr_mat <- cor(data_log)
ggcorrplot(corr_mat)

### remove correlated
# find correlated vars
find.cor <- findCorrelation(corr_mat, cutoff = abs(0.8))
print(find.cor)

# remove correlated vars from data
reduced_data <- data_log[, -c(find.cor)]
print(reduced_data)

### run PCA with varimax rotation
# run PCA
pca.run <- prcomp(reduced_data, center = T, scale. = T) 
summary(pca.run)

# check loading
pca.run$x[, c(1,2)]

# are the PC centroids significantly different from one another?
wilcox.test(pca.run$x[, 1], pca.run$x[, 2])

# varimax rotation
raw.loadings <- pca.run$rotation[, 1:9] %*% diag(pca.run$sdev, 9, 9)
scores <- scale(pca.run$x[, 1:9]) %*% varimax(raw.loadings)$rotmat %>% as.data.frame()
colnames(scores) = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9')
print(scores)

# get eigenvalue
get_eigenvalue(pca.run)

# get contribution value to each PC
pca.var <- get_pca_var(pca.run)
get.contrib <- as.data.frame(pca.var$contrib)
print(get.contrib)

# get 2 highest contributing variables for each PC
sort(get.contrib$Dim.1) # HW ... SL
sort(get.contrib$Dim.2) # TAL/TOL... TAL

# get loading table
pca.run$rotation[, 1:2]

### visualize eigen
fviz_eig(pca.run, addlabels = T)

### plot viz
# bind
viz.data <- cbind(morpho$Voucher, morpho$Species, scores)
colnames(viz.data) = c('voucher', 'species', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9')

head(viz.data)

# generate polygon
hull <- viz.data %>% 
  group_by(species) %>%
  slice(chull(PC1, PC2))

### plot
# reorder species
viz.data$species = factor(viz.data$species, levels = c('G.ussuriensis', 'G.brevicauda', 'G.intermedius'))

# plot
viz.data %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = species), size = 3.5, stroke = 1.0, color = 'black', pch = 21) +
  geom_polygon(data = hull, aes(color = species, fill = species), alpha = 0.2, linewidth = NA, show.legend = F) +
  scale_fill_manual(values = c('#A9D18E', '#8FAADC', '#FFD966')) +
  xlab('PC1 - (44.0%)') + ylab('PC2 - (32.9%)') +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 1.0, color = 'darkgrey') +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 1.0, color = 'darkgrey') +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16, face = 'italic'),
        legend.position = 'top',
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 14))

# export
ggsave('PCA_results_varimax.png', width = 30, height = 22, dpi = 1000, units = 'cm')

## PCA sample label
autoplot(pca.run, data = viz.data, label = T)

