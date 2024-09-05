####  try DAPC with morpho data
library(adegenet)
library(ggordiplots)

# first, look at raw data
print(morpho)

# now look at the log-transformed data
print(data_log)

# group membership
memb <- as.factor(morpho$Species)
memb <- factor(memb, levels = c('G.ussuriensis', 'G.brevicauda', 'G.intermedius'))
print(memb)

# run DAPC
run.dapc <- dapc(x = data_log, grp = memb, n.pca = 11, n.da = 11, center = T, scale = T, pca.select = 'nbEig', perc.pca = 99)

# check results
print(run.dapc)
summary(run.dapc)

# plot
scatter(run.dapc, col = c('#A9D18E', '#8FAADC', '#FFD966'), legend = F, pch = 20, solid = 1,
        cstar = 1, clabel = 0, cex = 3, scree.da = F, scree.pca = T, posi.pca = 'bottomright')

assignplot(run.dapc)


##############################   plot customization
#### try ggplot2 style
print(run.dapc$ind.coord)

### gg_ordiplot
ellips.plot <- gg_ordiplot(run.dapc$ind.coord, groups = run.dapc$grp, spiders = T, ellipse = T, kind = 'sd', conf = 0.62) 

### further customization
# ellipse data 
df_ellips <- ellips.plot$df_ellipse

ussuri_ellips <- subset(df_ellips, Group == 'G.ussuriensis', droplevels = T) 
brevi_ellips <- subset(df_ellips, Group == 'G.brevicauda', droplevels = T)
inter_ellips <- subset(df_ellips, Group == 'G.intermedius', droplevels = T)

# spider data
df_spider <- ellips.plot$df_spiders

# plot
ellips.plot$plot +
  geom_path(data = ussuri_ellips, aes(x = x, y = y), color = '#A9D18E', linewidth = 1.3) +
  geom_path(data = brevi_ellips, aes(x = x, y = y), color = '#8FAADC', linewidth = 1.3) +
  geom_path(data = inter_ellips, aes(x = x, y = y), color = '#FFD966', linewidth = 1.3) +
  geom_segment(aes(x = df_spider$cntr.x, y = df_spider$cntr.y, xend = df_spider$x, yend = df_spider$y, 
                   color = df_spider$Group), linewidth = 1.3) +
  geom_point(aes(x = ellips.plot$df_ord$x, y = ellips.plot$df_ord$y, color = ellips.plot$df_ord$Group, size = 2)) +
  scale_color_manual(values = c('#A9D18E', '#8FAADC', '#FFD966')) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 1.0, color = 'darkgrey') +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 1.0, color = 'darkgrey') +
  theme_bw() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 14)) 

# export
ggsave('dapc_custom_plot.png', width = 30, height = 22, dpi = 1000, units = 'cm')
