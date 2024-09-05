#### stress test DAPC

library(adegenet)
library(ggordiplots)

## run 1: 19GNK001,19GNK004 == G.ussuriensis /// 19GNK003 == G.brevicauda /// 19GNK002 == G.intermedius
## run 2: 19GNK001,19GNK004 == G.brevicauda /// 19GNK003 == G.ussuriensis /// 19GNK002 == G.brevicauda
## run 3: 19GNK001,19GNK004 == G.intermedius /// 19GNK003 == intermedius /// 19GNK002 == G.ussuriensis

print(morpho$Voucher)
print(data_log)

##### run2
# membership
memb2 <- c('G.brevicauda', 'G.brevicauda', rep('G.ussuriensis', 11), rep('G.intermedius', 2),
           'G.brevicauda', rep('G.intermedius', 9), 'G.ussuriensis', rep('G.brevicauda', 5))

memb2 <- as.factor(memb2)
memb2 <- factor(memb2, levels = c('G.ussuriensis', 'G.brevicauda', 'G.intermedius'))


# DAPC run2
run.dapc2 <- dapc(x = data_log, grp = memb2, n.pca = 11, n.da = 11, center = T, scale = T, pca.select = 'nbEig')

# plot
scatter(run.dapc2, col = c('#A9D18E', '#8FAADC', '#FFD966'), legend = F, pch = 20, solid = 1,
        cstar = 1, clabel = 0, cex = 3, scree.da = F, scree.pca = T, posi.pca = 'bottomright')

assignplot(run.dapc2)

### run3 
# memebership
memb3 <- c('G.intermedius', 'G.intermedius', rep('G.ussuriensis', 11), rep('G.intermedius', 2),
           'G.ussuriensis', rep('G.intermedius', 9), 'G.intermedius', rep('G.brevicauda', 5))

memb3  <- as.factor(memb3)
memb3 <- factor(memb3, levels = c('G.ussuriensis', 'G.brevicauda', 'G.intermedius'))

# DAPC run3
run.dapc3 <- dapc(x = data_log, grp = memb3, n.pca = 11, n.da = 11, center = T, scale = T, pca.select = 'nbEig')

# plot
scatter(run.dapc3, col = c('#A9D18E', '#8FAADC', '#FFD966'), legend = F, pch = 20, solid = 1,
        cstar = 1, clabel = 0, cex = 3, scree.da = F, scree.pca = T, posi.pca = 'bottomright')

assignplot(run.dapc3)
