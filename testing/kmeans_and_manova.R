# remotes::install_github('mcomas/coda.base')
# remotes::install_github('mcomas/coda.plot')

library(coda.base)
library(coda.plot)

alimentation = read_cdp('~/Software/CoDaPack/data/alimentation.cdp')
X = alimentation[,2:10]

library(fpc)
kmr = kmeansruns(coordinates(X), krange = 2:6, runs = 25)

ggplot(data = data.table(nclust = 2:6, ch = kmr$crit[2:6]),
       aes(x = nclust, y = ch)) +
  geom_line() +
  annotate(geom = 'point', x = kmr$bestk, y = kmr$crit[kmr$bestk],
           col = 'blue', size = 4)

clr_biplot(X, group = alimentation$NorthMed)
sbp = sbp_basis(b1 = C+N~RM+WM+E+M+F+S+FV,
                b2 = C~N,
                b3 = F+FV~RM+WM+E+M+S,
                b4 = F~FV,
                b5 = WM~RM+E+M+S,
                b6 = RM~E+M+S,
                b7 = E~M+S,
                b8 = M~S, data = X)
plot_balance(sbp)
H = coordinates(X, sbp)
dplot = melt(cbind(H, cluster = paste("Cluster", kmr$cluster)), id.vars = "cluster")
ggplot(data=dplot) +
  geom_boxplot(aes(x = value, y = cluster, fill = cluster), show.legend = FALSE) +
  facet_wrap(~variable)

rm(list = ls())
pollen = fread('~/Software/CoDaPack/data/pollen.txt')
pollen = transform(pollen, group = paste0('g', group))

model = lm(coord(pinus, abies, quercus)~group, data = pollen, y = TRUE)
summary(model)
anova(model)
# Pillai's test is set by default, see other tests at ?anova.mlm
Yobs = model$y
Yexp = fitted(model)
Y0 = matrix(colMeans(Yobs), ncol = ncol(Yobs), nrow = nrow(Yobs), byrow = TRUE)

# Between sum of squares matrix
B = t(Yexp-Y0) %*% (Yexp-Y0)
B

# Within sum of squares matrix
W = t(Yobs-Yexp) %*% (Yobs-Yexp)
W

# Overall sum of squares matrix
t(Yobs-Y0) %*% (Yobs-Y0)
B + W

geometric_mean_bar_plot(pollen[,1:3], pollen$group)
