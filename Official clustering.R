# Clustering with PAM #
install.packages("ClusterR"); library(ClusterR)
install.packages("cluster"); library(cluster)

gower_dist_11 <- daisy(X1.c, metric = "gower", type = list(symm = c(11, 13:16)))
gower_dist_22 <- daisy(X2.c, metric = "gower", type = list(symm = c(11, 13:16)))

diss11 = as.matrix(gower_dist_11)
diss22 = as.matrix(gower_dist_22)


gower_dist_11cc <- daisy(X1.cc, metric = "gower", type = list(symm = c(12, 14:17)))
gower_dist_22cc <- daisy(X2.cc, metric = "gower", type = list(symm = c(12, 14:17)))

diss11cc = as.matrix(gower_dist_11cc)
diss22cc = as.matrix(gower_dist_22cc)


kmd = function(x,k){ list(cluster= pam(x, k, diss = TRUE)$clustering) } # with k-med
gapstat.1c = gapp(diss11, FUN = kmd, K.max = 20, B = 500, d.power = 1, spaceH0 = 'scaledPCA')
gapstat.1cc = gapp(diss11cc, FUN = kmd, K.max = 20, B = 500, d.power = 2, spaceH0 = 'scaledPCA')

print(gapstat.1c.trial, method = "Tibs2001SEmax")
G1 <- maxSE(gapstat.1c.trial$Tab[, "gap"], gapstat.1c.trial$Tab[, "SE.sim"], method="Tibs2001SEmax")
plot(gapstat.1cc, main = "Gap Statistic on D1")
abline(v=G1, lty=3, lwd=2, col="Blue")

#mtds = c('firstSEmax', 'Tibs2001SEmax', 'globalSEmax', 'firstmax', 'globalmax')

gapstat.2cc = gapp(diss22cc, FUN = kmd, K.max = 20, B = 500, spaceH0 = 'scaledPCA')
(G2 <- maxSE(gapstat.2cc$Tab[, "gap"], gapstat.2c$Tab[, "SE.sim"], method="globalmax"))


par(mfrow=c(2,1))
plot(gapstat.1c, main = "Gap Statistic on D1")
abline(v=G1, lty=3, lwd=2, col="Blue")
plot(gapstat.2cc, main = "Gap Statistic on D2")
abline(v=G2, lty=3, lwd=2, col="Blue")

# Remarks - X1.c and X2.c are the data matrices with all continuous variables standardised.