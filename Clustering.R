### 3. Clustering ###
#install.packages('cluster')
library(cluster)
#install.packages('data.table'); library(data.table)

# Elbow Method #
par(mfrow=c(2,1))
K = 1:9
wss1 = sapply(K, function(k){kmeans(X1, k, nstart=25 ,iter.max = 15 )$tot.withinss})
wss1
plot(K, wss1,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters G in D1",
     ylab="Total within-clusters sum of squares", main = 'Within-cluster Sum of Squares vs No. of Clusters')

wss2 = sapply(K, function(k){kmeans(X2, k, nstart=25 ,iter.max = 15 )$tot.withinss})
wss2
plot(K, wss2,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters G in D2",
     ylab="Total within-clusters sum of squares", main = 'Within-cluster Sum of Squares vs No. of Clusters')
# End : Elbow Method #

# Silhouette Method #
silhouette.1 = function(k){
  df = X1
  #km = kmeans(df, centers = k, nstart=25) # with K-means
  kmd = Cluster_Medoids(df, k, swap_phase = TRUE, verbose = F) # with K-medoid 
  ss = silhouette(kmd$clusters, dist(df))
  #ss = silhouette(kmd, dist(df))
  mean(ss[, 3])
}

silhouette.2 = function(k){
  df = X2
  km = kmeans(df, centers = k, nstart=25) # with K-means
  #kmd = pam(df, k, cluster.only = TRUE) # with K-medoid 
  ss = silhouette(km$cluster, dist(df))
  #ss = silhouette(kmd, dist(df))
  mean(ss[, 3])
}

par(mfrow=c(2,1))
K <- 2:9
avg_sil_1 <- sapply(K, silhouette.1)
avg_sil_2 <- sapply(K, silhouette.2)

plot(K, type='b', avg_sil_1, main = 'Mean Silhouette Score on D1', xlab='Number of clusters G', ylab='Average Silhouette Scores', frame=FALSE)
plot(K, type='b', avg_sil_2, main = 'Mean Silhouette Score on D2', xlab='Number of clusters G', ylab='Average Silhouette Scores', frame=FALSE)
# End of Silouette Method #

# Gap Statistics #
km = function(x,k){ list(cluster= kmeans(x, k, nstart = 25)) } # with k-means
gapstat.1 = clusGap(as.matrix(X1), FUN = km, K.max = 9, B = 60)
gapstat.2 = clusGap(as.matrix(X2), FUN = km, K.max = 9, B = 60)

G1 <- maxSE(gapstat.1$Tab[, "gap"], gapstat.1$Tab[, "SE.sim"], method="Tibs2001SEmax")
G2 <- maxSE(gapstat.2$Tab[, "gap"], gapstat.2$Tab[, "SE.sim"], method="Tibs2001SEmax")

par(mfrow=c(2,1))
plot(gapstat.1c, main = "Gap Statistic on D1")
abline(v=G1, lty=3, lwd=2, col="Blue")
plot(gapstat.2c, main = "Gap Statistic on D2")
abline(v=G2, lty=3, lwd=2, col="Blue")
# End: Gap statistic #



# Clustering with PAM #
gower_dist_1 <- daisy(X1.c, metric = "gower", type = list(symm = c(11, 13:16)))
gower_dist_2 <- daisy(X2.c, metric = "gower", type = list(symm = c(11, 13:16)))
diss1 = as.matrix(gower_dist_1)
diss2 = as.matrix(gower_dist_2)

cm = Cluster_Medoids(X1.c, clusters = 2, swap_phase = TRUE, verbose = F)

sil = function(x, Kmax){
  s = NULL
  for (G in 2:Kmax){
    cm = Cluster_Medoids(x, clusters = G, swap_phase = TRUE, verbose = F)
    s = cbind(s, mean(cm$silhouette_matrix[,5]))
  }
  s
}

s1 = sil(diss1, 9)
s2 = sil(diss2, 9)

par(mfrow=c(2,1))
plot(2:9, type='b', s1, main = 'Mean Silhouette Score on D1', xlab='Number of clusters G', ylab='Average Silhouette Scores', frame=T)
abline(v = match(max(s1), s1)+1, lty=3, lwd=2, col="Blue")

plot(2:9, type='b', s2, main = 'Mean Silhouette Score on D2', xlab='Number of clusters G', ylab='Average Silhouette Scores', frame=T)
abline(v = match(max(s2), s2)+1, lty=3, lwd=2, col="Blue")

gower_dist_11 <- daisy(X1.c, metric = "gower", type = list(symm = c(11, 13:16)))
gower_dist_22 <- daisy(X2.c, metric = "gower", type = list(symm = c(11, 13:16)))
diss11 = as.matrix(gower_dist_11)
diss22 = as.matrix(gower_dist_22)

s11 = sil(diss11, 9)
s22 = sil(diss22, 9)

par(mfrow=c(2,1))
plot(2:9, type='b', s11, main = 'Mean Silhouette Score on D1', xlab='Number of clusters G', ylab='Average Silhouette Scores', frame=T)
abline(v = match(max(s11), s11)+1, lty=3, lwd=2, col="Blue")

plot(2:9, type='b', s22, main = 'Mean Silhouette Score on D2', xlab='Number of clusters G', ylab='Average Silhouette Scores', frame=T)
abline(v = match(max(s22), s22)+1, lty=3, lwd=2, col="Blue")

install.packages('clusterGenomics'); library(clusterGenomics)

# dissim matrix + k-med #
kmd = function(x,k){ list(cluster= pam(x, k)) } # with k-med
gapstat.1c = clusGap(diss11, FUN = kmd, K.max = 9, B = 60)
gapstat.2c = clusGap(diss22, FUN = kmd, K.max = 9, B = 60)

G1 <- maxSE(gapstat.1c$Tab[, "gap"], gapstat.1c$Tab[, "SE.sim"], method="Tibs2001SEmax")
G2 <- maxSE(gapstat.2c$Tab[, "gap"], gapstat.2c$Tab[, "SE.sim"], method="Tibs2001SEmax")

gs1 = gap(diss11, Kmax = 9, B=60)
gs2 = gap(diss22, Kmax = 20, B=60)



par(mfrow=c(2,1))
plot(gapstat.1c, main = "Gap Statistic on D1")
abline(v=G1, lty=3, lwd=2, col="Blue")
plot(gapstat.2c, main = "Gap Statistic on D2")
abline(v=G2, lty=3, lwd=2, col="Blue")

#pam1 = function(x,k){ list(cluster = Cluster_Medoids(x, k, swap_phase = TRUE, verbose = F)) } # with data themselves
#pam.diss = function(x,k){ list(cluster = pam(x, k, diss = T)) } # with dissimilarity matrix
##########


install.packages('factoextra')
fviz_dist(as.dist(diss1), gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))