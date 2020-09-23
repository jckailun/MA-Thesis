# From Data Official clustering.R, we take G = (2, 3)
# Hence, PAM is applied to cluster observations accordingly

##### Case 1 #####
theta.01 = rep(0, 49)
est01 = brd(theta.01, list(X1.c), k=16, cohort = 1, dummy = dummy)
v01 = -inv(est01$hess)
z01 = est01$zero/sqrt(diag(v01))
p01 = (1-pnorm(abs(z01), 0, 1))*2
low01 = est01$zero - 1.96*sqrt(diag(v01))
up01 = est01$zero + 1.96*sqrt(diag(v01))
out01 = cbind(est01$zero, sqrt(diag(v01)), p01)
write.csv(out01, file='out01.csv')

theta.02 = rep(0, 49)
est02 = brd(theta.02, k=16, list(X2.c), cohort = 2, dummy)
v02 = -inv(est02$hess)
z02 = est02$zero/sqrt(diag(v02))
p02 = (1-pnorm(abs(z02), 0, 1))*2
low02 = est02$zero - 1.645*sqrt(diag(v02))
up02 = est02$zero + 1.645*sqrt(diag(v02))
out02 = cbind(est02$zero, sqrt(diag(v02)), p02)
write.csv(out02, file='out02.csv')
##### End of Case 1 #####

##### Case 2 #####
G1 = 2; G2 = 3; n1 = nrow(X1); n2 = nrow(X2)

gower_dist_11 <- daisy(X1.c, metric = "gower", type = list(symm = c(11, 13:16)))
gower_dist_22 <- daisy(X2.c, metric = "gower", type = list(symm = c(11, 13:16)))
diss11 = as.matrix(gower_dist_11)
diss22 = as.matrix(gower_dist_22)

#c1 = Cluster_Medoids(diss11, clusters = G1, swap_phase = TRUE, verbose = TRUE, fuzzy = TRUE)
c1 = pam(diss11, k = 2, diss= TRUE)
#c2 = Cluster_Medoids(diss22, clusters = G2, swap_phase = TRUE, verbose = TRUE, fuzzy = TRUE)
c2 = pam(diss22, k = 3, diss= TRUE)
apply(c1$fuzzy_probs, 2, mean) # mean fuzzy prob
apply(c2$fuzzy_probs, 2, mean) # mean fuzzy prob

write.csv(c1$clustering_stats, file = '2_clusters_D1.csv')
write.csv(c2$clustering_stats, file = '3_clusters_D2.csv')

#cl1 = c1$clusters
cl1 = c1$clustering
#cl2 = c2$clusters
cl2 = c2$clustering

### Grouping On D1 - Case 2 ###
idx.g1.d1 = NULL; idx.g2.d1 = NULL
for (i in 1:n1){
  idx.g1.d1 = c(idx.g1.d1, ifelse(cl1[i]==1, i, 0))
  idx.g2.d1 = c(idx.g2.d1, ifelse(cl1[i]==2, i, 0))
}
idx.g1.d1 = idx.g1.d1[idx.g1.d1!=0]
idx.g2.d1 = idx.g2.d1[idx.g2.d1!=0]

g1.d1 = NULL
for (i in idx.g1.d1){
  g1.d1 = rbind(g1.d1, X1.c[i,])
}

g2.d1 = NULL
for (i in idx.g2.d1){
  g2.d1 = rbind(g2.d1, X1.c[i,])
}

G1 = list(g1.d1, g2.d1)
idx.d1 = list(idx.g1.d1, idx.g2.d1)
### End for D1 ###

### Grouping On D2 - Case 2 ###
idx.d2 = list(0, 0, 0)

for (i in 1:n2){
  for (j in 1:3){
    if (cl2[i] == j) {idx.d2[[j]] = c(idx.d2[[j]], i)}
  }
}

for (j in 1:3){
  idx.d2[[j]] = idx.d2[[j]][-1]}

G2 = list(0, 0, 0)
for (j in 1:3){
  for (i in idx.d2[[j]]){
    G2[[j]] = rbind(G2[[j]], X2.c[i,])
  }
  
  G2[[j]] = G2[[j]][-1, ]
}
### End for D2 - Case 2###

theta.11 = rep(0, 50)
est11 = brd(theta.11, G1, k=16, cohort = 1, dummy)
v11 = -inv(est11$hess)
z11 = est11$zero/sqrt(diag(v11))
low11 = est11$zero - 1.645*sqrt(diag(v11))
up11 = est11$zero + 1.645*sqrt(diag(v11))
out11 = cbind(est11$zero, sqrt(diag(v11)), z11, low11, up11)
write.csv(out11, file='out11.csv')

theta.12 = rep(0, 51)
est12 = brd(theta.12, G2, k=16, cohort = 2, dummy)
v12 = -inv(est12$hess)
z12 = est12$zero/sqrt(diag(v12))
low12 = est12$zero - 1.645*sqrt(diag(v12))
up12 = est12$zero + 1.645*sqrt(diag(v12))
out12 = cbind(est12$zero, sqrt(diag(v12)), z12, low12, up12)
write.csv(out12, file='out12.csv')

##### End of Case 2 #####

##### Formatting the output table ####
out01.format = NULL
for (i in 1:16){ out01.format = rbind(out01.format, rbind(out01[i,], out01[i+16,], out01[i+32,]))}
write.csv(out01.format, file = 'out01_formatted.csv')

out02.format = NULL
for (i in 1:16){ out02.format = rbind(out02.format, rbind(out02[i,], out02[i+16,], out02[i+32,]))}
write.csv(out02.format, file = 'out02_formatted.csv')

out11.format = NULL
for (i in 1:16){ out11.format = rbind(out11.format, rbind(out11[i,], out11[i+16,], out11[i+32,]))}
write.csv(out11.format, file = 'out11_formatted.csv')

out12.format = NULL
for (i in 1:16){ out12.format = rbind(out12.format, rbind(out12[i,], out12[i+16,], out12[i+32,]))}
write.csv(out12.format, file = 'out12_formatted.csv')
##### End of Formatting #####

##### LR test on all 4 models #####
lr.01 = brd.lr(0, cohort = 1, dummy)
lr.02 = brd.lr(0, cohort = 2, dummy)

lr.11 = brd.lr(c(0,0), cohort = 1, dummy)
lr.12 = brd.lr(c(0,0,0), cohort = 2, dummy)

LR01 = -2*(lr.01$loglik[length(lr.01$loglik)] - est01$loglik[length(est01$loglik)])
r01 = 1 - (est01$loglik[length(est01$loglik)] / lr.01$loglik[length(lr.01$loglik)])
LR02 = -2*(lr.02$loglik[length(lr.02$loglik)] - est02$loglik[length(est02$loglik)])
r02 = 1 - (est02$loglik[length(est02$loglik)] / lr.02$loglik[length(lr.02$loglik)])

LR11 = -2*(lr.11$loglik[length(lr.11$loglik)] - est11$loglik[length(est11$loglik)])
r11 = 1 - (est11$loglik[length(est11$loglik)] / lr.11$loglik[length(lr.11$loglik)])
LR12 = -2*(lr.12$loglik[length(lr.12$loglik)] - est12$loglik[length(est12$loglik)])
r12 = 1 - (est12$loglik[length(est12$loglik)] / lr.12$loglik[length(lr.12$loglik)])

##### End of LR test on all 4 models #####