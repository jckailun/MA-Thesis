##### BS LR test #####
### Cohort 2 ###

g1 = 2; g2 = 3; n1 = nrow(X1); n2 = nrow(X2)
est12.bs = NULL; LR12.bs = NULL; LR02.bs = NULL; est02.bs = NULL; LR02.null.bs = NULL; LR12.null.bs = NULL
set.seed(1000)

for (b in 1:B){
  
  idx.bs = sample(1:n2, n2, replace = TRUE)
  dummy.bs = dummy
  dummy.bs[[2]] = dummy.bs[[2]][idx.bs, ]
  
  # case 1
  theta.02 = rep(0, 49)
  X2.bs = X2.c[idx.bs, ]
  out = brd(theta.02, list(X2.bs), k=16, cohort = 2, dummy.bs)
  est02.bs = rbind(est02.bs, out$zero)
  LR02.bs = c(LR02.bs, ll(out$zero, k=16, list(X2.bs), cohort = 2, dummy.bs))
  
  out.null = brd.lr(0, cohort = 2 , dummy = dummy.bs)
  LR02.null.bs = c(LR02.null.bs, out.null$loglik[length(out.null$loglik)])
  
  # case 2
  theta.12 = rep(0, 51)
  gower_dist_22 = daisy(X2.bs, metric = "gower", type = list(symm = c(11, 13:16)))
  diss22 = as.matrix(gower_dist_22)
  #c1 = Cluster_Medoids(diss11, clusters = g1, swap_phase = TRUE, verbose = TRUE, fuzzy = TRUE)
  c2 = pam(diss22, k = 3, diss = TRUE)
  
  cl2 = c2$clustering
  
  idx.d2 = list(0, 0, 0)
  
  for (i in 1:n2){
    for (j in 1:3){
      if (cl2[i] == j) {idx.d2[[j]] = c(idx.d2[[j]], i)}
    }
  }
  
  for (j in 1:3){
    idx.d2[[j]] = idx.d2[[j]][-1]}
  
  G2.bs = list(0, 0, 0)
  for (j in 1:3){
    for (i in idx.d2[[j]]){
      G2.bs[[j]] = rbind(G2.bs[[j]], X2.bs[i,])
    }
    
    G2.bs[[j]] = G2.bs[[j]][-1, ]
  }
  
  #G.bs = list(g1.d1, g2.d1)
  #idx.d1 = list(idx.g1.d1, idx.g2.d1)
  out = brd(theta.12, G2.bs, k=16, cohort = 2, dummy.bs)
  est12.bs = rbind(est12.bs, out$zero)
  LR12.bs = c(LR12.bs, ll(out$zero, G2.bs, k=16, cohort = 2, dummy.bs))
  
  out.null = brd.lr(c(0,0,0), cohort = 2 , dummy = dummy.bs)
  LR12.null.bs = c(LR12.null.bs, out.null$loglik[length(out.null$loglik)])
  
}

LR02.dist = -2*(lr.02$loglik[length(lr.02$loglik)] - LR02.bs)
plot(density(LR02.dist), main = "Bootstrap distribution of LR", sub = 'Case 1, Cohort 2')
LR02.ts = LR02
LR02.pval = mean(LR02.null.bs >= LR02.ts) ###

LR12.dist = -2*(lr.12$loglik[length(lr.12$loglik)] - LR12.bs)
plot(density(LR12.dist), main = "Bootstrap distribution of LR", sub = 'Case 2, Cohort 2')
LR12.ts = LR12
LR12.pval = mean(LR12.null.bs >= LR12.ts) ###

cbc.bs.cohort2 = -2*(LR02.bs - LR12.bs) #- ts.cbc.cohort2 ?
plot(density(cbc.bs.cohort2), main = "Bootstrap distribution of LR")
ts.cbc.cohort2 = 2*(est12$loglik[length(est12$loglik)] - est02$loglik[length(est02$loglik)])
pval.cbc.cohort2 = mean(cbc.bs.cohort2 >= ts.cbc.cohort2)

### End Cohort 2 ###

### Case 1 ###
cov02.bs = cov(est02.bs)
mu02.bs = apply(est02.bs, 2, mean)
se02.bs = sqrt(diag(cov02.bs))

pval02.bs = NULL
for (v in 1:49){
  ref = est02.bs[, v] - est02$zero[v] ###
  summand = mean(abs(ref) >= abs(est02$zero[v]))
  pval02.bs = c(pval02.bs, summand)
}

CI02.bs = NULL
for (v in 1:49){
  l = quantile(est02.bs[,v] - est02$zero[v], 0.025)
  u = quantile(est02.bs[,v] - est02$zero[v], 0.975)
  CI02.bs = rbind(CI02.bs, c(est02$zero[v]-u, est02$zero[v]-l))
}

out02.bs = cbind(est02$zero, se02.bs, CI02.bs, pval02.bs)
# write.csv(out02.bs, file='out02bs.csv')
out02.format.bs = NULL
for (i in 1:16){ out02.format.bs = rbind(out02.format.bs, rbind(out02.bs[i,], out02.bs[i+16,], out02.bs[i+32,]))}
write.csv(out02.format.bs, file = 'out02_formatted_bs.csv')
### End Case 1 ###

### Case 2 ###
cov12.bs = cov(est12.bs)
mu12.bs = apply(est12.bs, 2, mean)
se12.bs = sqrt(diag(cov12.bs))

pval12.bs = NULL
for (v in 1:51){
  ref = est12.bs[, v] - est12$zero[v] ###
  summand = mean(abs(ref) >= abs(est12$zero[v]))
  pval12.bs = c(pval12.bs, summand)
}

CI12.bs = NULL
for (v in 1:51){
  l = quantile(est12.bs[,v] - est12$zero[v], 0.025)
  u = quantile(est12.bs[,v] - est12$zero[v], 0.975) 
  CI12.bs = rbind(CI12.bs, c(est12$zero[v]-u, est12$zero[v]-l))
}

#plot(density(est12.bs[,32]))
#abline(v = est12$zero[32])

#z12.bs = est12$zero/se12.bs
out12.bs = cbind(est12$zero, se12.bs, CI12.bs, pval12.bs)
out12.format.bs = NULL
for (i in 1:16){ out12.format.bs = rbind(out12.format.bs, rbind(out12.bs[i,], out12.bs[i+16,], out12.bs[i+32,]))}
write.csv(out12.format.bs, file = 'out12_formatted_bs.csv')
### End Case 2 ###



