##### Prediction with the models #####
### Predict using the highest probability evaluated at the estimates ###
### Insignificant variables are dropped ###

## Readme:
# level of significance = 0.05
# 01: Cohort 1 - case 1
# 02: Cohort 2 - case 1
# 11: Cohort 1 - case 2
# 12: Cohort 2 - case 2
##

## Before running this script, re-run the clustering part for X1.c and X2.c. 

# Cohort 1 - Case 1
theta.01 = rep(0, 49)
est01 = brd(theta.01, list(X1.c), cohort = 1)
v01 = -inv(est01$hess)
z01 = est01$zero/sqrt(diag(v01))
order01 = match(est01$zero[abs(z01)>1.96], est01$zero)
order01.1 = order01[order01<=16]
order01.2 = order01[(order01<=32) & (order01>16)]
order01.3 = order01[order01>33]
coef01.1 = est01$zero[order01.1]
coef01.2 = est01$zero[order01.2]
coef01.3 = est01$zero[order01.3]

d1 = dummy.1[,1]; d2 = dummy.1[,2]; d3 = dummy.1[,3]; d4 = dummy.1[,4]; 
x01.1 = X1.c[, order01.1]
x01.2 = X1.c[, order01.2 - 16]
x01.3 = X1.c[, order01.3 - 32]

P01 = matrix(0, nc=4, nr=n1)
for (i in 1:n1){
  
  eta1 = as.numeric(x01.1[i,] %*% coef01.1)
  eta2 = as.numeric(x01.2[i,] %*% coef01.2)
  eta3 = as.numeric(x01.3[i,] %*% coef01.3)
  
  pi = c((1 - plogis(eta1)), (plogis(eta1)*(1 - plogis(eta2))), (plogis(eta1) * plogis(eta2) * (1 - plogis(eta3))), (plogis(eta1) * plogis(eta2) * plogis(eta3)))

  P01[i,] = pi
  
}

pred01 = rep(0, n1); correct01 = pred01; actu01 = rep(0, n1)
for (i in 1:n1){
  pred01[i] = match(max(P01[i,]), P01[i,])
  actu01[i] = match(max(dummy.1[i,]), dummy.1[i,])
  correct01[i] = ifelse(pred01[i] == match(dummy.1[i,][dummy.1[i,]==1], dummy.1[i,]), 1, 0)
}
validate01 = cbind(P01, dummy.1, pred01, correct01)

tab01 = cbind(actu01, pred01)
confusion01 = matrix(0, nc=4, nr=4)
for (j in 1:4){
  for (k in 1:4){
    njk = 0
    for (i in 1:n1){
      njk = njk + ifelse(((tab01[i,1]==j) & (tab01[i,2]==k)), 1, 0)
    }
    confusion01[j,k] = njk
  }
}
               

# Cohort 1 - Case 2 #
theta.11 = rep(0, 50)
est11 = brd(theta.11, G1, cohort = 1)
#v11 = -inv(est11$hess)
#z11 = est11$zero/sqrt(diag(v11))
order11 = match(est11$zero[pval11.bs<=0.05], est11$zero)
order11.1 = order11[order11<=16]
order11.2 = order11[(order11<=32) & (order11>16)]
order11.3 = order11[(order11>33) & (order11<=48)]
coef11.1 = est11$zero[order11.1]
coef11.2 = est11$zero[order11.2]
coef11.3 = est11$zero[order11.3]

d1 = dummy.1[,1]; d2 = dummy.1[,2]; d3 = dummy.1[,3]; d4 = dummy.1[,4]; 
x11.1 = X1.c[, order11.1]
x11.2 = X1.c[, order11.2 - 16]
x11.3 = X1.c[, order11.3 - 32]
a11 = est11$zero[49:50]

P11 = matrix(0, nc=4, nr=n1)
pred11 = rep(0, n1); correct11 = pred11; actu11 = rep(0, n1)

for (i in 1:n1){
  
  if (i %in% idx.d1[[1]]){ group = 1 }
  if (i %in% idx.d1[[2]]){ group = 2 }
  
  a = ifelse(abs(pval11.bs[ifelse(group == 1, 49, 50)]) <= 0.05, est11$zero[ifelse(group == 1, 49, 50)], 0)
  
  eta1 = as.numeric(x11.1[i,] %*% coef11.1 + a)
  eta2 = as.numeric(x11.2[i,] %*% coef11.2 + a)
  eta3 = as.numeric(x11.3[i,] %*% coef11.3 + a)
  pi = c((1 - plogis(eta1)), (plogis(eta1)*(1 - plogis(eta2))), (plogis(eta1) * plogis(eta2) * (1 - plogis(eta3))), (plogis(eta1) * plogis(eta2) * plogis(eta3)))
  P11[i,] = pi
  
  pos = match(i, idx.d1[[group]])
  ind = idx.d1[[group]][pos]
  pred11[i] = match(max(P11[ind,]), P11[ind,])
  
  actu11[i] = match(max(dummy.1[ind,]), dummy.1[ind,]) #####
  correct11[i] = ifelse(pred11[i] == actu11[i], 1, 0)
  
  #match(dummy.1[ind,][dummy.1[ind,]==1], dummy.1[ind,])
}

# validate11 = cbind(P11, dummy.1, pred11, correct11)

tab11 = cbind(actu11, pred11)
confusion11 = matrix(0, nc=4, nr=4)
for (j in 1:4){
  for (k in 1:4){
    njk = 0
    for (i in 1:n1){
      njk = njk + ifelse(((tab11[i,1]==j) & (tab11[i,2]==k)), 1, 0)
    }
    confusion11[j,k] = njk
  }
}

diag(confusion11 / apply(confusion11, 1, sum)) %*% apply(dummy.1, 2, mean)
confusion11 / apply(confusion11, 1, sum)


# Cohort 2 - Case 1
theta.02 = rep(0, 49)
est02 = brd(theta.02, list(X2.c), cohort = 2)
v02 = -inv(est02$hess)
z02 = est02$zero/sqrt(diag(v02))

order02 = match(est02$zero[abs(z02)>1.96], est02$zero)
order02.1 = order02[order02<=16]
order02.2 = order02[(order02<=32) & (order02>16)]
order02.3 = order02[(order02>33) & (order02<=48)]
coef02.1 = est02$zero[order02.1]
coef02.2 = est02$zero[order02.2]
coef02.3 = est02$zero[order02.3]

# d1 = dummy.2[,1]; d2 = dummy.2[,2]; d3 = dummy.2[,3]; d4 = dummy.2[,4]; 
x02.1 = X2.c[, order02.1]
x02.2 = X2.c[, order02.2 - 16]
x02.3 = X2.c[, order02.3 - 32]
a02 = est02$zero[length(est02$zero)]

P02 = matrix(0, nc=4, nr=n2)
for (i in 1:n2){
  
  eta1 = as.numeric(x02.1[i,] %*% coef02.1 + a02)
  eta2 = as.numeric(x02.2[i,] %*% coef02.2 + a02)
  eta3 = as.numeric(x02.3[i,] %*% coef02.3 + a02)
  
  pi = c((1 - plogis(eta1)), (plogis(eta1)*(1 - plogis(eta2))), (plogis(eta1) * plogis(eta2) * (1 - plogis(eta3))), (plogis(eta1) * plogis(eta2) * plogis(eta3)))
  
  P02[i,] = pi
  
}

pred02 = rep(0, n2); correct02 = pred02; actu02 = rep(0, n2)
for (i in 1:n2){
  pred02[i] = match(max(P02[i,]), P02[i,])
  actu02[i] = match(max(dummy.2[i,]), dummy.2[i,])
  correct02[i] = ifelse(pred02[i] == match(dummy.2[i,][dummy.2[i,]==1], dummy.2[i,]), 1, 0)
}
validate02 = cbind(P02, dummy.2, pred02, correct02)

tab02 = cbind(actu02, pred02)
confusion02 = matrix(0, nc=4, nr=4)
for (j in 1:4){
  for (k in 1:4){
    njk = 0
    for (i in 1:n2){
      njk = njk + ifelse(((tab02[i,1]==j) & (tab02[i,2]==k)), 1, 0)
    }
    confusion02[j,k] = njk
  }
}

# Cohort 2 - Case 2
theta.12 = rep(0, 51)
est12 = brd(theta.12, G2, cohort = 2)
v12 = -inv(est12$hess)
z12 = est12$zero/sqrt(diag(v12))

order12 = match(est12$zero[pval12.bs<=0.05], est12$zero)
order12.1 = order12[order12<=16]
order12.2 = order12[(order12<=32) & (order12>16)]
order12.3 = order12[(order12>33) & (order12<=48)]
coef12.1 = est12$zero[order12.1]
coef12.2 = est12$zero[order12.2]
coef12.3 = est12$zero[order12.3]

# d1 = dummy.2[,1]; d2 = dummy.2[,2]; d3 = dummy.2[,3]; d4 = dummy.2[,4]; 
x12.1 = X2.c[, order12.1]
x12.2 = X2.c[, order12.2 - 16]
x12.3 = X2.c[, order12.3 - 32]
a12 = est12$zero[49:51]

P12 = matrix(0, nc=4, nr=n2)
pred12 = rep(0, n2); correct12 = pred12; actu12 = rep(0, n2)

for (i in 1:n2){
  
  if (i %in% idx.d2[[1]]){ group = 1 }
  if (i %in% idx.d2[[2]]){ group = 2 }
  if (i %in% idx.d2[[3]]){ group = 3 }
  
  a = ifelse(abs(pval12.bs[ifelse(group == 1, 49, ifelse(group == 2, 50, 51))]) <= 0.05, est12$zero[ifelse(group == 1, 49, ifelse(group == 2, 50, 51))], 0)
  
  #pval11.bs[ifelse(group == 1, 49, 50)]) <= 0.05]
  eta1 = as.numeric(x12.1[i,] %*% coef12.1 + a)
  eta2 = as.numeric(x12.2[i,] %*% coef12.2 + a)
  eta3 = as.numeric(x12.3[i,] %*% coef12.3 + a)
  
  pi = c((1 - plogis(eta1)), (plogis(eta1)*(1 - plogis(eta2))), (plogis(eta1) * plogis(eta2) * (1 - plogis(eta3))), (plogis(eta1) * plogis(eta2) * plogis(eta3)))
  P12[i,] = pi
  
  pos = match(i, idx.d2[[group]])
  ind = idx.d2[[group]][pos]
  print(ind)
  pred12[i] = match(max(P12[ind,]), P12[ind,])
  actu12[i] = match(max(dummy.2[ind,]), dummy.2[ind,])
  # when it comes to prediction, we compare the pred to actual class to which they belong.
  
  correct12[i] = ifelse(pred12[i] == actu12[i], 1, 0)
  #validate12 = cbind(P12, dummy.2, pred12, correct12)
}

tab12 = cbind(actu12, pred12)
confusion12 = matrix(0, nc=4, nr=4)
for (j in 1:4){
  for (k in 1:4){
    njk = 0
    for (i in 1:n2){
      njk = njk + ifelse(((tab12[i,1]==j) & (tab12[i,2]==k)), 1, 0)
    }
    confusion12[j,k] = njk
  }
}

conf01 = confusion01 / apply(confusion01, 1, sum)
conf02 = confusion02 / apply(confusion02, 1, sum)
conf11 = confusion11 / apply(confusion11, 1, sum)
conf12 = confusion12 / apply(confusion12, 1, sum)

write.csv(conf02*100, file = 'conf02.csv')
write.csv(conf12*100, file = 'conf12.csv')

diag(confusion02 / apply(confusion02, 1, sum)) %*% apply(dummy.2, 2, mean)
