##### Gap Statistic - trial for different scaling methods #####

d1.c = read.csv("NLSY79.csv")
d2.c = read.csv("NLSY97.csv")

vars1 = c("ASVAB_gs", "ASVAB_ar" , "ASVAB_wk" , "ASVAB_pc" , "ASVAB_no" , "ASVAB_cs", "ASVAB_mk", "ASVAB_mc", "ASVAB_ei", 
          "med", paste('acgrd', 1:13, sep = '_'), 'male', 'income', 'nuclear', 'urban', 'race')
vars1 = match(vars1, names(d1.c))
d1.c = d1.c[, vars1]
black1 = ifelse(d1.c$race==1, 1, 0); hisp1 = ifelse(d1.c$race==2, 1, 0); others1 = ifelse(d1.c$race == 3, 1, 0)
d1.c = as.data.frame(cbind(d1.c, black1, hisp1))
d1.c = d1.c[,-28]

rowstodelete1 = NULL
for (i in 1:nrow(d1.c)){
  if (sum(is.na(d1.c[i,])) != 0){
    rowstodelete1 = c(rowstodelete1, i)}
}

d1.c = d1.c[-rowstodelete1,]

vars2 = c("ASVAB_GS", "ASVAB_AR" , "ASVAB_WK" , "ASVAB_PC" , "ASVAB_NO" , "ASVAB_CS", "ASVAB_MK", "ASVAB_MC", 'ASVAB_EI', 
          "med", paste('acgrd', 1:13, sep = '_'), 'sex', 'income', 'nuclear', 'urban', 'race')
vars2 = match(vars2, names(d2.c))
d2.c = d2.c[, vars2]
d2.c$race = ifelse(d2.c$race==4, 3, d2.c$race); d2.c$sex = ifelse(d2.c$sex==1, 1, 0)
black2 = ifelse(d2.c$race==1, 1, 0); hisp2 = ifelse(d2.c$race==2, 1, 0); others2 = ifelse(d2.c$race == 3, 1, 0)
d2.c = as.data.frame(cbind(d2.c, black2, hisp2))
d2.c = d2.c[,-28]

rowstodelete2 = NULL
for (i in 1:nrow(d2.c)){
  if (sum(is.na(d2.c[i,])) != 0){
    rowstodelete2 = c(rowstodelete2, i)}
}

d2.c = d2.c[-rowstodelete2,]

# Next, let's compute the rates of having 0s' in ACGRDs
zeros1 = NULL
for (j in 11:23){ zeros1 = c(zeros1, sum((d1.c[,j]==0))) }
rate1 = zeros1/nrow(d1.c); margin1 = diff(zeros1)/nrow(d1.c)
# perhaps we shall choose asgrd_9 as the optimal grade transition variable.

zeros2 = NULL
for (j in 11:23){ zeros2 = c(zeros2, sum((d2.c[,j]==0))) }
rate2 = zeros2/nrow(d2.c); margin2 = diff(zeros2)/nrow(d2.c)
# perhaps we shall choose asgrd_9 as the optimal grade transition variable.

rowstodelete1 = NULL
for (i in 1:nrow(d1.c)){
  if (d1.c$acgrd_9[i] == 0){
    rowstodelete1 = c(rowstodelete1, i)}
}

rowstodelete2 = NULL
for (i in 1:nrow(d2.c)){
  if (d2.c$acgrd_10[i] == 0){
    rowstodelete2 = c(rowstodelete2, i)}
}

d1.c = d1.c[-rowstodelete1,]
d2.c = d2.c[-rowstodelete2,]

rowstodelete = NULL
for (i in 1:nrow(d1.c)){
  if (d1.c$med[i] == 0){
    rowstodelete = c(rowstodelete, i)}
}
d1.c = d1.c[-rowstodelete, ]

n1.c = nrow(d1.c)
n2.c = nrow(d2.c)

X1.c = as.matrix(d1.c); X1.c = X1.c[,-(11:23)]
X2.c = as.matrix(d2.c); X2.c = X2.c[,-(11:23)]
