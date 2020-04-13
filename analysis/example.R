library(dplyr)
Rawdat <- read.csv("recurdeath.csv")

dat = data.frame(G = as.character(Rawdat$study), E = (Rawdat$stage), H = (Rawdat$risk_group_nonodes), y = Rawdat$fu_time, cen = Rawdat$yesdied,
                 age = Rawdat$age , bmi = Rawdat$bmi_group)
dat[Rawdat$risk_group_nonodes==3,]$H = 1
dat[Rawdat$risk_group_nonodes==1,]$H = 2
dat[Rawdat$risk_group_nonodes==4,]$H = 3
dat[Rawdat$risk_group_nonodes==2,]$H = 4

dat[which((Rawdat$Tumor_size_cat_imp==1) & (Rawdat$node_pos_imp==0)),]$E = 1
dat[which(((Rawdat$Tumor_size_cat_imp==2) & (Rawdat$node_pos_imp==0)) | ((Rawdat$Tumor_size_cat_imp==1) & (Rawdat$node_pos_imp==1))),]$E = 2
dat[which(((Rawdat$Tumor_size_cat_imp==3) & (Rawdat$node_pos_imp==0)) | ((Rawdat$Tumor_size_cat_imp==2) & (Rawdat$node_pos_imp==1))),]$E = 3
dat[which(((Rawdat$Tumor_size_cat_imp==3) & (Rawdat$node_pos_imp==1)) | (Rawdat$node_pos_imp==2)),]$E = 4

dat = filter(dat,G=='N9831'|G=='9741'|G=='49907')
dat$G = factor(as.character(dat$G))
dat = filter(dat, y>0)
dat = na.omit(dat)
dat = arrange(dat, y)

## recoding data: cen, time
time = dat$y
cen = dat$cen

## covariate
cov = as.matrix(select(dat, age , bmi))
cov = cbind(cov,model.matrix(~G, dat)[,-c(1)])
colnames(cov)[3] = '9741:49907'
colnames(cov)[4] = 'N9831:49907'

## build staging variable matrix, 01 dummy variable
n = dim(dat)[1]; p = 4; q = 4
Z = matrix(0 , nrow = n , ncol = p*q)
cellName = rep(NA, p*q)
for (i in 1:p)
  for (j in 1:q ){
    cellName[(i-1)*q + j] = paste0('E',i,'H',j)
    Z[which(dat$E==i & dat$H==j),(i-1)*q + j]=1
  }
colnames(Z) = cellName

## build patial order constraints for this cancer staging problem
## each row stands for a paritial order constraints
## For example,
##      E1H1 E1H2 E1H3 E1H4 E2H1 E2H2 E2H3 E2H4 E3H1 E3H2 E3H3 E3H4 E4H1 E4H2 E4H3 E4H4
##[1,]   -1    0    0    0    1    0    0    0    0    0    0    0    0    0    0    0
## stands for E1H1 <= E2H1

PO <- matrix(0, 2*p*q - p - q, p*q)
if(p > 1)
  for(i in 1:(q*(p - 1))) {
    PO[i, i] <- -1
    PO[i, i + q] <- 1
  }
if(q > 1)
  for(i in 1:(p*(q - 1))) {
    j <- (i - 1)%/%(q-1) + i
    PO[i + q*(p - 1), j:(j + 1)] <- c(-1, 1)
  }
colnames(PO) <- cellName

## solve opera model
result = opera.solve(time,cen,Z,PO,cov)




