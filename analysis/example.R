library(dplyr)
dat <- read.csv('example.csv')
dat <- arrange(dat, y)
head(dat)

## recoding data: cen, time
time = dat$y
cen = dat$cen

## covariate
cov = as.matrix(select(dat, age , bmi))

## build staging variable matrix, 01 dummy variable
n = dim(dat)[1]; p = 4; q = 4
Z = matrix(0 , nrow = n , ncol = p*q)
cellName = rep(NA, p*q)
for (i in 1:p)
  for (j in 1:q ){
    cellName[(i-1)*q + j] = paste0('A',i,'B',j)
    Z[which(dat$A==i & dat$B==j),(i-1)*q + j]=1
  }
colnames(Z) = cellName
head(Z)

## build patial order constraints for this cancer staging problem
## input the Hasse diagram of poset Z utilizing igraph package
library(igraph)
Hasse_AB = graph(edges = c('A1B1','A1B2',  'A1B2','A1B3',  'A1B3','A1B4',
                           'A2B1','A2B2',  'A2B2','A2B3',  'A2B3','A2B4',
                           'A3B1','A3B2',  'A3B2','A3B3',  'A3B3','A3B4',
                           'A4B1','A4B2',  'A4B2','A4B3',  'A4B3','A4B4',
                           'A1B1','A2B1',  'A2B1','A3B1',  'A3B1','A4B1',
                           'A1B2','A2B2',  'A2B2','A3B2',  'A3B2','A4B2',
                           'A1B3','A2B3',  'A2B3','A3B3',  'A3B3','A4B3',
                           'A1B4','A2B4',  'A2B4','A3B4',  'A3B4','A4B4'),
                 directed = T)
plot(Hasse_AB)

## transform the Hasse diagram to coefficients inequalities
PO = as_inequalities(Hasse_AB,Z)
PO

## solve opera model
result = opera.solve(time,cen,Z,PO,cov)
result

## visualize on Hasse diagram
col = c('red', 'orange' , 'yellow' ,'green', 'blue', 'cyan', 'purple')
plot(Hasse_AB,vertex.color = col[result[[1]]])


