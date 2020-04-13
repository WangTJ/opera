#' iterative step to calculate cox model derivatives
#' @description calculate the first and second order derivatives for IRLS of Cox model.
#' @param Z a matrix of staging term
#' @param cov a matrix of covariates
#' @param time a vector of survival outcome
#' @param cen a vector of censoring status
#' @param beta a vector of staging coefficients
#' @param theta a vector of covariates coefficients
#' @param indvec an indicator of risk set
#' @return a list of derivatives
iter.step <- function(Z, cov, time, cen, beta, theta, indvec){
  n = length(time)

  if (is.null(dim(Z))) eta = Z*beta + cov %*% theta else eta <- Z %*% beta + cov %*% theta

  ## calculate eta, u, d, z (derivatives w.r.t. eta for logPL)
  atrisk =c(cumsum_rev(exp(eta)) *  cen)
  atrisk = 1/atrisk; atrisk[atrisk == Inf]=0
  u = c(cen - exp(eta)*cumsum(atrisk))
  d = c(u - cen + exp(2*eta)*cumsum(atrisk^2))
  z = ifelse(d!=0,eta - u/d,0)
  return(list(eta=eta, u=u, d=-d, z=z))
}

#' calculate negative log likelihood of Cox model
#' @param  Z a matrix of staging term
#' @param cov a matrix of covariates
#' @param cen a vector of censoring status
#' @param beta a vector of staging coefficients
#' @param theta a vector of covariates coefficients
#' @return negative log likelihood of Cox model
logPL <- function(Z, cov, cen, beta, theta)
{
  eta = Z %*% beta + cov %*% theta
  expeta = exp(eta)

  sum = sum(log(cumsum_rev(expeta))*cen)

  sum = sum - as.numeric(cen %*% eta)
  return(sum)
}


#' solve opera odel
#' @description Ordering poset element by recursive amalgamation.
#' Stratifying the staging variables with respect to the survival outcomes.
#' For computing efficiency, before inputting to functions, data should be sorted by follow-up time.
#' @param time a vector of follow-up time
#' @param cen a vector of censoring status
#' @param Z a matrix of staging variables
#' @param PO a vector of paritial order relationship staging variables
#' @param cov a matrix of covariates
#' @param maxiter maximum number of iterations in each approximation step
#' @param eps accuracy level
#' @param GIC the parameter of GIC law, the default value is 2, which is AIC law
#'
#' @return staging result of staging variables
#' @import quadprog
#' @import dplyr
#' @import Matrix
opera.solve = function(time, cen, Z, PO, cov, maxiter=10, eps=0.0001, GIC = 2){
  original_gamma = dim(Z)[2]
  categories = colnames(Z)
  result <- rep(NA,original_gamma)
  names(result) <- categories
  Z = cbind(1,Z)
  colnames(Z)[1] = 'mu'
  cellName = colnames(Z)
  covName = colnames(cov)

  S = 1 # current stage
  beta = c(1,rep(1,original_gamma))
  names(beta)=cellName
  theta = rep(0, length(covName))
  names(theta)=covName


  # prepare full constraints
  Apo0 <- cbind(0,PO)
  colnames(Apo0) <- cellName
  C.full = rbind(Apo0,cbind(0,diag(1,original_gamma)))

  while (sum(is.na(result))>0) {
    # estimate
    best_aic = 100000
    best_beta = beta
    best_theta = theta

    # find lambda_max

    lambda_now = 10000
    lambda_max = lambda_now
    while (TRUE){
      lambda = c(0,rep(lambda_now/2,original_gamma),rep(0,S-1), rep(0,length(covName)))
      n.iter=0
      diff.beta = 1000
      diff.theta = 1000
      while ((n.iter < maxiter) & (diff.beta > eps) & (diff.theta > eps)) {
        ## compute eta, u, d and z based on the current value of beta
        beta0=beta
        theta0=theta
        diff.beta1= 1000
        diff.theta1 = 1000

        iter <- iter.step(Z,cov, time, cen, beta0, theta0)

        dvec <- t(t(iter$z) %*% diag(iter$d) %*% cbind(Z,cov) - lambda )
        Dmat <- t(cbind(Z,cov)) %*% diag(iter$d) %*% cbind(Z,cov)
        Dmat <- Dmat + diag(0.0001, dim(Dmat)[1])
        QP <- solve.QP(Dmat, dvec, Amat=t(cbind(C.full, matrix(0,dim(C.full)[1],length(covName)))), bvec=rep(0,dim(C.full)[1]))
        beta <- QP$solution[1:length(cellName)]
        names(beta)=cellName
        theta <- QP$solution[(length(cellName)+1):(length(cellName)+length(covName))]
        names(theta) = covName

        diff.beta <- max(abs(beta - beta0))
        diff.theta <- max(abs(theta - theta0))
        n.iter <- n.iter + 1
      }
      beta = round(beta,3)
      if (sum(beta[2:(original_gamma+1)]!=0)>0) {
        lambda_max = lambda_now
        break
      }
      lambda_now = lambda_now / 2
    }

    lambda_min = lambda_now * 0.001
    lambda_val = c(0, exp(seq(log(lambda_min),log(lambda_max),length=40)))
    for (i in 1:length(lambda_val)){
      lambda = c(0,rep(lambda_val[i],original_gamma),rep(0,S-1), rep(0,length(covName)))
      n.iter=0
      diff.beta = 1000
      diff.theta = 1000
      while ((n.iter < maxiter) & (diff.beta > eps) & (diff.theta > eps)) {
        ## compute eta, u, d and z based on the current value of beta
        beta0=beta
        theta0=theta
        diff.beta1= 1000
        diff.theta1 = 1000

        iter <- iter.step(Z,cov, time, cen, beta0, theta0)

        dvec <- t(t(iter$z) %*% diag(iter$d) %*% cbind(Z,cov) - lambda )
        Dmat <- t(cbind(Z,cov)) %*% diag(iter$d) %*% cbind(Z,cov)
        Dmat <- Dmat + diag(0.0001, dim(Dmat)[1])
        QP <- solve.QP(Dmat, dvec, Amat=t(cbind(C.full, matrix(0,dim(C.full)[1],length(covName)))), bvec=rep(0,dim(C.full)[1]))
        beta <- QP$solution[1:length(cellName)]
        names(beta)=cellName
        theta <- QP$solution[(length(cellName)+1):(length(cellName)+length(covName))]
        names(theta) = covName

        diff.beta <- max(abs(beta - beta0))
        diff.theta <- max(abs(theta - theta0))
        n.iter <- n.iter + 1
      }
      # aic
      beta = round(beta,3)
      theta = round(theta, 3)
      aic_now = round(2*logPL(Z,cov,cen,beta,theta) - GIC*sum(beta[2:(original_gamma+1)]==0),3)
      if (aic_now < best_aic){
        best_aic = aic_now
        best_beta = beta
        best_theta = theta
        best_lambda = lambda_val[i]
      }
      if (n.iter==maxiter) print('not converge')
    }
    beta = best_beta
    print(beta)
    print(theta)
    print(paste0('lambda=',best_lambda))

    #shrink
    beta = round(beta,3)
    result[names(which(beta[2:(original_gamma+1)]==0))] = S

    newApo0 <- cbind(Apo0[, c('mu',names(which(beta[2:(original_gamma+1+S-1)]!=0)))],
                     rowSums(Apo0[,names(which(beta[2:(original_gamma+1)]==0))]))
    colnames(newApo0)[dim(newApo0)[2]] = paste0('S',S)
    newApo0 <- as.matrix(dplyr::distinct(data.frame(newApo0)))
    newApo0 <- newApo0[rowSums(newApo0 != 0) != 0, ]
    Apo0 <- newApo0

    newZ <- cbind(Z[, c('mu',names(which(beta[2:(original_gamma+1+S-1)]!=0)))],
                  rowSums(Z[,names(which(beta[2:(original_gamma+1)]==0))]))
    colnames(newZ)[dim(newZ)[2]] = paste0('S',S)
    Z <- newZ

    Spo0 <- matrix(0,S,S)
    for (i in 1:S-1){
      Spo0[i,i:(i+1)] = c(-1,1)
    }
    Spo0[S,S] = -1
    colnames(Spo0) = sapply(1:S , function(x) paste0('S',x))

    original_gamma = sum(beta[2:(original_gamma+1)]!=0)
    C.full <- rbind(Apo0 , cbind(0,as.matrix(Matrix::bdiag(diag(1,original_gamma),Spo0))))
    C.full <- as.matrix(dplyr::distinct(data.frame(C.full)))

    cellName <- colnames(Z)
    beta <- rep(1 , 1 + original_gamma + S)
    names(beta) <- cellName

    S = S + 1
  }

  result <- result[categories]
  return(list(result,theta))
}

