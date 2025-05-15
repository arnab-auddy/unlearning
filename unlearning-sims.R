library(glmnet)


logp_       = seq(0,5)
plogstep    = 10**(0.2)
p_          = floor(500* plogstep**(logp_))
delta       = 0.5 # n/p
rho         = 1 # k/p
#alpha_elnet = 0.5

if(delta == 0.5){  lamb0 <- 1}
if(delta == 0.6){  lamb0 <- .5}
if(delta == 1){  lamb0 <- 0.5}


# Simulation

loocv_logistic <- function(X, y, lambda0){
  
  n <- nrow(X)
  
  pred_i <-function(i){
    
    X_i <- X[-i,]
    y_i <- y[-i]
    
    myfit_i <- glmnet(X_i, y_i, lambda = lambda0, alpha = 0, intercept = FALSE, family= "binomial", standardize = FALSE)
    
    beta_hat_i <- coef(myfit_i)@x
    
    return(c(beta_hat_i))
    
  }
  return(sapply(seq(1:n), pred_i))
}


sampler = function(p){
  
  n = floor(p * delta)
  k = floor(p * rho)
  
  
  beta_star = rep(0,p)
  beta_star[1:k] = rnorm(k)#5#abs(rexp(k)*(2*rbinom(k,1,0.5)-1) * sqrt(n/(2*k)))
  X = matrix(rnorm(n*p)/sqrt(n), n, p)
  py = plogis((X%*%beta_star)[,1])
  y = rep(0, n)
  
  for(i in 1:n){y[i] = rbinom(1, size = 1, prob = py[i])}
  
  lambdaS = lamb0/n
  
  myfit <- glmnet(X, y, lambda = lambdaS, alpha = 0, intercept = FALSE, family= "binomial", standardize = FALSE)
  
  beta_hat <- coef(myfit)@x
  
  fitted_values <- predict(myfit, X, type = "response")
  
  # Compute the weight matrix W (diagonal elements are p(1 - p))
  p_hat <- c(fitted_values)
  W <- diag(p_hat * (1 - p_hat))
  
  XtW <- t(X) %*% W
  XtWX <- XtW %*% X
  XWX_lambda <- X %*% solve(XtWX + lambdaS * n* diag(ncol(X))) %*% XtW
  dof <- sum(diag(XWX_lambda))
  
  print(dof/p)
  
  Hess_full <- XtWX + lambdaS*n*diag(p)
  Hinv <- solve(Hess_full)
  
  HinvtX <- Hinv%*%t(X)
  
  xHx <- rep(0,n)
  
  diagfun <- function(i){sum(c(Hinv)*c(X[i,]%*%t(X[i,])))}
  
  #xHx <- sapply(seq(1:n), diagfun)
  
  xHx <- diag(X%*%HinvtX)
  
  py_hat = plogis((X%*%beta_hat)[,1])
  
  beta_1new_mat <- beta_hat + t(t(HinvtX)*(py_hat-y)/(-xHx*diag(W)+1))
  
  #beta_2new_mat <- matrix(0,p,n)
  #
  #for(i in 1:n){
  #  W <- diag(p_hat * (1 - p_hat))
  #  
  #  XtWi <- t(X[-i,]) %*% W[-i,-i]
  #  XtWXi <- XtWi %*% X[-i,]
  #  Hess_i <- XtWXi + lambdaS*diag(p)*n
  #  
  #  grad_i <- -t(X[-i,])%*%(y[-i]-py_hat[-i])+ lambdaS*n*beta_hat
  #  grad_i2 <- (X*(y-py_hat))[i,]
  #  
  #  beta_2new_mat[,i] = beta_hat - solve(Hess_i, grad_i)
  #  if(i%%10 ==0){print(i)}
  #}
  
  lo_pred = loocv_logistic(X, y, lambda0 = lambdaS)
  
  #lo = - mean(y * log(lo_pred) + (1 - y) * log(1 - lo_pred))
  
  #py_hat = plogis((X%*%beta_hat)[,1])
  
  #grad_mat <- t(X*(py_hat-y))
  
  err <- max(sqrt(colSums((beta_1new_mat-lo_pred)^2)))
  #err2 <- max(sqrt(colSums((beta_2new_mat-lo_pred)^2)))
  
  #eps <- 1
 
  ged_fun <- function(eps){
   
  ngauss <- rnorm(p)*err/eps
  
  nlap <- rgamma(1,shape = p)*ngauss*err/(eps*sqrt(sum(ngauss^2)))
  
  lo_noisy <- lo_pred
  beta_unlearn_lap <- beta_1new_mat + nlap
  beta_unlearn_gauss <- beta_1new_mat + ngauss
  
  ntest <- 100
  
  xtest <- matrix(rnorm(ntest*p)/sqrt(n), ntest, p)
  pytest = plogis((xtest%*%beta_star)[,1])
  
  ytest <- rep(0, ntest)
  
  for(i in 1:ntest){
    ytest[i] <- rbinom(1,1,prob = pytest[i])}
  
  yprobs_test_lo <- plogis(xtest%*%lo_noisy)
  yprobs_test_gauss <- plogis(xtest%*%beta_unlearn_gauss)
  yprobs_test_lap <- plogis(xtest%*%beta_unlearn_lap)
  
  ged_gauss <- rep(0,ntest)
  ged_lap <- rep(0,ntest)
  
  for(i in 1:ntest){
    nlik_lo <- -(ytest[i]*log(yprobs_test_lo[i,])+(1-ytest[i])*log(1-yprobs_test_lo[i,]))
    nlik_gauss <- -(ytest[i]*log(yprobs_test_gauss[i,])+(1-ytest[i])*log(1-yprobs_test_gauss[i,]))
    nlik_lap <- -(ytest[i]*log(yprobs_test_lap[i,])+(1-ytest[i])*log(1-yprobs_test_lap[i,]))
    
  ged_gauss[i] <- max(abs(nlik_lo-nlik_gauss))
  ged_lap[i] <- max(abs(nlik_lo-nlik_lap))}
  
  inerr_gauss <- rep(0,n)
  inerr_lap <- rep(0,n)
  
  for(i in 1:n){
    nlik_lo <- -(y[i]*log(plogis(sum(X[i,]*lo_pred[,i])))+(1-y[i])*log(1-plogis(sum(X[i,]*lo_pred[,i]))))
    nlik_gauss <- -(y[i]*log(plogis(sum(X[i,]*beta_unlearn_gauss[,i])))+(1-y[i])*log(1-plogis(sum(X[i,]*beta_unlearn_gauss[,i]))))
    nlik_lap <- -(y[i]*log(plogis(sum(X[i,]*beta_unlearn_lap[,i])))+(1-y[i])*log(1-plogis(sum(X[i,]*beta_unlearn_lap[,i]))))
    
    inerr_gauss[i] <- max(abs(nlik_lo-nlik_gauss))
    inerr_lap[i] <- max(abs(nlik_lo-nlik_lap))}
  
  return(c(mean(ged_lap),mean(ged_gauss),mean(inerr_lap),mean(inerr_gauss)))
}
  
  #nsamp <- 1000
  #gauss_errs <- matrix(rnorm(nsamp*p),p,nsamp)*err/eps
  #lap_errs <- matrix(0,p,nsamp)
  #for(i in 1:nsamp){
  #  ngauss <- rnorm(p)
  #  lap_errs[,i] <- rgamma(1,shape = p)*ngauss*err/(eps*sqrt(sum(ngauss^2)))
  #}
  
  #vdir <- (beta_1new_mat - lo_pred)[,sample(1:n,1)]
  
  #alpha_gauss <- rep(0, nsamp)
  #beta_gauss <- rep(0, nsamp)
  #alpha_lap <- rep(0, nsamp)
  #beta_lap <- rep(0, nsamp)
  
  #for(i in 1:nsamp){
  #  alpha_gauss[i] <- sum((gauss_errs[,i]-vdir)^2) - sum(gauss_errs[,i]^2)
  #  alpha_lap[i] <- sqrt(sum((lap_errs[,i]-vdir)^2)) - sqrt(sum(lap_errs[,i]^2))
    
  #  beta_gauss[i] <- sum((gauss_errs[,i]+vdir)^2) - sum(gauss_errs[,i]^2)
  #  beta_lap[i] <- sqrt(sum((lap_errs[,i]+vdir)^2)) - sqrt(sum(lap_errs[,i]^2))
  #}
  
  #alpha_gauss_sort <- sort(alpha_gauss)
  #alpha_lap_sort <- sort(alpha_lap)
  
  #beta_gauss_sort <- ecdf(beta_gauss)(alpha_gauss_sort)
  #beta_lap_sort <- ecdf(beta_lap)(alpha_lap_sort)
  
  return(c(p,err,
           0.25, ged_fun(0.25),
           0.5,ged_fun(0.5),
           0.75,ged_fun(0.75),
           1,ged_fun(1),1.25,ged_fun(1.25),1.5,ged_fun(1.5)))
}

sampler_all = function(){return(c(sampler(p_[1]), sampler(p_[2]), 
                                  sampler(p_[3]), sampler(p_[4]), 
                                  sampler(p_[5]), sampler(p_[6])))}