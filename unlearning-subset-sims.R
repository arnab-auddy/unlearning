library(glmnet)


logp_       = seq(0,5)
plogstep    = 10**(0.2)
p_          = floor(500* plogstep**(logp_))
delta       = 1 # n/p
rho         = 1 # k/p
#alpha_elnet = 0.5

if(delta == 0.5){  lamb0 <- 1}
if(delta == 0.6){  lamb0 <- .5}
if(delta == 1){  lamb0 <- 0.5}

nforget <- 1000

m = 4

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
  
  Hess_full <- XtWX + lambdaS*n*diag(p)
  Hinv <- solve(Hess_full)
  
  XWX_lambda <- X %*% Hinv %*% XtW
  
  dof <- sum(diag(XWX_lambda))
  
  print(dof/p)
  
  HinvtXW <- (Hinv%*%t(X))%*% W
  
  py_hat = plogis((X%*%beta_hat)[,1])
  
  Dvec1 <- (py_hat-y)/diag(W)
  
  idxs <- matrix(NA,nforget,m)
  
  for(i in 1:nforget){ idxs[i,] = sample(1:n,m)}
  
  beta_1new_mat <- matrix(NA, p, nforget)
  beta_1new_mat2 <- matrix(NA, p, nforget)
  
  lo_pred <- matrix(NA, p, nforget)
  
  for(i in 1:nforget){
    #beta_1new_mat[,i] = beta_hat + (HinvtXW[,idxs[i,]])%*%solve(diag(m)+XWX_lambda[idxs[i,],idxs[i,]])%*%Dvec1[idxs[i,]]
    
    Umat <- t(X[idxs[i,],])%*%W[idxs[i,],idxs[i,]]
    Vmat <- -X[idxs[i,],]
    
    Lvec <- (py_hat-y)
    
    Dvec2 <- Dvec1[idxs[i,]]
    
    beta_1new_mat[,i] = beta_hat + c(Hinv%*%Umat%*%solve(diag(m)+Vmat%*%Hinv%*%Umat)%*%c(Dvec2))
    
    
    myfit_i <- glmnet(X[-idxs[i,],], y[-idxs[i,]], lambda = lambdaS, alpha = 0, intercept = FALSE, family= "binomial", standardize = FALSE)
    
    lo_pred[,i] <- coef(myfit_i)@x
    
  }
    
  err <- max(sqrt(colSums((beta_1new_mat-lo_pred)^2)))
  #err2 <- max(sqrt(colSums((beta_2new_mat-lo_pred)^2)))
  
  #eps <- 1
  
  ged_fun <- function(eps){
    
    ngauss <- rnorm(p)*err
    
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
    
    yprobs_test_lap[yprobs_test_lap>0.999] = 0.999
    yprobs_test_lap[yprobs_test_lap<0.001] = 0.001
    
    ged_gauss <- rep(0,ntest)
    ged_lap <- rep(0,ntest)
    
    for(i in 1:ntest){
      nlik_lo <- -(ytest[i]*log(yprobs_test_lo[i,])+(1-ytest[i])*log(1-yprobs_test_lo[i,]))
      nlik_gauss <- -(ytest[i]*log(yprobs_test_gauss[i,])+(1-ytest[i])*log(1-yprobs_test_gauss[i,]))
      nlik_lap <- -(ytest[i]*log(yprobs_test_lap[i,])+(1-ytest[i])*log(1-yprobs_test_lap[i,]))
      
      ged_gauss[i] <- max(abs(nlik_lo-nlik_gauss))
      ged_lap[i] <- max(abs(nlik_lo-nlik_lap))}
    
    inerr_gauss <- rep(0,nforget)
    inerr_lap <- rep(0,nforget)
    
    for(i in 1:nforget){
      nlik_lo <- -(y[idxs[i,]]*log(plogis(X[idxs[i,],]%*%lo_pred[,i]))+(1-y[idxs[i,]])*log(1-plogis(X[idxs[i,],]%*%lo_pred[,i])))
      nlik_gauss <- -(y[idxs[i,]]*log(plogis(X[idxs[i,],]%*%beta_unlearn_gauss[,i]))+(1-y[idxs[i,]])*log(1-plogis(X[idxs[i,],]%*%beta_unlearn_gauss[,i])))
      
      lap_probs <- plogis(X[idxs[i,],]%*%beta_unlearn_lap[,i])
      lap_probs[lap_probs > 0.999] = 0.999
      lap_probs[lap_probs < 0.001] = 0.001
      
      nlik_lap <- -(y[idxs[i,]]*log(lap_probs)+(1-y[idxs[i,]])*log(1-lap_probs))
      
      inerr_gauss[i] <- max(abs(sum(nlik_lo-nlik_gauss)))/m
      inerr_lap[i] <- max(abs(sum(nlik_lo-nlik_lap)))/m}
    
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