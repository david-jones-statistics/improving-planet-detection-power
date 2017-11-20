set_priors <- function(){
  
  library(msm)
  
  log_model_prior <- function(x,index){sum(dnorm(x[index],0,sqrt(10),log=TRUE))}
  log_cov_prior <- function(x,bounds){sum(dtnorm(x,mean=apply(cbind(bounds[[2]],bounds[[1]]),1,mean),sd=1,lower=bounds[[1]],upper=bounds[[2]],log=TRUE))}
  log_extra_cov_prior <- function(x,bounds){sum(dtnorm(x,mean=apply(cbind(bounds[[2]],bounds[[1]]),1,mean),sd=1,lower=bounds[[1]],upper=bounds[[2]],log=TRUE))}
  log_nugget_prior <- function(x){sum(dnorm(x,0,sqrt(10),log=TRUE))}
  log_mean_prior <- function(x,bounds){-sum(log(bounds[[2]]-bounds[[1]]))}
  
  value <- list(log_model_prior=log_model_prior,log_cov_prior=log_cov_prior,log_extra_cov_prior=log_extra_cov_prior,log_nugget_prior=log_nugget_prior,log_mean_prior=log_mean_prior)
  return(value)
  
}