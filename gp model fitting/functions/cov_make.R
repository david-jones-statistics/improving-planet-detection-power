cov_make <- function(t1, t2, cov_scalar,params){
  outer(t1, t2, cov_scalar, params=params)
}