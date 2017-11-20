cov_fun <- function(t1,t2,params){
  value <- 0
  t.diff <- t1-t2
  for (k in 1:(length(params)/2)){
    beta <- exp(params[(k-1)*2+1])
    lambda <- exp(params[2*k])
    value <- value + beta^2*exp(-t.diff^2/(2*lambda^2))
  }
  return(value)
}