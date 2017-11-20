dt_cov_fun <- function(t1,t2,params){
  period <- exp(params[1])
  length.scale <- exp(params[2])
  t.diff <- t1-t2
  phi <- 2*t.diff*pi/period
  value <- cov_fun(t1,t2,params)*(pi*sin(phi)/(2*period*length.scale^2))
  return(value)
}