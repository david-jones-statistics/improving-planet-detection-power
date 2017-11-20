dt_cov_fun <- function(t1,t2,params){
  period <- exp(params[1])
  r <- exp(params[2])/(1+exp(params[2]))
  t.diff <- t1-t2
  const <- 2*pi/period
  phi <- const*t.diff
  term2 <- 1-2*r*cos(phi)+r^2
  value <- 2*r*const*(sin(phi)/term2)*(1+cov_fun(t1,t2,params))
  return(value)
}