dt_cov_fun <- function(t1,t2,params){
  period <- exp(params[1])
  t.diff <- t1-t2
  const <- 2*pi/period
  phi <- const*t.diff
  value <- const*sin(phi)
  return(value)
}