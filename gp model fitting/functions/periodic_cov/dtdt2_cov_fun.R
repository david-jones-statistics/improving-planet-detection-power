dtdt2_cov_fun <- function(t1,t2,params){

  period <- exp(params[1])
  length.scale <- exp(params[2])
  t.diff <- t1-t2
  phi <- 2*t.diff*pi/period
  term1 <- pi*sin(phi)/(2*period*length.scale^2) 
  term2 <- pi^2*cos(phi)/(period^2*length.scale^2) 
  term3 <- pi^3*sin(phi)/(period^3*length.scale^2)
  value <- dtdt_cov_fun(t1,t2,params)*term1 + 2*dt_cov_fun(t1,t2,params)*term2 + 2*cov_fun(t1,t2,params)*term3
  
  return(value)
}