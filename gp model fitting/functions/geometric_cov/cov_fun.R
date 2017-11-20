cov_fun <- function(t1,t2,params){
  period <- exp(params[1])
  r <- exp(params[2])/(1+exp(params[2]))
  t.diff <- t1-t2
  const <- 2*pi/period
  phi <- const*t.dif
  term1 <- 2*r*(cos(phi)-r)
  term2 <- 1-2*r*cos(phi)+r^2
  value <- term1/term2
  return(value)
}