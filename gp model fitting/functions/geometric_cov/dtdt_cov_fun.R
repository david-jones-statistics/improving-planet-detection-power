dtdt_cov_fun <- function(t1,t2,params){
  period <- exp(params[1])
  r <- exp(params[2])/(1+exp(params[2]))
  t.diff <- t1-t2
  const <- 2*pi/period
  phi <- const*t.diff
  term2 <- 1-2*r*cos(phi)+r^2
  part1 <- 2*r*const^2*(cos(phi)-(2*r)*sin(phi)^2/term2)*(1+cov_fun(t1,t2,params))/term2
  part2 <- 2*r*const*sin(phi)*(-dt_cov_fun(t1,t2,params))/term2
  value <- part1 + part2
  return(value)
}