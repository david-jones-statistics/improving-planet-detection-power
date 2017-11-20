d2td2t_cov_fun <- function(t1,t2,params){
  period <- exp(params[1])
  r <- exp(params[2])/(1+exp(params[2]))
  t.diff <- t1-t2
  const <- 2*pi/period
  phi <- const*t.diff
  term2 <- 1-2*r*cos(phi)+r^2
  const2 <- 2*r/term2
  part1 <- const2*const^4*(cos(phi)+const2*(3*cos(phi)^2-4*sin(phi)^2)-12*const2^2*cos(phi)*sin(phi)^2+6*const2^3*sin(phi)^4)*(1+cov_fun(t1,t2,params))
  part2 <- const2*const^2*(-(2+const)*sin(phi)-2*const2*const*cos(phi)-7*const2*const*cos(phi)*sin(phi)+2*const2^2*const*sin(phi)^2+4*const2^2*const*sin(phi)^3)*dt_cov_fun(t1,t2,params)/term2
  part3 <- const2*const^2(3*cos(phi)-3*const2*sin(phi)^2)*dtdt_cov_fun(t1,t2,params)/term2
  part4 <- -const2*const*sin(phi)*dtd2t_cov_fun(t1,t2,params)
  value <- part1 + part2 + part3 + part4
  return(value)
}