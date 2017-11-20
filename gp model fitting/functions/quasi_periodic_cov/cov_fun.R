cov_fun <- function(t1,t2,params){
  # Quasi-periodic kernel
  period <- exp(params[1])
  length.scale <- exp(params[2])
  evol.time.scale <- exp(params[3])
  t.diff <- t1-t2
  phi <- t.diff*pi/period
  value <- exp(-(sin(phi))^2/(2*length.scale^2) - t.diff^2/(2*evol.time.scale^2))
  return(value)
}