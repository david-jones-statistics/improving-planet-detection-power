cov_fun <- function(t1,t2,params){
  b <- exp(params[1])
  p <- matern_p
  v <- p+0.5
  const <- sqrt(2*v)
  coeffs <- numeric(p+1)
  for (i in 0:p){
    coeffs[p+1-i] <- (gamma(p+1)/gamma(2*p+1))*(factorial(p+i)/(factorial(i)*factorial(p-i)))*(sqrt(8*v)/b)^(p-i)
  }
  abs.t.diff <- abs(t1-t2)
  value <- 0
  for (i in 0:p){
    value <- value + coeffs[i+1]*abs.t.diff^i
  }
  value <- value*exp(-const*abs.t.diff/b)
  return(value)
}