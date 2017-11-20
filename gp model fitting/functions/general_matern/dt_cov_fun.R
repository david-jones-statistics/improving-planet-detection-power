dt_cov_fun <- function(t1,t2,params){
  
  b <- exp(params[1])
  p <- matern_p
  v <- p+0.5
  const <- sqrt(2*v)
  coeffs <- numeric(p+1)
  for (i in 0:p){
    coeffs[p+1-i] <- (gamma(p+1)/gamma(2*p+1))*(factorial(p+i)/(factorial(i)*factorial(p-i)))*(sqrt(8*v)/b)^(p-i)
  }
  
  abs.t.diff <- abs(t1-t2)
  t.diff <- t2-t1 # derivative process should be first
  sign_diff <- sign(t.diff)
  
  term1 <- 0
  if (p >= 1){
    for (i in 1:p){
      term1 <- term1 + i*coeffs[i+1]*sign_diff*abs.t.diff^(i-1)
    }
    term1 <- term1*exp(-const*abs.t.diff/b)
  }
  
  value <- -const*sign_diff*cov_fun(t1,t2,params)/b + term1
  
  return(value)
  
}