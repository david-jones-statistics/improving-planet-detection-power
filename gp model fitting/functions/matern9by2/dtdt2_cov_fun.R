dtdt2_cov_fun <- function(t1,t2,params){
  
  b <- exp(params[1])
  p <- 4
  v <- 9/2
  coeffs <- numeric(p+1)
  for (i in 0:p){
    coeffs[p+1-i] <- (gamma(p+1)/gamma(2*p+1))*(factorial(p+i)/(factorial(i)*factorial(p-i)))*(sqrt(8*v)/b)^(p-i)
  }
  
  abs.t.diff <- abs(t1-t2)
  t.diff <- t2-t1 # Derivative process should be first here
  sign_diff <- sign(t.diff)
  term1 <- 0
  for (i in 1:p){
    term1 <- term1 + i*coeffs[i+1]*sign_diff*abs.t.diff^(i-1)   
    # Changes form slightly because need to actually multiply by sign before evaluating
    # becuase of the diff=0 case
  }
  term1 <- term1*exp(-3*abs.t.diff/b)
  
  term2 <- 0
  for (i in 2:p){
    term2 <- term2 + i*(i-1)*coeffs[i+1]*sign_diff*abs.t.diff^(i-2)
  }
  term2 <- term2*exp(-3*abs(t.diff)/b)
  
  term3 <- 0
  for (i in 3:p){
    term3 <- term3 + i*(i-1)*(i-2)*coeffs[i+1]*sign_diff*abs.t.diff^(i-3)
  }
  term3 <- term3*exp(-3*abs(t.diff)/b)
  
  value <- -(3/b)^2*dt_cov_fun(t1,t2,params) - (18/b^2)*term1 + 9*term2/b - term3
  
  return(value)
  
}