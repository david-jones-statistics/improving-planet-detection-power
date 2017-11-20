cov_fun <- function(t1,t2,params){
  b <- exp(params[1])
  p <- 2
  v <- 5/2
  coeffs <- numeric(p+1)
  for (i in 0:p){
    coeffs[p+1-i] <- (gamma(p+1)/gamma(2*p+1))*(factorial(p+i)/(factorial(i)*factorial(p-i)))*(sqrt(8*v))^(p-i)
  }
  t.diff <- abs(t1-t2)
  value <- 0
  for (i in 0:p){
    value <- value + coeffs[i+1]*(t.diff/b)^i
  }
  value <- value*exp(-sqrt(5)*t.diff/b)
  return(value)
}