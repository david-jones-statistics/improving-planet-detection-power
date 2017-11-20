dtd2t_cov_fun <- function(t1,t2,params){
  period <- exp(params[1])
  const <- 2*pi/period
  value <- const^2*dt_cov_fun(t1,t2,params)
  return(value)
}