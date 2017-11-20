dt2_cov_fun <- function(t1,t2,params){

  value <- -dtdt_cov_fun(t1,t2,params)
  
  return(value)
}