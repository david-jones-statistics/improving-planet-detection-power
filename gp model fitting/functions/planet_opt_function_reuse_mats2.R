planet_opt_function_reuse_mats2 <- function(log_normalizing,shrink,kepler_paras,non_rv_mu_vec,for_like,rv_errs){
  
  planet_rv <- v_regress(rv_errs,for_like$time,kepler_paras)
  mu_vec <- c(planet_rv[[1]],non_rv_mu_vec)
  
  value <- -as.numeric((log_normalizing - 0.5*t(for_like$y-mu_vec)%*%shrink%*%(for_like$y-mu_vec)))
  
  return(value)
}