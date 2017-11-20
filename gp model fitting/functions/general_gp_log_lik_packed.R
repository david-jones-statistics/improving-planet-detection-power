general_gp_log_lik_packed <- function(for_like){
  
  time=for_like$time;y=for_like$y;sd=for_like$sd;cov_computer=for_like$cov_computer;model_paras=for_like$model_paras;
  cov_paras=for_like$cov_paras;extra_cov_paras=for_like$extra_cov_paras;extra_cov_fun=for_like$extra_cov_fun;
  index_extra_cov=for_like$index_extra_cov;add.nugget=for_like$add.nugget;log_nugget_vec=for_like$log_nugget_vec;
  max_values=for_like$max_values;mean_paras=for_like$mean_paras;mean_spec=for_like$mean_spec;sub_index=for_like$sub_index
  
  num_outputs <- length(model_paras)/4
  
  k_full <- cov_computer(time,time,model_paras,cov_paras,extra_cov_paras,extra_cov_fun,index_extra_cov)
  obs_error <- diag(sd^2) 
  cov_plus_nugget <- k_full + obs_error
  
  if (is.null(sub_index)){
    sub_index <- 1:(dim(k_full)[1])
  }
  
  # Nugget for modelling purposes
  if (add.nugget==TRUE){
    nugget_vec <- c()
    for (i in 1:num_outputs){
      nugget_vec <- c(nugget_vec,rep(exp(log_nugget_vec[i]),length(time)))
    }
    nugget_mat <- diag(nugget_vec)
    cov_plus_nugget <-  cov_plus_nugget + nugget_mat
  }
  
  if (prod(is.finite(cov_plus_nugget)) == 1){
    shrink <- 1
    count <- 1
    nuggest_value <- 0.01
    while (is.null(dim(shrink))){ 
      #if (count > 1){  # The stability nugget is needed so this if is omitted
      # Nugget to increase numerical stability
      cov_plus_nugget <- (1-nuggest_value)*cov_plus_nugget + nuggest_value*diag(diag(cov_plus_nugget))
      #}
      count <- count + 1
      shrink <- tryCatch(solve(cov_plus_nugget, tol=1e-30),error=function(cond){
        print("Matrix inversion error")
      })
    }
    
    log_normalizing <- -0.5*determinant(cov_plus_nugget[sub_index,sub_index], logarithm=TRUE)$modulus  
    mu_vec <- mean_model(time,y,max_values,mean_paras,mean_spec,num_outputs)[sub_index]
    quadratic <- -0.5*t(y[sub_index]-mu_vec)%*%shrink[sub_index,sub_index]%*%(y[sub_index]-mu_vec)
    key_quantites <- list(log_normalizing=as.numeric(log_normalizing),quadratic=as.numeric(quadratic),shrink=shrink,y=y,t=time) # Used to save computation for mean unpdates
    value <- as.numeric(log_normalizing + quadratic)
    #if (value > 0){
    #  value^(1/5)
    #}
  } else {
    
    print(cov_paras)
    print(cov_plus_nugget)
    key_quantites <- NULL
    value <- -743.7469 # smallest value allowed by R (converted to log scale)
    
  }
  
  return(list(value=value,key_quantites=key_quantites))
}