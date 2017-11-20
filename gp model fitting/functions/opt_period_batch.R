opt_period_batch <- function(para_type,basic_period_opt,para_opt_index,period_index,period_vec,plot_out_spec,bounds,for_like){
  
  para_type_index <- which(names(for_like)==para_type)
  if (length(para_opt_index)==0){
    para_opt_index <- 1:length(for_like[[para_type_index]])
  }
  paras_old <- for_like[[para_type_index]]
  period_vec <- c(for_like[[para_type_index]][period_index],period_vec)
  
  # Optimize period
  like_by_period <- numeric(length(period_vec))
  paras_store <- matrix(NA,length(period_vec),length(para_opt_index))
  if (basic_period_opt){
    for (per_num in 1:length(period_vec)){
      for_like[[para_type_index]][period_index] <- period_vec[per_num]
      like_by_period[per_num] <- like_fun2(for_like)[[1]]
    }
    for_like[[para_type_index]][period_index] <- period_vec[which(like_by_period==max(like_by_period,na.rm=TRUE))[1]]
  } else {
    for (per_num in 1:length(period_vec)){
      for_like[[para_type_index]][period_index] <- period_vec[per_num]
      initial_guess <- for_like[[para_type_index]][para_opt_index]
      if (length(para_opt_index) > 1){
        out <- optim(initial_guess,opt_function,para_type=para_type,para_opt_index=para_opt_index,for_like=for_like,
                     lower=bounds[[1]][para_opt_index],upper=bounds[[2]][para_opt_index],method="L-BFGS-B")
        paras_store[per_num,] <- out$par
        like_by_period[per_num] <- out$value
      } else {
        out <- optimize(opt_function,para_type=para_type,para_opt_index=para_opt_index,for_like=for_like,
                        lower=bounds[[1]][para_opt_index],upper=bounds[[2]][para_opt_index])
        like_by_period[per_num] <- out$objective
        paras_store[per_num] <- out$minimum
      }
    }
    for_like[[para_type_index]][para_opt_index] <- paras_store[which(like_by_period==max(like_by_period,na.rm=TRUE))[1],]
  }
  
  return(list(for_like,paras_store,cbind(period_vec,like_by_period)))
}