opt_batch_reuse_mats <- function(para_type,para_opt_index,plot_out_spec,bounds,for_like,log_normalizing,shrink){
  
  para_type_index <- which(names(for_like)==para_type)
  if (length(para_opt_index)==0){
    para_opt_index <- 1:length(for_like[[para_type_index]])
  }
  initial_guess <- for_like[[para_type_index]][para_opt_index]
  paras_old <- for_like[[para_type_index]]
  
  if (is.null(bounds)){
    out <- optim(initial_guess,opt_function_reuse_mats,para_type=para_type,para_opt_index=para_opt_index,
                 log_normalizing=log_normalizing,shrink=shrink,for_like=for_like)
    for_like[[para_type_index]][para_opt_index] <- out$par
  } else {
    if (length(para_opt_index) > 1){
      out <- optim(initial_guess,opt_function_reuse_mats,para_type=para_type,para_opt_index=para_opt_index,
                   log_normalizing=log_normalizing,shrink=shrink,for_like=for_like,
                   lower=bounds[[1]][para_opt_index],upper=bounds[[2]][para_opt_index],method="L-BFGS-B")
      for_like[[para_type_index]][para_opt_index] <- out$par
    } else {
      out <- optimize(opt_function_reuse_mats,para_type=para_type,para_opt_index=para_opt_index,
                   log_normalizing=log_normalizing,shrink=shrink,for_like=for_like,
                   lower=bounds[[1]][para_opt_index],upper=bounds[[2]][para_opt_index])
      for_like[[para_type_index]][para_opt_index] <- out$minimum
    }
  }
  
  paras_diff <- max(abs(paras_old-for_like[[para_type_index]]))
  if (plot_out_spec$show_progress){
    print(paste("Max ",para_type," change: ",round(paras_diff,5),sep=""))
  }
  
  return(for_like)
}