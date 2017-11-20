opt_function_reuse_mats <- function(opt_paras,para_type,para_opt_index,log_normalizing,shrink,for_like){
 
  num_outputs <- for_like$model_spec$num_outputs
  
  para_type_index <- which(names(for_like)==para_type)
  if (length(para_opt_index)==0){
    para_opt_index <- 1:length(for_like[[para_type_index]])
  }
  for_like[[para_type_index]][para_opt_index] <- opt_paras
  
  time <- for_like$time
  y <- for_like$y
  max_values <- for_like$max_values
  mean_paras <- for_like$mean_paras
  mean_spec <- for_like$mean_spec
  mu_vec <- mean_model(time,y,max_values,mean_paras,mean_spec,num_outputs)
  
  value <- -as.numeric((log_normalizing - 0.5*t(y-mu_vec)%*%shrink%*%(y-mu_vec)))
  
  return(value)
}