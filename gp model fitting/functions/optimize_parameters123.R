optimize_parameters <- function(time,y,sd,initial_value,model_spec,cov_spec,extra_cov_spec,compute_settings,plot_out_spec,fix_pred,max_values,mean_spec){
  
  if (1==2){
    initial_value <- initial_value_spec$initial_value
  }
  
  # Get model details
  model_mat <- model_spec$model_mat
  model_paras <- model_spec$model_paras_initial 
  num_outputs <- model_spec$num_outputs
  opt_index_model_paras <- which(model_paras!=0)
  add.nugget <- model_spec$add.nugget
  opt.nugget <- model_spec$opt.nugget
  opt.mean <- mean_spec$opt.mean
  mean_paras <- mean_spec$mean_paras
  opt_index_mean_paras <- mean_spec$opt_index
  log_nugget_vec <- model_spec$fixed_log_nugget
  key_planet_paras_index <- mean_spec$key_planet_paras_index
  
  # Get covariance details 
  num_cov_paras <- cov_spec$num_cov_paras
  num_extra_cov_paras <- extra_cov_spec$num_cov_paras
  paras_names <- cov_spec$paras_names
  extra_paras_names <- extra_cov_spec$paras_names
  opt_index_cov_paras <- cov_spec$opt_index_cov_paras
  cov_computer <- cov_spec$cov_computer
  gp_para_bounds <- model_spec$gp_para_bounds
  cov_para_bounds <- gp_para_bounds$cov_para_bounds
  initial_cov_para_bounds <- model_spec$initial_cov_para_bounds
  extra_cov_para_bounds <- gp_para_bounds$extra_cov_para_bounds
  opt_index_extra_cov_paras <- extra_cov_spec$opt_index_cov_paras
  index_extra_cov <- which(model_mat[,4]!=0)
  extra_cov_fun <- extra_cov_spec$extra_cov_fun
  if (model_spec$opt_stellar_period){
    stellar_period_increment <- 0.25
    stellar_period_vec <- log(seq(exp(cov_para_bounds[[1]][1]),exp(cov_para_bounds[[2]][1]),stellar_period_increment))
    extra_stellar_period_vec <- log(seq(exp(extra_cov_para_bounds[[1]][1]),exp(extra_cov_para_bounds[[2]][1]),stellar_period_increment))
  }
  
  # Initial covariance parameter values
  cov_paras <- initial_value[1:num_cov_paras]
  if (length(index_extra_cov)>0){
    extra_cov_paras <- initial_value[(num_cov_paras+1):(num_cov_paras+num_extra_cov_paras)]
  } else {
    extra_cov_paras <- rep(0,num_extra_cov_paras)
  }
  mean_paras[opt_index_mean_paras] <- initial_value[(num_cov_paras+num_extra_cov_paras+1):(num_cov_paras+num_extra_cov_paras+length(opt_index_mean_paras))]
  
  names(cov_paras) <- paras_names
  names(extra_cov_paras) <- extra_paras_names
  names(log_nugget_vec) <- rownames(model_mat)
  
  # Plot settings
  if (is.null(plot_out_spec)){
    mag_levels <- rep(1,4)
  } else {
    mag_level <- plot_out_spec$plot_settings$mag_levels
  }
  
  # Pack likelihood arguments
  for_like_full <- list(time=time,y=y,sd=sd,cov_computer=cov_computer,model_paras=model_paras,
                        cov_paras=cov_paras,extra_cov_paras=extra_cov_paras,extra_cov_fun=extra_cov_fun,
                        index_extra_cov=index_extra_cov,add.nugget=add.nugget,log_nugget_vec=log_nugget_vec,
                        max_values=max_values,mean_paras=mean_paras,mean_spec=mean_spec,model_spec=model_spec)
  
  # Step-wise optimization 
  schedule <- list()
  for (i in 1:(num_outputs-1)){
    schedule[[i]] <- 2:(i+1)
  }
  schedule[[i+1]] <- 1:num_outputs
  num_mean_paras <- length(mean_paras_bounds[[1]]) 
  end_paras_index <- (num_mean_paras-(num_outputs-1)):num_mean_paras
  
  true_location_paras_index <- location_paras_index

  
  for (opt_step in 1:length(schedule)){
    
    # # # Let exponential part dominate at first
    # if (opt_step < length(schedule) & cov_spec$cov_choice=="quasi_periodic_"){
    #   #for_like$cov_paras[2] <- runif(1,cov_para_bounds[[1]][2],cov_para_bounds[[2]][2])
    #   #for_like$cov_paras[2] <- runif(1,cov_para_bounds[[1]][2],cov_para_bounds[[2]][2])
      cov_para_bounds <- initial_cov_para_bounds
    # } else {
    # cov_para_bounds <- gp_para_bounds$cov_para_bounds
    #   extra_cov_para_bounds <- gp_para_bounds$extra_cov_para_bounds
    # }
    
    model_spec$num_outputs <- length(schedule[[opt_step]])
    
    # Data subset
    y <- c(matrix(c(for_like_full$y),ncol=num_outputs)[,schedule[[opt_step]]])
    sd <- c(matrix(c(for_like_full$sd),ncol=num_outputs)[,schedule[[opt_step]]])
    
    # Model paras subset
    model_paras_indexes <- c()
    for (i in schedule[[opt_step]]){
      model_paras_indexes <- c(model_paras_indexes,((i-1)*4+1):(i*4))
    }
    model_paras <- for_like_full$model_paras[model_paras_indexes]
    model_paras[which(model_paras != 0)] <- runif(length(which(model_paras != 0)),-1,1)
    model_mat <- t(matrix(model_paras,nrow=4))
    index_extra_cov <- which(model_mat[,4]!=0)
    opt_index_model_paras <- which(model_paras!=0)
    
    # Nugget vec subset
    log_nugget_vec <- for_like_full$log_nugget_vec[schedule[[opt_step]]]
    
    # Mean paras subset 
    num_mean_paras <- length(mean_paras_bounds[[1]]) 
    if (opt_step == num_outputs){
      mean_paras_subset <- 1:num_mean_paras
      key_planet_paras_index <- for_like_full$mean_spec$key_planet_paras_index
      other_paras_index <- setdiff(opt_index_mean_paras,key_planet_paras_index)
      mean_spec$planet_indicator <- for_like_full$mean_spec$planet_indicator
      #compute_settings$min_iter <- 2*compute_settings$min_iter
      location_paras_index <- true_location_paras_index
      #cov_para_bounds[[2]][2:3] <- c(2,2) 
      cov_para_bounds <- gp_para_bounds$cov_para_bounds
    } else {
      mean_paras_subset <- end_paras_index[schedule[[opt_step]]]
      key_planet_paras_index <- c()
      other_paras_index <- 1:model_spec$num_outputs #c() # 
      location_paras_index <- 1:model_spec$num_outputs
      mean_spec$planet_indicator <- FALSE
    }
    mean_paras <- for_like_full$mean_paras[mean_paras_subset]
    mean_paras_bounds_temp <- list()
    mean_paras_bounds_temp[[1]] <- mean_spec$mean_paras_bounds[[1]][mean_paras_subset]
    mean_paras_bounds_temp[[2]] <- mean_spec$mean_paras_bounds[[2]][mean_paras_subset]
    
    # Pack likelihood arguments
    for_like <- list(time=time,y=y,sd=sd,cov_computer=cov_computer,model_paras=model_paras,
                     cov_paras=for_like_full$cov_paras,extra_cov_paras=for_like_full$extra_cov_paras,extra_cov_fun=extra_cov_fun,
                     index_extra_cov=index_extra_cov,add.nugget=add.nugget,log_nugget_vec=log_nugget_vec,
                     max_values=max_values,mean_paras=mean_paras,mean_spec=mean_spec,model_spec=model_spec)
    
    loop_out <- main_optimize_loop(for_like,compute_settings,plot_out_spec,inflate_sd,opt_index_model_paras,opt_index_cov_paras,opt_index_extra_cov_paras,key_planet_paras_index,other_paras_index,cov_para_bounds,extra_cov_para_bounds,mean_paras_bounds_temp,stellar_period_vec,log_per_guesses,index_extra_cov,opt.nugget,opt.mean,fix_pred,location_pars_index)
    
    # Update paras
    for_like <- loop_out$for_like
    for_like_full$model_paras[model_paras_indexes] <- for_like$model_paras
    for_like_full$mean_paras[mean_paras_subset] <- for_like$mean_paras
    for_like_full$cov_paras <- for_like$cov_paras
    for_like_full$extra_cov_paras <- for_like$extra_cov_paras
    for_like_full$log_nugget_vec <- for_like$log_nugget_vec[schedule[[opt_step]]]
    
    print(paste("Opt step: ",opt_step,sep=""))
    
    print(for_like$cov_paras)
    
    print(for_like$mean_paras)
    
    store_log_like <- as.matrix(loop_out$store_log_like)
    
  }
  
  fix_pred_temp <- list(t=for_like$time,y=for_like$y,sd=for_like$sd,indicator=TRUE,temp_only=TRUE)
  plot_out <- general_gp_plot_packed(for_like,plot_out_spec,fix_pred_temp)
  rv_errs <- matrix(for_like$y,ncol=num_outputs)[,1] - plot_out[[4]][,1]
  quiet_index <- which(plot_out[[4]][,2] < -0.9)
  quiet_rv_errs <- abs(matrix(for_like$y,ncol=num_outputs)[quiet_index,1] - plot_out[[4]][quiet_index,1])
  max_err_location <- for_like$t[quiet_index[which.max(quiet_rv_errs)]]
  max_err <- max(quiet_rv_errs)
  
  model_mat <- t(matrix(for_like$model_paras,nrow=4))
  rownames(model_mat) <- rownames(model_spec$model_mat)
  colnames(model_mat) <- colnames(model_spec$model_mat)
  
  return(list(model_mat=model_mat,cov_paras=for_like$cov_paras,extra_cov_paras=for_like$extra_cov_paras,mean_paras=for_like$mean_paras,log_nugget_vec=for_like$log_nugget_vec,log_like_store=as.matrix(store_log_like),max_err_location=max_err_location,max_err=max_err,rv_errs=rv_errs,for_like=for_like))
  
}