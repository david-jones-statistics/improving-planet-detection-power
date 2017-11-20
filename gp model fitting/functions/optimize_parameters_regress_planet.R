optimize_parameters <- function(time,y,sd,initial_value,model_spec,cov_spec,extra_cov_spec,compute_settings,plot_out_spec,fix_pred,max_values,mean_spec){
  
  if (1==2){
    initial_value <- initial_value_spec$initial_value
  }
  
  # Get model details
  model_mat <- model_spec$model_mat
  model_paras <- model_spec$model_paras_initial 
  num_outputs <- dim(model_mat)[1]
  opt_index_model_paras <- which(model_paras!=0)
  add.nugget <- model_spec$add.nugget
  opt.nugget <- model_spec$opt.nugget
  opt.mean <- mean_spec$opt.mean
  mean_paras <- mean_spec$model_paras_initial
  opt_index_mean_paras <- mean_spec$opt_index
  log_nugget_vec <- model_spec$fixed_log_nugget
  key_planet_paras_index <- mean_spec$key_planet_paras_index
  if (length(key_planet_paras_index) > 0){
    other_paras_index <- setdiff(opt_index_mean_paras,key_planet_paras_index)
    planet_period_index <- mean_spec$planet_period_index
    if (mean_spec$planet_period_search){
      log_per_guesses <- mean_spec$log_per_guesses
    }
  } else {
    other_paras_index <- opt_index_mean_paras
  }
  
  # Get covariance details 
  num_cov_paras <- cov_spec$num_cov_paras
  num_extra_cov_paras <- extra_cov_spec$num_cov_paras
  paras_names <- cov_spec$paras_names
  extra_paras_names <- extra_cov_spec$paras_names
  opt_index_cov_paras <- cov_spec$opt_index_cov_paras
  cov_computer <- cov_spec$cov_computer
  gp_para_bounds <- model_spec$gp_para_bounds
  cov_para_bounds <- gp_para_bounds$cov_para_bounds
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
  for_like <- list(time=time,y=y,sd=sd,cov_computer=cov_computer,model_paras=model_paras,
                   cov_paras=cov_paras,extra_cov_paras=extra_cov_paras,extra_cov_fun=extra_cov_fun,
                   index_extra_cov=index_extra_cov,add.nugget=add.nugget,log_nugget_vec=log_nugget_vec,
                   max_values=max_values,mean_paras=mean_paras,mean_spec=mean_spec,model_spec=model_spec)
  
  # Initialize
  count <- 1
  log_like_diff <- 100
  store_log_like <- list()
  store_log_like[[1]] <- like_fun2(for_like)[[1]]
  
  # For annealing-type approach
  true_sd <- sd
  
  # Main loop
  while ((count < compute_settings$min_iter | log_like_diff > compute_settings$log_like_tol) &
         count < compute_settings$max_iter){
    
    count <- count + 1 
    
    # Let exponential part dominate at first
    # if (count < compute_settings$min_iter & cov_spec$cov_choice=="quasi_periodic_"){
    #   cov_para_bounds[[1]][2:3] <- c(2,-1.5)
    #   extra_cov_para_bounds[[1]][2:3] <- c(2,-1.5)
    # } else {
    #   cov_para_bounds <- gp_para_bounds$cov_para_bounds
    #   extra_cov_para_bounds <- gp_para_bounds$extra_cov_para_bounds
    # }
    
    # Annealing-type approach: inflate standard deviations
    sd <- true_sd*compute_settings$inflate_sd[count]
    
    if (count==compute_settings$max_iter){
      print("Warning: max. iters reached - check convergence")
    }
    
    # Model coefficients
    for_like <- opt_batch("model_paras",opt_index_model_paras,plot_out_spec,NULL,for_like)
    
    # GP covaraince paramters 
    basic_period_opt <- TRUE
    period_out <- opt_period_batch("cov_paras",basic_period_opt,opt_index_cov_paras,1,stellar_period_vec,plot_out_spec,cov_para_bounds,for_like)
    for_like <- period_out[[1]]
    if (basic_period_opt){
      for_like <- opt_batch("cov_paras",opt_index_cov_paras,plot_out_spec,cov_para_bounds,for_like)
    }
    
    # Extra GP covaraince paramters 
    if (length(index_extra_cov)>0){
      basic_period_opt <- TRUE
      period_out <- opt_period_batch("extra_cov_paras",basic_period_opt,opt_index_extra_cov_paras,1,extra_stellar_period_vec,plot_out_spec,extra_cov_para_bounds,for_like)
      for_like <- period_out[[1]]
      if (basic_period_opt){
        for_like <- opt_batch("extra_cov_paras",opt_index_extra_cov_paras,plot_out_spec,extra_cov_para_bounds,for_like)
      }
    }
    
    # Nuggets (for modeling)
    if (opt.nugget){
      log_nugget_bounds <- list()
      log_nugget_bounds[[1]] <- rep(-compute_settings$opt_bound,num_outputs)
      log_nugget_bounds[[2]] <- rep(compute_settings$opt_bound,num_outputs)
      for_like <- opt_batch("log_nugget_vec",c(),plot_out_spec,log_nugget_bounds,for_like)
    }
    
    # Mean parameters (including planet)
    if (opt.mean){

      key_quantities <- like_fun2(for_like)[[2]]
      shrink <- key_quantities$shrink
      log_normalizing <- key_quantities$log_normalizing
      if (length(key_planet_paras_index) > 0){
  
        temp_log_scale_index <- c()
        temp_logit_scale_index <- c()
        temp_mean_paras_bounds <- list() 
        temp_mean_paras_bounds[[1]] <- rep(-2,num_outputs) # bounds have to be on optimization scale
        temp_mean_paras_bounds[[2]] <- rep(2,num_outputs) # bounds have to be on optimization scale
        temp_mean_paras <-  rep(0,num_outputs)
        names(temp_mean_paras) <- paste("mean_proxy",1:num_outputs,sep="")
        temp_mean_spec <- list(opt.mean=TRUE,opt_index=1:num_outputs,planet_indicator=FALSE,planet_period_index=NULL,
                          key_planet_paras_index=c(),mean_paras_bounds=temp_mean_paras_bounds,
                          planet_para_names=NULL,log_scale_index=temp_log_scale_index,logit_scale_index=temp_logit_scale_index,
                          mean_paras=temp_mean_paras,planet_period_search=FALSE)
        temp_for_like <- for_like
        temp_for_like$mean_spec <- temp_mean_spec 
        # Corrected 1st May 2017 (previously only had indexes 5 and 8, not 7):
        #temp_for_like$mean_paras <- for_like$mean_paras[c(5,5+num_outputs)]
        temp_for_like$mean_paras <- for_like$mean_paras[c(5,(length(for_like$mean_paras)-num_outputs+2):length(for_like$mean_paras))]
        
        fix_pred_temp <- list(t=for_like$time,y=for_like$y,sd=for_like$sd,indicator=TRUE,temp_only=TRUE)
        plot_out <- general_gp_plot_packed(temp_for_like,plot_out_spec,fix_pred_temp)
        rv_errs <- matrix(temp_for_like$y,ncol=num_outputs)[,1] - plot_out[[4]][,1]
        planet_out <- planet_optimizer_reuse_mats2("mean_paras",key_planet_paras_index,3,2,log_normalizing,shrink,mean_paras_bounds,for_like,rv_errs)
        for_like$mean_paras[key_planet_paras_index] <- planet_out$out$par
      }
      
      if (length(other_paras_index)>=1){
        for_like <-  opt_batch_reuse_mats("mean_paras",other_paras_index,plot_out_spec,mean_paras_bounds,for_like,log_normalizing,shrink)
      }
    }
    
    # Store likelihood
    store_log_like[[count]] <- like_fun2(for_like)[[1]]
    log_like_diff <- abs(store_log_like[[count]]-store_log_like[[count-1]])
    if (plot_out_spec$show_log_like_progress){
      log_like_diff_now <- store_log_like[[count]]-store_log_like[[count-1]]
      print(paste("Log like change: ",round(log_like_diff_now,2-log10(compute_settings$log_like_tol)),sep=""))
    }
    
    if (plot_out_spec$show_progress){
      print("Cov paras: ")
      para_type_index <- which(names(for_like)=="cov_paras")
      print(round(for_like[[para_type_index]],2))
    }
    
    if (plot_out_spec$plot_progress){
      log_like_store <- as.matrix(store_log_like)
      plot_out <- general_gp_plot_packed(for_like,plot_out_spec,fix_pred)[[1]]
      plot(log_like_store,type="l",ylab="Log-Likelihood",xlab="Iteration",cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4])
    }
    
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