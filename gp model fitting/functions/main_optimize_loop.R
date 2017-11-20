main_optimize_loop <- function(for_like,compute_settings,plot_out_spec,inflate_sd,opt_index_model_paras,opt_index_cov_paras,opt_index_extra_cov_paras,key_planet_paras_index,other_paras_index,cov_para_bounds,extra_cov_para_bounds,mean_paras_bounds,stellar_period_vec,log_per_guesses,index_extra_cov,opt.nugget,opt.mean,fix_pred,location_pars_index){

  # For annealing-type approach
  true_sd <- for_like$sd
  
  true_mean_paras_bounds <- mean_paras_bounds
  
  # Initialize
  num_outputs <- for_like$model_spec$num_outputs
  count <- 1
  log_like_diff <- 100
  store_log_like <- list()
  store_log_like[[1]] <- like_fun2(for_like)[[1]]

  # Main loop
  while ((count < compute_settings$min_iter | log_like_diff > compute_settings$log_like_tol) &
         count < compute_settings$max_iter){
    
    # if (count <= compute_settings$min_iter/2 & length(key_planet_paras_index) > 0){
    #   mean_paras_bounds[[1]][1] <- log(0.1*rv_range/2)
    #   mean_paras_bounds[[1]][2] <- log(0.75*rv_range/2)
    # } else {
    #   mean_paras_bounds <- true_mean_paras_bounds
    # }
    
    # Annealing-type approach: inflate standard deviations
    for_like$sd <- true_sd*compute_settings$inflate_sd[count]
    
    count <- count + 1 
    
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
      period_out <- opt_period_batch("extra_cov_paras",basic_period_opt,opt_index_extra_cov_paras,1,stellar_period_vec,plot_out_spec,extra_cov_para_bounds,for_like)
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
    if (length(other_paras_index)==0){
      fix_pred_temp <- list(t=for_like$time,y=for_like$y,sd=for_like$sd,indicator=TRUE,temp_only=TRUE)
      plot_out <- general_gp_plot_packed(for_like,plot_out_spec,fix_pred_temp)
      if (log_like_diff > 0.5){
        if (num_outputs > 1){
          curr_means <- apply(plot_out[[4]],2,mean)
          data_means <- apply(matrix(for_like$y,ncol=num_outputs),2,mean)
          for_like$mean_paras[location_paras_index] <- for_like$mean_paras[location_paras_index] + (data_means-curr_means)
        } else {
          curr_means <- mean(plot_out[[4]])
          data_means <- mean(for_like$y)
          for_like$mean_paras[location_paras_index] <- for_like$mean_paras[location_paras_index] + (data_means-curr_means)
        }
      } else {
        other_paras_index <- location_paras_index
      }
    }
    
    if (opt.mean){
      basic_period_opt <- FALSE
      key_quantities <- like_fun2(for_like)[[2]]
      shrink <- key_quantities$shrink
      log_normalizing <- key_quantities$log_normalizing
      if (length(key_planet_paras_index) > 0){
        period_out <- opt_period_batch_reuse_mats("mean_paras",basic_period_opt,key_planet_paras_index,3,log_per_guesses,log_normalizing,shrink,plot_out_spec,mean_paras_bounds,for_like)
        for_like <- period_out[[1]]
      }
      if (length(other_paras_index)>=1){
        for_like <-  opt_batch_reuse_mats("mean_paras",other_paras_index,plot_out_spec,mean_paras_bounds,for_like,log_normalizing,shrink)
      }
      if (length(key_planet_paras_index) > 0){
        print(for_like$mean_paras)
        print(for_like$cov_paras)
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
  
  return(list(for_like=for_like,store_log_like=store_log_like))
  
}