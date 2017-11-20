fit_activity_gp_model <- function(data,model_spec,mean_spec,select_cov="quasi_periodic_",custom_cov_spec=NULL,plot_out_spec=NULL,missing.data.index=NULL,initial_value_spec=NULL,cross.validation=FALSE){
  
  # Testing 
  if (1==2){
    select_cov="expos_"
    select_cov="quasi_periodic_"
    select_cov="matern_"
    custom_cov_spec=NULL
    missing.data.index <- c()  #missing.data.index.store[[1]]missing.data.index.store[[data_loop]]
  }
  
  stopifnot(dim(data$y)[2]==dim(model_spec$model_mat)[1])
  #num_outputs <- dim(data$y)[2]
  
  # Set covariance function
  if (is.null(custom_cov_spec)){
    cov_spec <- get_preset_cov_funs(select_cov,model_spec$opt_stellar_period,no_expos=1)
  } else {
    cov_spec <- custom_cov_spec
    source("cov_spec$functions")
  }
  
  # Set covariance function for addtional independent GPs
  # Here its just assumed to be the same as the above covariance function
  # Could easily be modified to allow a differnt covariance function
  # Slightly more work would be required to allow 
  # multiple covariance functions for the independent GPs
  extra_cov_spec <- cov_spec # Only need parameter numbers, indexes, and names
  extra_cov_spec$extra_cov_fun <- cov_fun
  
  # Set build function - uses the derviative functions loaded above
  cov_spec$cov_computer <- build_full_cov_mat
  
  # Extract data amd normalize for numerical stability and then vectorize
  stellar_period <- data$stellar_period
  no_phases <- data$no_phases # number of phases i.e. 25 in our case
  time <- data$time
  if (num_outputs>1){
    mid_point_vec <- apply(apply(data$y,2,range),2,mean)
  } else {
    mid_point_vec <- mean(range(data$y))
  }
  y <- data$y-t(matrix(rep(mid_point_vec,length(time)),nrow=num_outputs))
  sd <- data$sd
  means <- data$means-mid_point_vec
  n <- dim(y)[1]
  if (num_outputs>1){
    max_values <- apply(abs(y),2,max)
  } else {
    max_values <- max(abs(y))
  }
  if (num_outputs>1){
    for (i in 1:num_outputs){
      y[,i] <- y[,i]/max_values[i]
      sd[,i] <- sd[,i]/max_values[i]
      means[,i] <- means[,i]/max_values[i]
      # NOTE: if using GPCA values covariance matrix instead of just SDs then will need to 
      # normalize its entries 
    }
  } else {
    y <- y/max_values
    sd <- sd/max_values
    means <- means/max_values
  }
  y <- c(y)
  sd <- c(sd)
  means <- c(means)
  
  # Extract missing data for cross validation 
  n_missing <- 0
  time_missing <- NULL
  y_missing <- NULL
  sd_missing <- NULL
  fix_pred <- NULL
  if (is.null(missing.data.index)==FALSE & cross.validation==TRUE){
    if (length(missing.data.index) > 0){
      n_missing <- data$n_missing
      time_missing <- data$time_missing
      y_missing <- data$y_missing -t(matrix(rep(mid_point_vec,length(time_missing)),nrow=num_outputs))
      sd_missing <- data$sd_missing
      for (i in 1:num_outputs){
        y_missing[,i] <- y_missing[,i]/max_values[i]
        sd_missing[,i] <- sd_missing[,i]/max_values[i]
      }
      fix_pred <- list()
      fix_pred$t <- time_missing
      fix_pred$y <- y_missing 
      fix_pred$sd <- sd_missing
      fix_pred$indicator <- TRUE
      fix_pred$temp_only <- FALSE
    }
  } 
  
  # Compute settings - 100s mean we don't care about those parameters only the likelihood
  max_iter <- 100
  min_iter <- 5 #28
  # Annealing-type optimization set-up (sd inflated to begin with)
  #min_final_iters <- 10
  #num_per_step <- 3
  #num_sd_steps <- (min_iter-min_final_iters)/num_per_step
  #temp_inflate_sd <- c(t(matrix(rep(rev(seq((1/num_sd_steps)*0.1/min(sd),0.1/min(sd),length.out=num_sd_steps)),num_per_step),num_sd_steps)))
  #temp_inflate_sd <- c(t(matrix(rep(rev(seq(2,10,length.out=num_sd_steps)),num_per_step),num_sd_steps)))
  inflate_sd <- rep(1,max_iter) #c(temp_inflate_sd,rep(1,max_iter-min_iter+min_final_iters)) # 
  cov_para_bounds <- c()
  compute_settings <- list(log_like_tol=0.001,model_tol=100,cov_tol=100,extra_cov_tol=100,mean_tol=100,log_nugget_tol=100,
                           min_iter=min_iter,max_iter=max_iter,opt_bound=15,cov_paras_max=20,extra_cov_paras_max=20,
                           inflate_sd=inflate_sd)
  
  # Optimization
  # If there is only a single run then the initialization sepcified in the run file e.g. example_run_code.R
  # If multiple multi_initializations==TRUE then default random initializations will be used. 
  # Multiple initializations are often needed to ensure the global optimum is found. 
  if (initial_value_spec$multi_initializations==TRUE){
    opt_store <- list()
    max_log_like_store <- numeric(initial_value_spec$num_initializations)
    periodic <- cov_spec$periodic
    initial_value_length <- length(cov_spec$opt_index_cov_paras)+ifelse(periodic,1,0)
    for (initial_value_num in 1:initial_value_spec$num_initializations){
      
      print(paste("Initialization: ",initial_value_num,sep=""))
      
      # Initialize
      opt_index <- mean_spec$opt_index
      mean_initial <- runif(length(opt_index),min=mean_spec$mean_paras_bounds[[1]][opt_index],max=mean_spec$mean_paras_bounds[[2]][opt_index])
      # More restrictive initialization to stop numerical errors arising from extreme disagreement with the data
      if (mean_spec$planet_indicator){
        mean_initial[1] <- runif(1,min=mean_spec$mean_paras_bounds[[1]][1],max=log(0.25*rv_range/2))
      }
      # final initialization vector
      initial_value <- c(runif(3,c(2,-1,1),c(2.5,0,2)),runif(3,c(2,-1,1),c(2.5,0,2)),mean_initial)
      #mean_initial <- runif(length(opt_index),min=mean_spec$mean_paras_bounds[[1]],max=mean_spec$mean_paras_bounds[[2]])
      #initial_value_spec$initial_value <- c(log(data$stellar_period),runif(2,-2,2),log(data$stellar_period),runif(2,-2,2),mean_initial)
      #initial_value <- c(initial_value,mean_initial)
      
      # Optimize
      opt_out <- optimize_parameters(time,y,sd,initial_value,model_spec,cov_spec,extra_cov_spec,compute_settings,plot_out_spec,fix_pred,max_values,mean_spec)
      opt_store[[initial_value_num]] <- opt_out
      max_log_like_store[initial_value_num] <- max(as.numeric(opt_out$log_like_store))
    }
    best_index <- which.max(max_log_like_store)
    opt_out <- opt_store[[best_index]]
  } else {
    initial_value <- initial_value_spec$initial_value
    opt_out <- optimize_parameters(time,y,sd,initial_value,model_spec,cov_spec,extra_cov_spec,compute_settings,plot_out_spec,fix_pred,max_values,mean_spec)
  }
  
  # Cross validation (constant mean assumed)
  cv_output <- NULL
  if (is.null(missing.data.index)==FALSE & cross.validation==TRUE){
    if (length(missing.data.index) > 0){
      
      # 1) Full data log-like
      for_like <- opt_out$for_like
      for_like$time <- c(time,time_missing)
      for_like$y <- c(rbind(matrix(y,ncol=3),y_missing)) # c(y,c(y_missing))
      for_like$sd <- c(rbind(matrix(sd,ncol=3),sd_missing)) #c(sd,c(sd_missing))
      like_fun2 <- general_gp_log_lik_packed
      cv_full_neg_log_like <- -like_fun2(for_like)[[1]]
      
      # 2) Missing data log-like
      for_like$time <- time_missing
      for_like$y <- c(y_missing)
      for_like$sd <- c(sd_missing)
      cv_neg_log_like <- -like_fun2(for_like)[[1]]
      
      # 3) Missing data sum of squared residuals 
      index_extra_cov <- which(model_mat[,4]!=0)
      plot_out <- general_gp_plot_packed(opt_out$for_like,plot_out_spec,fix_pred)
      pred_errors <- plot_out[[4]] - c(y_missing)
      cv_mean_sr <- mean(pred_errors^2)
      
      # 4) Missing data weighted sum of squared residuals
      cv_weighted_mean_sr <- mean((pred_errors/c(sd_missing))^2)
      
      # 5) Doppler missing data sum of squared residuals
      doppler_pred_errors <- pred_errors[1:n_missing]
      cv_doppler_mean_sr <- mean(doppler_pred_errors^2)
      
      # 6) Doppler missing data weighted sum of squared residuals
      cv_doppler_weighted_mean_sr <- mean((doppler_pred_errors/sd_missing[,1])^2)
      
      # 7) Conditional log-like
      # Initially used max of log_like_store
      log_like_vec <- as.numeric(opt_out$log_like_store)
      cv_cond_neg_log_like <- log_like_vec[length(log_like_vec)] + cv_full_neg_log_like
      
      # Store CV criteria
      cv_output <- list(cv_full_neg_log_like=cv_full_neg_log_like,cv_neg_log_like=cv_neg_log_like,cv_mean_sr=cv_mean_sr,
                        cv_weighted_mean_sr=cv_weighted_mean_sr,cv_doppler_mean_sr=cv_doppler_mean_sr,cv_doppler_weighted_mean_sr=cv_doppler_weighted_mean_sr,
                        pred_errors=pred_errors,doppler_pred_errors=doppler_pred_errors,cv_cond_neg_log_like=cv_cond_neg_log_like)
    }
  }
  
  # Plot fitted model
  # Save option
  if (plot_out_spec$save_plot){
    plot_settings <- plot_out_spec$plot_settings
    filename <- plot_out_spec$filename
    pdf(filename,width=plot_settings$pdf.width,height=plot_settings$pdf.height)
    par(mar=plot_settings$plot_gaps,mgp=plot_settings$axis_gaps,oma=plot_settings$image_gaps)
    plot_out <- general_gp_plot_packed(opt_out$for_like,plot_out_spec,fix_pred)
    mag_levels <- plot_out_spec$plot_settings$mag_levels
    plot(opt_out$log_like_store,type="l",ylab="Log-Likelihood",xlab="Iteration",cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4])
    dev.off()
  }
  if (initial_value_spec$multi_initializations==TRUE & plot_out_spec$plot_progress){
    plot_settings <- plot_out_spec$plot_settings
    par(mar=plot_settings$plot_gaps,mgp=plot_settings$axis_gaps,oma=plot_settings$image_gaps)
    plot_out <- general_gp_plot_packed(opt_out$for_like,plot_out_spec,fix_pred)
    plot(opt_out$log_like_store,type="l",ylab="Log-Likelihood",xlab="Iteration",cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4])
    print("Final fit now plotted")
  }
  
  # Remove unused model components to avoid confusion 
  if (model_spec$add.nugget){
    optim_out <- opt_out
  } else {
    optim_out <- opt_out[-5] 
  }
  if (sum(model_spec$model_mat[,4])>0){
    optim_out <- opt_out
  } else {
    optim_out <- opt_out[-3] 
  }
  
  # Final data after normalization
  final_data <- list(stellar_period=stellar_period,no_phases=no_phases,n=n,time=time,y=matrix(y,ncol=num_outputs),sd=matrix(sd,ncol=num_outputs),
                     n_missing=n_missing,time_missing=time_missing,y_missing=y_missing,sd_missing=sd_missing,max_values=max_values)
  
  value <- list(optim_out=optim_out,final_data=final_data,cov_spec=cov_spec,extra_cov_spec=extra_cov_spec,cv_output=cv_output)
  return(value)
}