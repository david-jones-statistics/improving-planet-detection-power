planet_optimizer_reuse_mats2 <- function(para_type,para_opt_index,period_index,ang1_index,log_normalizing,shrink,bounds,for_like,rv_errs){
  
  if (1==2){
    para_type <- "mean_paras"
    period_index <- 3
    ang1_index <- 2
    bounds <- mean_paras_bounds
    para_opt_index <- key_planet_paras_index
    #for_like_store <- for_like
    #for_like <- for_like_store
  }
  
  para_type_index <- which(names(for_like)==para_type)
  
  # Search vectors
  period_vec <- log(seq(exp(bounds[[1]][period_index]),exp(bounds[[2]][period_index]),0.1))
  angle_search_vec <- seq(0,2*pi,0.1)
  
  # Planet paras initialize
  mean_paras <- for_like$mean_paras
  mean_spec <- for_like$mean_spec
  mp_len <- length(mean_paras)
  kepler_paras <- as.list(mean_paras[1:(mp_len-(num_outputs-1))])
  for (k in mean_spec$log_scale_index){
    kepler_paras[[k]] <- exp(kepler_paras[[k]])
  }
  for (k in mean_spec$logit_scale_index){
    kepler_paras[[k]] <- exp(kepler_paras[[k]])/(1+exp(kepler_paras[[k]]))
  }
  names(kepler_paras) <- mean_spec$planet_para_names
  
  # Non-RV means 
  num_outputs <- for_like$model_spec$num_outputs
  non_rv_mu_vec <- c(t(matrix(rep(mean_paras[(mp_len-(num_outputs-2)):mp_len],length(for_like$time)),num_outputs-1,length(for_like$time))))
  
  # Optimize period
  lper <- length(period_vec)
  lang <- length(angle_search_vec)
  like_by_grid <- numeric(lper*lang)
  
  # 3 dimensional grid
  for (per_num in 1:lper){
    print(paste("Percentage complete: ",round(100*per_num/lper),sep=""))
    for (ang1_num in 1:lang){
      kepler_paras[[period_index]] <- exp(period_vec[per_num])
      kepler_paras[[ang1_index]] <- angle_search_vec[ang1_num]
      like_by_grid[(per_num-1)*lang + ang1_num] <- planet_opt_function_reuse_mats2(log_normalizing,shrink,kepler_paras,non_rv_mu_vec,for_like,rv_errs)
    }
  }
  
  # Identify highest likelihood
  best_paras_index <- which.min(like_by_grid)
  pindex <- floor(best_paras_index/lang)
  for_like[[para_type_index]][period_index] <- period_vec[pindex+1]
  a1_index <- best_paras_index - pindex*lang
  for_like[[para_type_index]][ang1_index] <- angle_search_vec[a1_index]
  kepler_paras[[period_index]] <- exp(for_like[[para_type_index]][period_index])
  kepler_paras[[ang1_index]] <- for_like[[para_type_index]][ang1_index]
  planet_rv <- v_regress(rv_errs,for_like$time,kepler_paras)
  kepler_paras$gamma <- planet_rv[[2]]$coefficients[1]
  kepler_paras$K <- sqrt(sum(planet_rv[[2]]$coefficients[2:3]^2))
  kepler_paras$w <- atan(-(planet_rv[[2]]$coefficients[3]/planet_rv[[2]]$coefficients[2]))
  for (k in 1:(mp_len-(num_outputs-1))){
    mean_paras[k] <- kepler_paras[[k]]
  }
  for (k in mean_spec$log_scale_index){
    mean_paras[k] <- log(mean_paras[k])
  }
  for (k in mean_spec$logit_scale_index){
    mean_paras[k] <- log(mean_paras[k]/(1-mean_paras[k]))
  }
  for_like$mean_paras <- mean_paras
  
  if (abs(max(-like_by_grid)-like_fun2(for_like)[[1]])>10^-4){
    mean_paras[6] <- mean_paras[6] + pi
    for_like$mean_paras <- mean_paras
    if (abs(max(-like_by_grid)-like_fun2(for_like)[[1]])>10^-4){
      print("WARNING GRID MESSED UP!")
    }
  }
  
  # Run normal optimizer
  initial_guess <- for_like[[para_type_index]][para_opt_index]
  out <- optim(initial_guess,opt_function_reuse_mats,para_type=para_type,para_opt_index=para_opt_index,
               log_normalizing=log_normalizing,shrink=shrink,for_like=for_like,
               lower=bounds[[1]][para_opt_index],upper=bounds[[2]][para_opt_index],method="L-BFGS-B")
  
  return(list(for_like=for_like,period_vec=period_vec,angle_search_vec=angle_search_vec,out=out,grid_max_log_like=-min(like_by_grid)))
}