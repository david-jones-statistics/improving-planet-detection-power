mean_model <- function(t,y,max_values,mean_paras,mean_spec,num_outputs){
  
  planet_indicator <- mean_spec$planet_indicator
  
  if (planet_indicator){
    mp_len <- length(mean_paras)
    kepler_paras <- as.list(mean_paras[1:(mp_len-(num_outputs-1))])
    for (k in mean_spec$log_scale_index){
      kepler_paras[[k]] <- exp(kepler_paras[[k]])
    }
    for (k in mean_spec$logit_scale_index){
      kepler_paras[[k]] <- exp(kepler_paras[[k]])/(1+exp(kepler_paras[[k]]))
    }
    names(kepler_paras) <- mean_spec$planet_para_names
    planet_rv <- v(t,kepler_paras) #/max_values[1] # Normalize planet signal as for the stellar activity
    value <- c(planet_rv,c(t(matrix(rep(mean_paras[(mp_len-(num_outputs-2)):mp_len],length(t)),num_outputs-1,length(t)))))
  } else {
    value <- c(t(matrix(rep(mean_paras,length(t)),num_outputs,length(t))))
  }
  
  return(value)
  
}