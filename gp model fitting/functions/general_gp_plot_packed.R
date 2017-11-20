general_gp_plot_packed <- function(for_like,plot_out_spec,fix_pred=NULL){
  
  # Unpack
  time=for_like$time;y=for_like$y;sd=for_like$sd;cov_computer=for_like$cov_computer;model_paras=for_like$model_paras;
  cov_paras=for_like$cov_paras;extra_cov_paras=for_like$extra_cov_paras;extra_cov_fun=for_like$extra_cov_fun;
  index_extra_cov=for_like$index_extra_cov;add.nugget=for_like$add.nugget;log_nugget_vec=for_like$log_nugget_vec;
  max_values=for_like$max_values;mean_paras=for_like$mean_paras;mean_spec=for_like$mean_spec;sub_index=for_like$sub_index

  num_outputs <- length(model_paras)/4
  mu.plot <- mean_model(time,y,max_values,mean_paras,mean_spec,num_outputs)
  
  # Prediction points
  t_high <- max(time) 
  t_low <- min(time) 
  test.pts <- seq(t_low, t_high, length.out=801) 
  
  # Covariance between prediction points
  COV.pred.pred <- cov_computer(test.pts,test.pts,model_paras,cov_paras,extra_cov_paras,extra_cov_fun,index_extra_cov)
  
  # Inverse of observaed data covariance matrix 
  var.plot <- sd^2 
  COV.obs.obs <- cov_computer(time,time,model_paras,cov_paras,extra_cov_paras,extra_cov_fun,index_extra_cov) + diag(var.plot)
  if (add.nugget){
    nugget_vec <- c()
    for (i in 1:num_outputs){
      nugget_vec <- c(nugget_vec,rep(exp(log_nugget_vec[i]),length(time)))
    }
    nugget_mat <- diag(nugget_vec)
    COV.obs.obs <- COV.obs.obs + nugget_mat
  }
  COV.obs.obs.inverse <-  solve(COV.obs.obs)
  
  # Covariance between prediction and observed points
  COV.pred.obs <- cov_computer(test.pts, time,model_paras,cov_paras,extra_cov_paras,extra_cov_fun,index_extra_cov)
  COV.obs.pred <- cov_computer(time, test.pts,model_paras,cov_paras,extra_cov_paras,extra_cov_fun,index_extra_cov)
  
  if (is.null(sub_index)){
    sub_index <- 1:length(y)
    pred_index <- 1:(num_outputs*length(test.pts))
  } else {
    start_indexes <- seq(1,length(y),length(time))
    output_indexes <- which(is.element(start_indexes,sub_index[which(sub_index%%length(time)==1)]))
    pred_index <- c()
    for (i in 1:length(output_indexes)){
      pred_index <- c(pred_index,((output_indexes[i]-1)*length(test.pts)+1):(output_indexes[i]*length(test.pts)))
    }
    num_outputs <-  length(output_indexes)
    y <- y[sub_index]
    sd <- sd[sub_index]
    mu.plot <- mu.plot[sub_index]
  }
  
  # Predicted values
  Ef <- COV.pred.obs[pred_index,sub_index] %*% COV.obs.obs.inverse[sub_index,sub_index] %*% (y-mu.plot)
  
  # Special predictions (for cross validation)
  if (is.null(fix_pred)==FALSE){
    if (fix_pred$indicator==TRUE){
      pred_index_special <- 1:(num_outputs*length(fix_pred$t))
      special.COV.pred.obs <- cov_computer(fix_pred$t,time,model_paras,cov_paras,extra_cov_paras,extra_cov_fun,index_extra_cov)
      special.Ef <- special.COV.pred.obs[pred_index_special,sub_index] %*% COV.obs.obs.inverse[sub_index,sub_index] %*% (y-mu.plot)
      if (fix_pred$temp_only == TRUE){
        fix_pred$indicator <- FALSE
        mu_mat <- matrix(mu.plot,ncol=num_outputs)
        special.Ef <- matrix(special.Ef,ncol=num_outputs)+mu_mat
      } else {
        # NOTE: Assumes constant values and shorter than observed data
        mu_mat <- c(matrix(mu.plot,ncol=num_outputs)[1:length(fix_pred$t),])
        special.Ef <- special.Ef + mu_mat
      }
    } else {
      special.Ef <- NULL
    }
  } else {
    special.Ef <- NULL
  }
  
  # Covaraince matrix for mean function
  Cf <- COV.pred.pred[pred_index,pred_index] - COV.pred.obs[pred_index,sub_index]  %*% COV.obs.obs.inverse[sub_index,sub_index] %*% COV.obs.pred[sub_index,pred_index]
  
  # Put in matrix and get error bars
  pred_by_curve <- matrix(Ef,length(test.pts),num_outputs)
  mu_mat <- matrix(mu.plot,ncol=num_outputs)
  for (i in 1:num_outputs){
    pred_by_curve[,i] <- spline(x=time,y=mu_mat[,i],xout=test.pts)$y + pred_by_curve[,i]
  }
  diag.all <- diag(Cf)
  diag.all[diag.all<0] <- 0
  error_by_curve <- matrix(diag.all,length(test.pts),num_outputs)
  
  # Plotting
  if (plot_out_spec$save_plot | plot_out_spec$plot_progress){
    # Plot settings
    if (is.null(plot_out_spec)){
      mag_levels <- rep(1,4)
    } else {
      mag_levels <- plot_out_spec$plot_settings$mag_levels
    }
    
    # Plot
    if (plot_out_spec$plot_log_like){
      plot_num <- num_outputs+1
    } else {
      plot_num <- num_outputs
    }
    par(mfrow=c(ceiling((plot_num)/ceiling(sqrt(plot_num))),ceiling(sqrt(plot_num)))) # Extra one is for log likelihood plot
    
    data.now <- matrix(y,ncol=num_outputs)
    sd.mat.up <- matrix(sd,ncol=num_outputs)
    for (i in 1:num_outputs){
      pred.now <- pred_by_curve[,i]
      error.now <- error_by_curve[,i]
      y.now <- data.now[,i]
      sd.now <- sd.mat.up[,i]
      y_region <- c(rev(pred.now-2*sqrt(error.now)), pred.now+2*sqrt(error.now))
      ylims <- c(min(c(pred.now,y.now))-3*max(c(abs(sd.now),abs(error.now))),max(c(pred.now,y.now))+3*max(c(abs(sd.now),abs(error.now))))
      if (length(plot_out_spec$xlim_plot)==0){
        if (i==1){
          if (mean_spec$planet_indicator){
            ylab1 <- "Corrupted RV (m/s)"
          } else {
            ylab1 <- "Apparent RV (m/s)"
          }
          plot(c(t_low,t_high),c(min(c(y.now,pred.now)),max(c(y.now,pred.now))),ylim=ylims,col="red",type="n",xlab="Time",ylab=ylab1,cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4])
        } else {
          plot(c(t_low,t_high),c(min(c(y.now,pred.now)),max(c(y.now,pred.now))),ylim=ylims,col="red",type="n",xlab="Time",ylab=paste(output_name,i-1," score",sep=""),cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4])
        }
      } else {
        if (i==1){
          if (mean_spec$planet_indicator){
            ylab1 <- "Corrupted RV (m/s)"
          } else {
            ylab1 <- "Apparent RV (m/s)"
          }
          plot(c(t_low,t_high),c(min(c(y.now,pred.now)),max(c(y.now,pred.now))),xlim=plot_out_spec$xlim_plot,ylim=ylims,col="red",type="n",xlab="Time",ylab=yab1,cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4])        
        } else {
          plot(c(t_low,t_high),c(min(c(y.now,pred.now)),max(c(y.now,pred.now))),xlim=plot_out_spec$xlim_plot,ylim=ylims,col="red",type="n",xlab="Time",ylab=paste(output_name,i-1," score",sep=""),cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4])
        }
      }
      polygon(c(rev(test.pts), test.pts), y_region, col = 'grey48', border = NA)
      points(time, y.now)
      segments(time, y.now-2*sd.now,time, y.now+2*sd.now)
      epsilon = 0.1
      segments(time-epsilon,y.now-2*sd.now,time+epsilon,y.now-2*sd.now)
      segments(time-epsilon,y.now+2*sd.now,time+epsilon,y.now+2*sd.now)
      if (is.null(fix_pred)==FALSE){
        if (fix_pred$indicator==TRUE){
          points(fix_pred$t,fix_pred$y[,i],pch=16,col=3)
          segments(fix_pred$t, fix_pred$y[,i]-2*fix_pred$sd[,i],fix_pred$t, fix_pred$y[,i]+2*fix_pred$sd[,i],col=3)
          segments(fix_pred$t-epsilon,fix_pred$y[,i]-2*fix_pred$sd[,i],fix_pred$t+epsilon,fix_pred$y[,i]-2*fix_pred$sd[,i],col=3)
          segments(fix_pred$t-epsilon,fix_pred$y[,i]+2*fix_pred$sd[,i],fix_pred$t+epsilon,fix_pred$y[,i]+2*fix_pred$sd[,i],col=3)
        }
      }
      lines(test.pts,pred.now,col="red")
      
    }
  }
  
  
  return(list(test.pts,Ef,Cf,special.Ef))
  
}
