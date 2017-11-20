get_data <- function(data_filename,n,num_outputs,missing.data.index=NULL,plot.data=TRUE,make_cov_mat=FALSE){
  
  true_planet_paras <- NULL
  load(data_filename)
  
  time <- tvec[1:n] 
  values <- activity_data[1:n,1:num_outputs]
  values_sd <- activity_sd[1:n,1:num_outputs]
  
  # In stead of just the SDs we could work with the covariance matrix for the 
  # observations.
  # Question: are observations correlated with all other measurements made at that phase
  # or only those measurements made at the same time? 
  # NOTE if we use this matrix we need to add normalization to fit_Activity_gp_model.R 
  if (make_cov_mat){
    no_phases <- no_phases*n/length(tvec)
    values_cov_mat <- matrix(0,num_outputs*n,num_outputs*n)
    diag(values_cov_mat) <- c(values_sd)^2
    for (i in 1:n){
      phase <- i - floor((i-1)/no_phases)*no_phases
      for (j in 1:num_outputs){
        for (k in setdiff(1:num_outputs,j)){
          values_cov_mat[(j-1)*n+i,(k-1)*n+i] <- cov_mat_store[[phase]][j,k]
        }
      }
    }
  } else {
    values_cov_mat <- NULL
    no_phases <- NULL
  }

  
  # Delete observaitons if set to missing 
  if (is.null(missing.data.index)==FALSE){
    if (length(missing.data.index) > 0){
      n_missing <- length(missing.data.index)
      n <- n - length(missing.data.index)
      time_missing <- time[missing.data.index]
      time <- time[-missing.data.index]
      values_missing <- values[missing.data.index,]
      values <- values[-missing.data.index,]
      values_sd_missing <- values_sd[missing.data.index,]
      values_sd <- values_sd[-missing.data.index,]
      means <- t(matrix(rep(apply(values[-missing.data.index,],2,mean),n),num_outputs,n))
      index <- rep(missing.data.index,num_outputs) + c(sapply(seq(0,(num_outputs-1)*n,n),rep_len,length.out=length(missing.data.index)))
      if (make_cov_mat){
        values_cov_mat <- values_cov_mat[-index,-index]
      }
    } else {
      means <- t(matrix(rep(apply(values,2,mean),n),num_outputs,n))
      n_missing <- 0
      time_missing <- NULL
      values_missing <- NULL
      values_sd_missing <- NULL
    }
  } else {
    if (num_outputs > 1){
      means <- t(matrix(rep(apply(values,2,mean),n),num_outputs,n))
    } else { 
      means <- mean(values)
    }
    n_missing <- 0
    time_missing <- NULL
    values_missing <- NULL
    values_sd_missing <- NULL
  }
  
  if (plot.data==TRUE){
    plot_num <- num_outputs
    pdf.width <- 10
    pdf.height <- 7
    plot_gaps <- c(6,5,0,2) # bottom, left, top, right
    image_gaps <- c(1,1,1,1) # bottom, left, top, right
    axis_gaps <- c(3,1,0) # label, number, tick (?)
    mag_levels <- c(1.5,1.75,1.5,2,1)  # main, lab, axis, line, pts
    if (save.data.plot){
      pdf(save.data.filename,width=pdf.width,height=pdf.height)
      par(mar=plot_gaps,mgp=axis_gaps,oma=image_gaps,mfrow=c(ceiling(sqrt(plot_num)),ceiling((plot_num)/ceiling(sqrt(plot_num)))))
    } else {
      par(mfrow=c(ceiling(sqrt(plot_num)),ceiling((plot_num)/ceiling(sqrt(plot_num)))))
    }
    if (output_name=="gpca" | output_name=="gpca_hd"){
      ylabel <- "PC"
    } else {
      ylabel <- "DM"
    }
    for (i in 1:num_outputs){
      if (num_outputs > 1){
        if (i == 1){
          ylabel.now <- "Apparent RV (m/s)"
        } else {
          ylabel.now <- paste(ylabel,i-1,sep="")
        }
        y.bar.limits <- c(values[,i]-2*values_sd[,i],values[,i]+2*values_sd[,i])
        plot(time,values[,i],ylab=ylabel.now,ylim=range(y.bar.limits),xlab="Days",pch=16,col=4,cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4],cex=mag_levels[5])
        segments(time, values[,i]-2*values_sd[,i],time, values[,i]+2*values_sd[,i],col=4,lwd=mag_levels[4])
        epsilon = 0.25
        segments(time-epsilon,values[,i]-2*values_sd[,i],time+epsilon,values[,i]-2*values_sd[,i],col=4,lwd=mag_levels[4])
        segments(time-epsilon,values[,i]+2*values_sd[,i],time+epsilon,values[,i]+2*values_sd[,i],col=4,lwd=mag_levels[4])
      } else {
        plot(time,values,"Apparent RV (m/s)",xlab="Days",pch=16,col=4,cex.main=mag_levels[1],cex.lab=mag_levels[2],cex.axis=mag_levels[3],lwd=mag_levels[4],cex=mag_levels[5])
        segments(time, values-2*values_sd,time, values+2*values_sd,col=4,lwd=mag_levels[4])
        epsilon = 0.25
        segments(time-epsilon,values-2*values_sd,time+epsilon,values-2*values_sd,col=4,lwd=mag_levels[4])
        segments(time-epsilon,values+2*values_sd,time+epsilon,values+2*values_sd,col=4,lwd=mag_levels[4])
      }
    }
    if (save.data.plot){
      dev.off()
    }
  }
  
  data <- list(stellar_period=stellar_period,no_phases=no_phases,n=n,time=time,y=values,sd=values_sd,means=means,values_cov_mat=values_cov_mat,
               n_missing=n_missing,time_missing=time_missing,y_missing=values_missing,sd_missing=values_sd_missing,true_planet_paras=true_planet_paras)
  
  return(data)
}