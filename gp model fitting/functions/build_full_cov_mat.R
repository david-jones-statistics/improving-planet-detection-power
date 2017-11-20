build_full_cov_mat <- function(t1,t2,model_paras,cov_paras,extra_cov_paras,extra_cov_fun,index_extra_cov){
  
  num_outputs <- length(model_paras)/4

  coeffs <- t(matrix(model_paras,ncol=num_outputs))
  
  # Get main and cross terms
  cov_GG <- cov_make(t1,t2,cov_fun,cov_paras)
  
  if (sum(abs(coeffs[,2:3]))>0){
    cov_GdG <- cov_make(t1,t2,dt_cov_fun,cov_paras) # this is the way the derivative is coded
    cov_dGG <- -cov_GdG 
    cov_dGdG <- cov_make(t1,t2,dtdt_cov_fun,cov_paras)
  } else {
    cov_GdG <- 0
    cov_dGG <- 0
    cov_dGdG <- 0
  }
  
  if (sum(abs(coeffs[,3]))>0){
    cov_Gd2G <- cov_make(t1,t2,dt2_cov_fun,cov_paras)
    cov_d2GG <- cov_Gd2G
    cov_dGd2G <- cov_make(t1,t2,dtdt2_cov_fun,cov_paras)
    cov_d2GdG <- -cov_dGd2G
    cov_d2Gd2G <- cov_make(t1,t2,dt2dt2_cov_fun,cov_paras)
  } else {
    cov_Gd2G <- 0
    cov_d2GG <- 0
    cov_dGd2G <- 0
    cov_d2GdG <- 0
    cov_d2Gd2G <- 0
  }
  
  if (sum(index_extra_cov)>0){
    extra_cov <- cov_make(t1,t2,extra_cov_fun,extra_cov_paras)
  }
   
  # Build covariance matrix
  len1 <- length(t1)
  len2 <- length(t2)
  mat <- matrix(NA,num_outputs*len1,num_outputs*len2)
  for (i in 1:num_outputs){
    if (is.element(i,index_extra_cov)){
      extra_cov_now <- extra_cov
    } else {
      extra_cov_now <- 0
    }
    mat[((i-1)*len1+1):(i*len1),((i-1)*len2+1):(i*len2)] <- general_var(coeffs[i,1],coeffs[i,2],coeffs[i,3],coeffs[i,4],
                                                                        cov_GG,cov_dGdG,cov_d2Gd2G,cov_dGG,cov_GdG,cov_d2GG,cov_Gd2G,cov_d2GdG,cov_dGd2G,extra_cov_now)
  }
  for (i in 1:num_outputs){
    for (j in setdiff(1:num_outputs,i)){
      mat[((i-1)*len1+1):(i*len1),((j-1)*len2+1):(j*len2)] <- general_cov(coeffs[i,1],coeffs[i,2],coeffs[i,3],
                                                                          coeffs[j,1],coeffs[j,2],coeffs[j,3],
                                                                          cov_GG,cov_dGdG,cov_d2Gd2G,cov_dGG,cov_GdG,cov_d2GG,cov_Gd2G,cov_d2GdG,cov_dGd2G)
    }
  }
  
  return(mat)
  
}