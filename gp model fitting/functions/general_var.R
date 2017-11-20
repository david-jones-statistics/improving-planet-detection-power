general_var <- function(a1,a2,a3,a4,GG,dGdG,d2Gd2G,dGG,GdG,d2GG,Gd2G,d2GdG,dGd2G,plus_cov){
  
  # For testing
  # a1 <- coeffs[i,1]
  # a2 <- coeffs[i,2]
  # a3 <- coeffs[i,3]
  # a4 <- coeffs[i,4]
  # GG <- cov_GG
  # dGdG <- cov_dGdG
  # d2Gd2G <- cov_d2Gd2G
  # dGG <- cov_dGG
  # GdG <- cov_GdG
  # d2GG <- cov_d2GG
  # Gd2G <- cov_Gd2G
  # d2GdG <- cov_d2GdG
  # dGd2G <- cov_dGd2G
  # plus_cov <- extra_cov_now

  
  value <- a1^2*GG + a2^2*dGdG + a3^2*d2Gd2G + a1*a2*(dGG + GdG) + a1*a3*(Gd2G + d2GG) + a2*a3*(dGd2G + d2GdG) + a4^2*plus_cov
  return(value)
}