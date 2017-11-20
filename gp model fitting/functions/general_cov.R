general_cov <- function(a1,a2,a3,b1,b2,b3,GG,dGdG,d2Gd2G,dGG,GdG,d2GG,Gd2G,d2GdG,dGd2G){
  value <- a1*b1*GG + a2*b2*dGdG + a3*b3*d2Gd2G + a2*b1*dGG + a1*b2*GdG + a1*b3*Gd2G + a3*b1*d2GG + a2*b3*dGd2G + a3*b2*d2GdG
  return(value)
}