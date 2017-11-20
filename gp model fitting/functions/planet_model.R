M <- function(t,kepler_paras){
  tau <- kepler_paras$tau
  M0 <- kepler_paras$M0
  value <- 2*pi*t/tau + M0
  return(value)
}

E <- function(t,kepler_paras){
  e <- kepler_paras$e
  Evalue <- matrix(NA,length(t),1)
  for (j in 1:length(t)){
    Evalue[j] <- uniroot(function(x){x-e*sin(x)-M(t[j],kepler_paras)},lower=M(t[j],kepler_paras)-1-abs(e),upper=abs(e)+M(t[j],kepler_paras)+1)$root
  }
  return(Evalue)
}

phi <- function(t,kepler_paras){
  e <- kepler_paras$e
  value <- 2*atan(((1+e)/(1-e))*tan(E(t,kepler_paras)/2))
  return(value)
}

v <- function(t,kepler_paras){
  e <- kepler_paras$e
  K <- kepler_paras$K
  gamma <- kepler_paras$gamma
  w <- kepler_paras$w
  value <- K*(e*cos(w) + cos(w+phi(t,kepler_paras))) + gamma
  return(value)
}

v_regress <- function(resids,t,kepler_paras){
  e <- kepler_paras$e
  x1 <- e + cos(phi(t,kepler_paras))
  x2 <- sin(phi(t,kepler_paras))
  resid_planet_model <- lm(resids ~ x1+x2)
  value <- resid_planet_model$fitted
  return(list(value,resid_planet_model))
}