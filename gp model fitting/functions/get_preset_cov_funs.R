get_preset_cov_funs <- function(cov_choice="periodic_and_expo_",opt_stellar_period,no_expos=2){
  
  store_dir <- getwd()
  
  if (cov_choice == "general_matern_"){
    setwd("./functions/general_matern/")
    num_cov_paras <- 1
    periodic <- FALSE
    paras_names <- c("rho")
  }
  if (cov_choice == "matern5by2_"){
    setwd("./functions/matern5by2/")
    num_cov_paras <- 1
    periodic <- FALSE
    paras_names <- c("rho")
  }
  if (cov_choice == "matern_"){
    setwd("./functions/matern9by2/")
    num_cov_paras <- 1
    periodic <- FALSE
    paras_names <- c("rho")
  }
  if (cov_choice == "quasi_periodic_"){
    setwd("./functions/quasi_periodic_cov")
    num_cov_paras <- 3
    periodic <- TRUE
    paras_names <- c("log.period","log.lambda.p","log.lambda.e")
  }
  if (cov_choice == "periodic_"){
    setwd("./functions/periodic_cov/")
    num_cov_paras <- 2
    periodic <- TRUE
    paras_names <- c("log.period","log.lambda.p")
  }
  if (cov_choice == "cosine_"){
    setwd("./functions/cos_cov/")
    num_cov_paras <- 1
    periodic <- TRUE
    paras_names <- c("log.period")
  }
  if (cov_choice == "geometric_"){
    setwd("./functions/geometric_cov/")
    num_cov_paras <- 2
    periodic <- TRUE
    paras_names <- c("log.period","log.factor")
  }
  if (cov_choice == "expos_"){
    setwd("./functions/squared_expos_cov/")
    num_cov_paras <- 2*no_expos
    periodic <- FALSE
    paras_names <- character(num_cov_paras)
    for (i in 1:no_expos){
      paras_names[(i-1)*2+1] <- paste("log.beta",i,sep="")
      paras_names[2*i] <- paste("log.lambda",i,sep="")
    }
  }

  if (periodic==TRUE & opt_stellar_period==FALSE){
    if (num_cov_paras > 1){
      opt_index_cov_paras <- 2:num_cov_paras
    } else {
      opt_index_cov_paras <- c()
    } 
  } else {
    opt_index_cov_paras <- 1:num_cov_paras
  }
  num_cov_paras_opt <- length(opt_index_cov_paras)
  
  source("cov_fun.R")
  source("dt_cov_fun.R")
  source("dtdt_cov_fun.R")
  source("dt2_cov_fun.R")
  source("dtdt2_cov_fun.R")
  source("dt2dt2_cov_fun.R")
  
  setwd(store_dir)
  
  return(list(num_cov_paras=num_cov_paras,
              opt_index_cov_paras=opt_index_cov_paras,
              paras_names=paras_names,
              periodic=periodic,cov_choice=cov_choice))
  
}