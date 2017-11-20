load_main_functions <- function(){
  
  source("fit_activity_gp_model2.R")
  
  store_dir <- getwd()
  
  setwd("./functions/")
  source("build_full_cov_mat.R")
  source("cov_make.R")
  source("general_cov.R")
  source("general_gp_log_lik_packed.R")
  source("general_gp_plot_packed.R")
  source("general_var.R")
  source("get_data.R")
  source("get_preset_cov_funs.R")
  source("main_optimize_loop.R")
  source("mean_model.R")
  source("opt_batch.R")
  source("opt_batch_reuse_mats.R")
  source("opt_function.R")
  source("opt_function_reuse_mats.R")
  source("opt_period_batch.R")
  source("optimize_parameters2.R")
  source("planet_model.R")
  source("planet_opt_function_reuse_mats2.R")
  source("set_priors.R")
  
  setwd(store_dir)
  
}