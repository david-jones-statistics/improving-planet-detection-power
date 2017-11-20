opt_function <- function(opt_paras,para_type,para_opt_index,for_like){
  para_type_index <- which(names(for_like)==para_type)
  for_like[[para_type_index]][para_opt_index] <- opt_paras

  value <- -like_fun2(for_like)[[1]]
  return(value)
}