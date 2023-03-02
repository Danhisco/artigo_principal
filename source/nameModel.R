f_nameModel <- function(md){
  paste(unlist(find_predictors(md)),collapse = " * ")
}