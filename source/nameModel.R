f_nameModel <- function(md){
  string <- unlist(find_predictors(md))
  string <- string[string!="SiteCode"]
  paste0("f(",paste(string,collapse = ","),")")
}