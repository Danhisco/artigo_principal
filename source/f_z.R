f_z <- function(x){
  m <- base::mean(x,na.rm=FALSE)
  sd <- sd(x,na.rm=FALSE)
  output <- (x-m)/sd
  return(output)
}