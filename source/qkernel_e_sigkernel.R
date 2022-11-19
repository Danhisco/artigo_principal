library(rmutil)
library(lamW)
qkernel<- function(sigma, kernel, p, density=20852/50, npoints){
  kernel <- match.arg(kernel, choices=c("normal","gaussian","laplace","uniform"))
  d_ind_MA  <- 100/sqrt(density)
  if(kernel=="laplace"){
    b_laplace <- sigma / sqrt(2)
    X_laplace <- d_ind_MA * round(rlaplace(npoints, s=b_laplace) / d_ind_MA)
    Y_laplace <- d_ind_MA * round(rlaplace(npoints, s=b_laplace) / d_ind_MA)
    dist_laplace <- sqrt(X_laplace^2+Y_laplace^2)
    result <- quantile(dist_laplace, p)
  }
  if(kernel=="normal"|kernel=="gaussian"){
    b_norm <- sigma
    X_norm <- d_ind_MA * round(rnorm(npoints, sd=b_norm) / d_ind_MA)
    Y_norm <- d_ind_MA * round(rnorm(npoints, sd=b_norm) / d_ind_MA)
    dist_norm <- sqrt(X_norm^2+Y_norm^2)
    result <- quantile(dist_norm, p)
  }
  if(kernel=="uniform"){
    b_unif <- sigma/2 # sigma correspondente ao range da distribuição uniforme
    X_unif <- d_ind_MA * round(runif(npoints, min = -b_unif, max = b_unif) / d_ind_MA)
    Y_unif <- d_ind_MA * round(runif(npoints, min = -b_unif, max = b_unif) / d_ind_MA)
    dist_unif <- sqrt(X_unif^2+Y_unif^2)
    result <- quantile(dist_unif, p)
  }
  return(unname(result))
}
sigkernel <- function(kernel = "laplace", p, density,
                      npoints =1e5, sigma.min = 1e-6, sigma.max= 1e6){
  distance <- 100/sqrt(density)
  f1 <- function(x) distance - qkernel(x, kernel, p, density, npoints)
  v_ <- uniroot( f1 , lower = sigma.min, upper = sigma.max)
  return(v_$root)
}