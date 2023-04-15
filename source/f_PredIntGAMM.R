f_PredInt.GAMM <- \(data, gamm, 
                    quantiles = c(0.05,0.25,0.5,0.75,0.95),nsim=10000){
  data$predito <- predict(gamm, type = "link",
                          exclude = c("s(k_cont_z,SiteCode)","s(SiteCode)"),
                          newdata=data,newdata.guaranteed=TRUE) |> as.vector()
  beta <- coef(gamm)
  V <- vcov.gam(gamm)
  num_beta_vecs <- nsim
  Cv <- chol(V)
  set.seed(1)
  nus <- rnorm(num_beta_vecs * length(beta))
  beta_sims <- beta + t(Cv) %*% matrix(nus, nrow = length(beta), ncol = num_beta_vecs)
  matrix_lprediction <- predict(gamm,type="lpmatrix",
                                exclude = c("s(k_cont_z,SiteCode)","s(SiteCode)"),
                                newdata=data,newdata.guaranteed=TRUE) 
  predict_link <- matrix_lprediction %*% beta_sims
  df_pred <- t(apply(predict_link,1,\(X) quantile(X,probs = quantiles))) %>% 
    as.data.frame()
  names(df_pred) <- paste0("Q_",quantiles)
  cbind(data,df_pred)
}