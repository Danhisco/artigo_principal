## gamm: um mgcv::gam
## data: df com colunas com as preditoras do gamm
## v_exclude: quais os splines devem ser desconsiderados?
## quantiles: intervalo de predição retornado
## nsim: número de amostras da distribuição normal multivariada dos coef estimados pelo gamm
#
# output: retorna um data frame com o intevalo de predição na escala da função de ligação e as preditoras
f_PredInt.GAMM <- \(gamm,
                    data, 
                    v_exclude = c("s(k_z,SiteCode)","s(SiteCode)"),
                    quantiles = c(0.05,0.25,0.5,0.75,0.95),
                    nsim=10000){
  data$predito <- predict(gamm, type = "link",
                          exclude = v_exclude,
                          newdata=data,
                          newdata.guaranteed=TRUE) |> as.vector()
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
