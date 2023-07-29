# f_TabSelGAMM
# tabela de seleção com deviance explained
# input: uma lista com os modelos candidatos nomeados
# output um data frame com o output de bbmle::AICctab + mgcv::summary$dev.exp
f_TabSelGAMM <- function(l_md){
  l_names <- names(l_md)
  df_aicctab <- AICctab(l_md,weights=TRUE,mnames = l_names) |> as.data.frame()
  df_aicctab$modelo <- row.names(df_aicctab)
  row.names(df_aicctab) <- NULL
  df_dev.exp <- ldply(l_md,.fun = \(x) summary(x)$dev.expl)
  names(df_dev.exp) <- c("modelo","dev.expl")
  df_dev.exp |> 
    inner_join(y=df_aicctab,"modelo") |> 
    arrange(dAICc) |> 
    select(modelo,dAICc:weight,dev.expl)
}
#
# f_validaGAMM
# input: um gam (md)
# uso: en um chunk indivíduo rode 'md |> f_validaGAMM()'
# outpu: dois outputs de console (k.check e summary) e dois conjuntos de gráficos diag do pacote gratia
f_validaGAMM <- \(md_name,l_md,size=0){
  if(size == 1){
    md <- l_md
    print(md_name)
  }else{
    print(paste0("Var. resposta: ",sub("efeito_","contraste ",md_name)))
    md <- l_md[[md_name]]
  }
  print("k.check:")
  print(k.check(md))
  print("summary:")
  print(summary(md))
  print("gratia::appraise:")
  print(gratia::appraise(md))
  print("gratia::draw:")
  print(gratia::draw(md))
}
 #
# f_PredInt.GAMM
## gamm: um mgcv::gam
## data: df com colunas com as preditoras do gamm
## v_exclude: quais os splines devem ser desconsiderados?
## quantiles: intervalo de predição retornado
## nsim: número de amostras da distribuição normal multivariada dos coef estimados pelo gamm
#
# output: retorna um data frame com o intevalo de predição na escala da função de ligação e as preditoras
# adaptado de https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
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
#
# f_PostPredPlotGAMM1d_obsEpred
## posterior prediction interval from a GAMM
# input: GAMM
# output: a ggplot2 scatter plot with the observed values and predictions;
# and predictions without variability between sites.
f_PostPredPlotGAMM1d_obsEpred <- \(gamm,
                                   site.posteriori ="SPigua1",
                                   length.out_x = 200,
                                   n.posteriori_samples = 10000){
  # take the GAMM objects:
  f_invlink <- gamm$family$linkinv
  v_devexp <- summary(gamm)$dev.exp
  if(gamm$family$family=="binomial"){
    df_obs <- gamm$model |> select(-1) |> 
      mutate(y = gamm$model[,1][,1]) |> relocate(y)
    names(df_obs)[1] <- colnames(gamm$model[,1])[1]  
  }else{
    df_obs <- gamm$model
  }
  # create the new data 
  y_var <- names(df_obs)[1]
  x_var <- names(df_obs)[2]
  df_newpred <- data.frame(x = seq(min(df_obs[,x_var]),
                                   max(df_obs[,x_var]),
                                   length.out=length.out_x)) |> 
    mutate(SiteCode = site.posteriori)
  names(df_newpred)[1] <- x_var
  # obtain the predictions
  ## new predictions without variability between sites
  df_pred_newdata <- f_PredInt.GAMM(data=df_newpred,
                                    gamm=gamm,
                                    nsim=n.posteriori_samples,
                                    v_exclude=c(paste0("s(",x_var,",SiteCode)"),
                                                "s(SiteCode)"))
  ## predictions for the observed data (with variability between sites)
  df_obs$predito <- predict(gamm)
  ## apply the inverse link function
  ### if the case is a binomial GAMM
  if(gamm$family$family=="binomial"){
    df_pred_newdata <- df_pred_newdata |> 
      mutate(across(predito:8,\(x) f_invlink(x) * sum(gamm$model[1,1])))
    df_obs <- df_obs |> 
      mutate(predito = f_invlink(predito) * sum(gamm$model[1,1]))  
  }else{
    df_pred_newdata <- df_pred_newdata |> 
      mutate(across(predito:8,f_invlink))
    df_obs <- df_obs |> 
      mutate(predito = f_invlink(predito))  
  }
  # create the final figure
  df_obs |> 
    ggplot(aes(x=.data[[x_var]],y=.data[[y_var]])) +
    geom_point(alpha=0.2) +
    geom_line(aes(y=predito,group=SiteCode),alpha=0.2) +
    geom_line(data=df_pred_newdata,
              aes(x=.data[[x_var]],y=Q_0.5),color="darkred") +
    geom_line(data = df_pred_newdata |> 
                select(-Q_0.5) |> 
                pivot_longer(starts_with("Q_0.")),
              aes(x=.data[[x_var]],y=value,group=name),color="darkblue") +
    labs(y = "congruência",
         subtitle = paste0("full model dev. exp. = ",round(v_devexp*100,digits = 1),"%")) +
    geom_vline(xintercept = quantile(df_obs[,x_var],probs = c(0.95,0.99)),linetype="dashed",alpha=0.3) +
    theme_bw() +
    theme(plot.caption = element_text(hjust=0))
}
#
# f_PostPredPlotGAMM1d1f_obsEpred
## função para post PI 1d 1f # pensar em formas 
f_PostPredPlotGAMM1d1f_obsEpred <- \(gamm,
                                     site.posteriori ="SPigua1",
                                     length.out_x = 200,
                                     n.posteriori_samples = 10000){
  # take the GAMM objects:
  f_invlink <- gamm$family$linkinv
  v_devexp <- summary(gamm)$dev.exp
  if(gamm$family$family=="binomial"){
    df_obs <- gamm$model |> select(-1) |> 
      mutate(y = gamm$model[,1][,1]) |> relocate(y)
    names(df_obs)[1] <- colnames(gamm$model[,1])[1]  
  }else{
    df_obs <- gamm$model
  }
  v_other.factor <- df_obs |> 
    select(-SiteCode) |> select(where(is.factor)) |> names()
  # create the new data 
  y_var <- names(df_obs)[1]
  x_var <- df_obs |> 
    select(-1) |> select(where(is.numeric)) |> names()
  if(!is_empty(v_other.factor)){
    df_newpred <- expand.grid(
      x = seq(min(df_obs[,x_var]),
              max(df_obs[,x_var]),
              length.out=length.out_x),
      factor = df_obs |> 
        select(-SiteCode) |> 
        select(where(is.factor)) |> 
        pull() |> 
        unique()
      ) |> 
      mutate(SiteCode = factor(site.posteriori))
    names(df_newpred) <- c(x_var,v_other.factor,"SiteCode")
    }else{
    df_newpred <- data.frame(x = seq(min(df_obs[,x_var]),
                                     max(df_obs[,x_var]),
                                     length.out=length.out_x)) |> 
      mutate(SiteCode = site.posteriori)
    names(df_newpred)[1] <- x_var  
  }
  # obtain the predictions
  ## new predictions without variability between sites
  df_pred_newdata <- f_PredInt.GAMM(data=df_newpred,
                                    gamm=gamm,
                                    nsim=n.posteriori_samples,
                                    v_exclude=c(paste0("s(",x_var,",SiteCode)"),
                                                "s(SiteCode)"))
  ## predictions for the observed data (with variability between sites)
  df_obs$predito <- predict(gamm)
  ## apply the inverse link function
  ### if the case is a binomial GAMM
  if(gamm$family$family=="binomial"){
    df_pred_newdata <- df_pred_newdata |> 
      mutate(across(predito:8,\(x) f_invlink(x) * sum(gamm$model[1,1])))
    df_obs <- df_obs |> 
      mutate(predito = f_invlink(predito) * sum(gamm$model[1,1]))  
  }else{
    df_pred_newdata <- df_pred_newdata |> 
      mutate(across(predito:8,f_invlink))
    df_obs <- df_obs |> 
      mutate(predito = f_invlink(predito))  
  }
  # create the final figure
  df_obs |> 
    ggplot(aes(x=.data[[x_var]],y=.data[[y_var]])) +
    geom_point(alpha=0.2) +
    geom_line(aes(y=predito,group=SiteCode),alpha=0.2) +
    geom_line(data=df_pred_newdata,
              aes(x=.data[[x_var]],y=Q_0.5),color="darkred") +
    geom_line(data = df_pred_newdata |> 
                select(-Q_0.5) |> 
                pivot_longer(starts_with("Q_0.")),
              aes(x=.data[[x_var]],y=value,group=name),color="darkblue") +
    geom_hline(yintercept = 0,color="darkred") +
    labs(y = "diiffS",
         subtitle = paste0("full model dev. exp. = ",round(v_devexp*100,digits = 1),"%")) +
    # geom_vline(xintercept = quantile(df_obs[,x_var],probs = c(0.95,0.99)),linetype="dashed",alpha=0.3) +
    theme_bw() +
    theme(plot.caption = element_text(hjust=0)) +
    facet_wrap(~.data[[v_other.factor]],nrow=1)
}
#
# f_l_mdGAMM_FiguraFinal_1dGAMM
## função para customizar o gráfico produzido por f_PostPredPlotGAMM1d_obsEpred
# adiciona título e rótulo do eixo x
f_l_mdGAMM_FiguraFinal_1dGAMM <- \(l_md){
  df_labels <- data.frame(pair = names(l_md)) |> 
    mutate(contraste = case_when(pair == "cont.ideal" ~ "contemporaneidade",
                                 pair == "cont.non_frag" ~ "fragmentação per se",
                                 TRUE ~ "área per se"),
           labels = case_when(pair == "cont.ideal" ~ "contemp. : ideal.",
                              pair == "cont.non_frag" ~ "contemp. : sem frag.",
                              TRUE ~ "sem frag. : ideal."))
  
  l_p <- llply(l_md,f_PostPredPlotGAMM1d_obsEpred)
  names(l_p) <- names(l_md)
  for(plot in names(l_p)){
    l_p[[plot]] <- l_p[[plot]] +
      labs(x=df_labels |> 
             filter(pair == plot) |> 
             pull(contraste),
           title=df_labels |> 
             filter(pair == plot) |> 
             pull(labels)) +
      theme(plot.title = element_text(hjust = 0.5,
                                      size = 12, 
                                      face = "bold"))
  }
  return(l_p)
}