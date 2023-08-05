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
## v_exclude: quais os splines devem ser desconsiderados? em geral desconsidero a varib. entre sítios
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
  # data$predito <- predict(gamm, type = "link",
  #                         exclude = v_exclude,
  #                         newdata=data,
  #                         newdata.guaranteed=TRUE) |> as.vector()
  beta <- coef(gamm)
  V <- vcov.gam(gamm)
  num_beta_vecs <- nsim
  Cv <- chol(V)
  set.seed(1)
  nus <- rnorm(num_beta_vecs * length(beta))
  beta_sims <- beta + t(Cv) %*% matrix(nus, nrow = length(beta), ncol = num_beta_vecs)
  matrix_lprediction <- predict(gamm,type="lpmatrix",
                                exclude = v_exclude,
                                newdata=data,newdata.guaranteed=TRUE) 
  predict_link <- matrix_lprediction %*% beta_sims
  df_pred <- t(apply(predict_link,1,\(X) quantile(X,probs = quantiles))) %>% 
    as.data.frame()
  names(df_pred) <- paste0("Q_",quantiles)
  cbind(data,df_pred)
}
#
# f_calcPI: 
# aplica f_PredInt.GAMM para um GAMM de dois tipos: ~ contraste_z ou ~ k,p
f_calcPI <- \(gamm,
              site.posteriori ="SPigua1",
              length_pred = 150,
              n.posteriori_samples = 10000){
    # take the GAMM objects:
    f_invlink <- gamm$family$linkinv
    if(gamm$family$family=="binomial"){
      df_obs <- gamm$model |> select(-1) |> 
        mutate(y = gamm$model[,1][,1]) |> relocate(y)
      names(df_obs)[1] <- colnames(gamm$model[,1])[1]  
    }else{
      df_obs <- gamm$model
    }
    # create the new data 
    if(any(names(df_obs) |> str_detect("contraste_z"))){
      y_var <- names(df_obs)[1]
      x_var <- names(df_obs)[2]
      df_newpred <- data.frame(x = seq(min(df_obs[,x_var]),max(df_obs[,x_var]),length.out=length_pred)) |> 
        mutate(SiteCode = site.posteriori)
      names(df_newpred)[1] <- x_var  
    }else{
      y_var <- names(df_obs)[1]
      x1_var <- names(df_obs)[2]
      x2_var <- names(df_obs)[3]
      df_newpred <- expand.grid(x1 = seq(min(df_obs[,x1_var]),max(df_obs[,x1_var]),length.out=length_pred),
                                x2 = seq(min(df_obs[,x2_var]),max(df_obs[,x2_var]),length.out=length_pred)) |> 
        mutate(SiteCode = site.posteriori)
      names(df_newpred)[1:2] <- c(x1_var,x2_var)
      x_var <- "k_z"
    }
    # obtain the predictions
    ## new predictions without variability between sites
    df_newpred <- f_PredInt.GAMM(data=df_newpred,gamm=gamm,
                                 nsim=n.posteriori_samples,
                                 v_exclude=c(paste0("s(",x_var,",SiteCode)"),"s(SiteCode)"))
    ## predictions for the observed data (with variability between sites)
    ## apply the inverse link function
    ### if the case is a binomial GAMM
    if(gamm$family$family=="binomial"){
      df_newpred <- df_newpred |> 
        mutate(across(starts_with("Q_"),\(x) f_invlink(x) * sum(gamm$model[1,1])))
    }else{
      df_newpred <- df_newpred |> 
        mutate(across(starts_with("Q_"),f_invlink))
    }
    return(df_newpred)
}
# f_ggplot_PI1d:
# padroniza um scatterplot simples para f_plotCongCont
f_ggplot_PI1d <- \(df_pred,df_obs,x_var,y_var,title_contraste){
  df_obs |> 
    ggplot(aes(x=.data[[x_var]],y=.data[[y_var]])) +
    geom_point(alpha=0.2) +
    geom_line(aes(y=predito,group=SiteCode),alpha=0.2) +
    geom_line(data=df_pred,
              aes(x=.data[[x_var]],y=Q_0.5),color="darkred") +
    geom_line(data = df_pred |> 
                select(-Q_0.5) |> 
                pivot_longer(starts_with("Q_0.")),
              aes(x=.data[[x_var]],y=value,group=name),color="darkblue") +
    labs(y = "congruência", x = "z(contraste da taxa U)") +
    ggtitle(title_contraste) +
    geom_vline(xintercept = quantile(df_obs[[x_var]],probs = c(0.95,0.99)),linetype="dashed",alpha=0.3) +
    theme_bw() +
    theme(plot.title.position = "plot")
}
# f_ggplot_PI2d:
# padroniza um heatamp para f_plotCongCont
f_ggplot_PI2d <- \(df,v_range,v_breaks=NULL,title_contraste){
  f_ggplot <- \(df,facets=2){
    ggplot(df,aes(x=p,y=k,z=pred,fill=pred)) +
      geom_raster() +
      geom_contour(aes(z=pred),
                   size=0.5,color="darkblue",
                   breaks=v_breaks) +
      geom_text_contour(aes(z=pred),
                        breaks=v_breaks,
                        label.placer = label_placer_flattest()) +
      # scale_fill_gradient2(midpoint = median(v_range),
      #                      limits=v_range) +
      scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
      scale_y_reverse() +
      theme_classic() +
      # guides(fill = guide_colourbar(title = ) +
      coord_cartesian(expand = FALSE) +
      labs(fill="Prob. Cong. Pred.") +
      facet_wrap(~label,ncol = facets,scales="free",
                 labeller = as_labeller(v_labels))
  }
  v_labels <- c("Q_0.05"="0.05",
                "Q_0.25"="0.25",
                "Q_0.5"="median",
                "Q_0.75"="0.75",
                "Q_0.95"="0.95")
  l_p <- list()
  l_p[[1]] <- df |> 
    filter(label == "Q_0.5") |>
    f_ggplot() + labs(x="")
  l_p[[2]] <- df |> 
    filter(label != "Q_0.5") |>
    f_ggplot()
  ggpubr::ggarrange(plotlist = l_p,nrow=2, common.legend = TRUE, legend="bottom") |> 
    ggpubr::annotate_figure(top = ggpubr::text_grob(title_contraste,hjust=0,x=0))
}
# f_titleContraste
# padroniza os títulos dos gráficos para f_plotCongCont
f_titleContraste <- \(string){
  if(string=="cont.ideal"){
    "a) Cont. : Ideal."
  }else if(string=="cont.non_frag"){
    "b) Cont. : Sem Frag."
  }else if(string=="non_frag.ideal"){
    "c) Sem Frag. : Ideal."
  }
}
# f_plotCongCont
# função que gera os gráficos finais da congruência das SADs preditas entre contrastes
# plota a predição dos GAMM mais plausíveis segundo tabela 1 dos resultados
f_plotCongCont <- \(df){
  df_pred <- read_csv(paste0("./dados/csv/df_PIcongCont",df$pair[1],".csv"))
  if(any(names(df_pred)=="contraste_z")){
    df$predito <- predict.gam(object = l_md_congContrastes[[df$pair[1]]][["f(contraste_z)"]],
                              type = "response") * 500
    title_contraste <- f_titleContraste(df$pair[1])
    f_ggplot_PI1d(df_pred,df,"contraste_z","nCong",title_contraste = title_contraste)
  }else{
    df_pred <- df_pred |> 
      mutate(p = p_z*sd(df_ad$p) + mean(df_ad$p),
             k = k_z*sd(df_ad$k) + mean(df_ad$k)) |> 
      pivot_longer(starts_with("Q_0."),values_to = "pred",names_to = "label")
    df_pred$label <- factor(df_pred$label,levels = unique(df_pred$label)[c(3,2,4,1,5)])
    # gráficos
    v_range <- range(df_pred$pred)
    title_contraste <- f_titleContraste(df$pair[1])
    f_ggplot_PI2d(df=df_pred,title_contraste = title_contraste)
  }
}
#
#
f_diffSBtsIntPred <- \(lme){
  df_obs <- lme@frame
  #new data
  df_newdat <- expand.grid(SiteCode="SPigua1", 
                           p = seq(min(df_obs$p),max(df_obs$p), length=150),
                           land_hyp = unique(df_obs$land_hyp))
  ## Passo 2: crie as função que devem ser calculadas dos modelos a cada simulação
  ## Previstos por efeitos fixos e aleatórios
  f1 <- function(.) predict(., newdata=df_newdat)
  ## Previstos por efeitos fixos (argumento re.form=~0)
  f2 <- function(.) predict(., newdata=df_newdat, re.form=~0)
  ## Os dois bootstraps. Ajuste o argumento ncpus para o numero de cores de seu computador
  b3 <- bootMer(lme, FUN = f1, nsim=1000, parallel="multicore", ncpus=2)
  b4 <- bootMer(lme, FUN = f2, nsim=1000, parallel="multicore", ncpus=2)
  ## calcula as médias e intervalos de confiança quantílicos para cada combinação de preditoras
  ## no novo conjunto de dados
  df_newdat$mean <- apply(b3$t,2,mean)
  df_newdat$IC.low <- apply(b3$t,2,quantile, 0.025)
  df_newdat$IC.upp <- apply(b3$t,2,quantile, 0.975)
  df_newdat$mean.fixed <- apply(b4$t,2,mean)
  df_newdat$IC.low.fixed <- apply(b4$t,2,quantile, 0.025)
  df_newdat$IC.upp.fixed <- apply(b4$t,2,quantile, 0.975)
  return(df_newdat)
}
f_diffSplotPI <- \(lme,df_newpred){
  df_obs <- lme@frame
  levels(df_obs$land_hyp) <- levels(df_obs$land_hyp)[c(1,3,2)] 
  #
  df_obs %>% 
    ggplot(aes(x=p,y=diffS)) + 
    geom_point(alpha=0.1) + 
    geom_hline(yintercept = 0,color="red",alpha=0.8) +
    geom_ribbon(aes(y = mean, ymin=IC.low, ymax=IC.upp), 
                data=df_newpred, fill="grey15", alpha=0.5) +
    geom_ribbon(aes(y=mean, ymin=IC.low.fixed, ymax=IC.upp.fixed), 
                data=df_newpred, fill="white", alpha=0.5) +
    geom_line(aes(x=p, y=mean.fixed), 
              data=df_newpred,color="lightgreen") +
    labs(x="p",y="erro na estimativa da riqueza") +
    facet_wrap(~land_hyp,ncol=3,
               labeller = as_labeller(c("cont" = "Cont.",
                                        "ideal" = "Ideal.",
                                        "non_frag" = "Sem Frag.")))
}
#
# deprecated funcions:
d_f_PostPredPlotGAMM1d_obsEpred <- \(gamm,
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
d_f_PostPredPlotGAMM1d1f_obsEpred <- \(gamm,
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
d_f_l_mdGAMM_FiguraFinal_1dGAMM <- \(l_md){
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
