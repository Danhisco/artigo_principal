## Função que calcula a estatística da somatória dos chi2
f_chi.statistic <- \(df_1,df_2){
  names(df_2) <- names(df_1)
  df_samples <- full_join(df_1,df_2,by=names(df_1)[1]) |> 
    mutate(across(Freq.x:Freq.y,\(x) ifelse(is.na(x),0,x))) |> 
    pivot_longer(2:3,names_to = "sample",values_to = "Freq") |> 
    rename(Abund=1) |> 
    mutate(sample=factor(sample,labels=c("x","y")),
           Abund=as.numeric(as.character(Abund)))
  # calculo do chi2.obs, df = df_samples[3:4,]
  f_chi2 <- \(df)  ((df$Freq[1] - df$Freq[2])^2)/(df$Freq[1] + df$Freq[2])
  df_chi.obs <- ddply(df_samples,"Abund",f_chi2)
  chi.value = sum(df_chi.obs$V1)
  return(chi.value)
} # não usado
## Função que calculo o p valor do chi quadrado para duas SADs amostras 
f_chi2Boot0 <- function(df_sample1,df_sample2,n_sim=1e4){
  df_1 <- df_sample1$N |> table() |> as.data.frame()
  df_2 <- df_sample2$N |> table() |> as.data.frame()
  # ci2 observado
  v_chi.obs <- f_chi.statistic(df_1,df_2)
  # calculo do chi2.rep
  ## df da onde serão sorteados as novas SADs amostras
  df_samples <- rbind(df_sample1,df_sample2) |> 
    mutate(Prop = N/sum(N))
  ## função para obter a SAD amostra
  f_samples <- \(){
    sample(df_samples[,1],
           size=sum(df_samples[,2])/2,replace = T,
           prob = df_samples[,3]) |> 
      table() |> table() |> as.data.frame()
  }
  ## função para criar a distribuição de chi2 sob hipótese nula
  f_chi.rep <- \(){
    df_1 <- f_samples()
    df_2 <- f_samples()
    f_chi.statistic(df_1,df_2)
  } 
  ## rotina
  v_chi.rep <- replicate(n_sim,f_chi.rep(),T)
  ## p valor: proporção de chi2 sob hipótese nula que são maiores ou iguais ao chi2 observado
  p_value <- length(v_chi.rep[v_chi.rep>=v_chi.obs]) / length(v_chi.rep)
  return(p_value)
} # não usado
