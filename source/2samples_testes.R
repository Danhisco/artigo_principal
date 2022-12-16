# Funções para comparar a SAD observada com as SADs réplicas usando as funções do pacote twosamples
#
# o input da função é uma linha de um data frame com:
## i) path para o csv com as réplicas
## ii) path para o csv com a SAD observada
#
# o output é um df com as colunas guias (SiteCode e k) e o valor p
#
# exemplo de uso:
## df_KSrep <- adply(df_SADrep,1,f_ksObsRep,.parallel = TRUE)
#
f_ksObsRep <- function(df){
  # dados
  SADobs <- read_csv(df$SADobs.path) |> arrange(N) |> pull(N)
  m_SADrep <- read.csv(df$SADrep.path)
  df_return <- m_SADrep[,1:2] 
  # subfunções
  f_KS <- function(SADrep, v_obs = SADobs, boots=1000){
    testeKS <- ks_test(a=v_obs,b=SADrep,nboots = boots)
    testeKS[2]
  }
  f_rep <- function(X){
    v_rep <- sort(table(X))
    f_KS(SADrep = v_rep)
  }
  # rotina
  df_return$p.value <- apply(as.matrix(m_SADrep[,-(1:2)]),1,f_rep)
  return(df_return)
}
f_DTSrep <- function(df){
  # dados
  SADobs <- read_csv(df$SADobs.path) |> arrange(N) |> pull(N)
  m_SADrep <- read.csv(df$SADrep.path)
  df_return <- m_SADrep[,1:2] 
  # subfunções
  f_DTS <- function(SADrep, v_obs = SADobs, boots=1000){
    testeDTS <- dts_test(a=v_obs,b=SADrep,nboots=boots,keep.boots=FALSE)
    testeDTS[2]
  }
  f_rep <- function(X){
    v_rep <- sort(table(X))
    f_DTS(SADrep = v_rep)
  }
  # rotina
  df_return$p.value <- apply(as.matrix(m_SADrep[,-(1:2)]),1,f_rep)
  return(df_return)
}
