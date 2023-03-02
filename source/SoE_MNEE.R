# f_SoE_MNEE calcula o p para aquela extensão espacial 
f_SoE_MNEE <- function(df,n_escalas=6,n_replicas=4){
  # estima U para diferentes extensões espaciais da paisagem local, pressupõe um único cenário de dispersão per capita
  # 
  # colunas de df:
  # S_obs: riqueza observada
  # disp_range: desvio padrão da função de dispersão
  # txt.file: arquivo com o mapa de cobertura vegetal usado em MNEE
  # lado_km: lado da extensão espacial informado no txt.file, a maior extensão usada
  # 
  f_simCoalescente <- function(land_name){
    v <- dinamica_coalescente(U = 1.25e-06, S=df$S_obs, disp_range = df$d, N_simul = 1,
                              landscape = land_name)
    v$U_est
  }
  m_full <- read.table(df$txt.file) |> as.matrix()
  l_df_U <- list()
  for(i in 1:n_escalas){
    index <- (nrow(m_full)-nrow(m_full)/(2^(i-1)) )/2
    m_sim <- m_full[(floor(index)+1):(nrow(m_full)-ceiling(index)),
                    (floor(index)+1):(ncol(m_full)-ceiling(index))]
    p <- length(m_sim[m_sim!=0])/length(m_sim)
    v_U <- replicate(n=n_replicas,f_simCoalescente(land_name = m_sim))
    l_df_U[[i]] <- cbind(data.frame(SiteCode = df$SiteCode,
                                    lado_km=16.02/(2^(i-1)),
                                    k=df$k,
                                    p=p),
                         matrix(v_U,nrow=1))
  }
  df_write <- rbind.fill(l_df_U)
  write_csv(df_write,file = paste0("../csv/SoE/",df$SiteCode,"_k",round(df$k*100,0),".csv"))
}
# a versão a seguir não calcula a proporção de cobertura vegetal da paisagem (p), pois é sempre 1
f_SoE_MNEE_null <- function(df,n_escalas=6,n_replicas=10){
  f_simCoalescente <- function(land_name){
    v <- dinamica_coalescente(U = 1.25e-06, S=df$S_obs, disp_range = df$d, N_simul = 1,
                              landscape = land_name)
    v$U_est
  }
  m_full <- read.table(df$txt.file) |> as.matrix()
  l_df_U <- list()
  for(i in 1:n_escalas){
    index <- ( nrow(m_full)-(nrow(m_full)/ (2^(i-1))) ) /2
    m_sim <- m_full[(floor(index)+1):(nrow(m_full)-ceiling(index)),
                    (floor(index)+1):(ncol(m_full)-ceiling(index))]
    v_U <- replicate(n=n_replicas,f_simCoalescente(land_name = m_sim))
    l_df_U[[i]] <- cbind(data.frame(SiteCode = df$SiteCode,
                                    lado_km=16.02/(2^(i-1)),
                                    k=df$k),
                         matrix(v_U,nrow=1))
  }
  df_write <- rbind.fill(l_df_U)
  write_csv(df_write,file = paste0("../csv/SoE/",df$SiteCode,"_k",round(df$k*100,0),".csv"))
}
#
#
#
# funções para criar as variáveis de interesse
f_lag.U <- function(df){
  v_Udiff <- diff(range(df$Umed))
  df$Ulag <- c((df$Umed[1:(nrow(df)-1)] - df$Umed[2:nrow(df)])/v_Udiff,NA)
  df <- df |> mutate(Ulag0 = ifelse(Ulag<0,0,Ulag))
  df$Ucsum <- cumsum(df$Ulag)
  df$Ucsum0 <- cumsum(df$Ulag0)
  df <- df |> mutate(Ucsum0 = ifelse(Ucsum0>1,1,Ucsum0))
  return(df)
}
f_k.SoE <- \(df,limiar=0.99){
  if(max(df$Ucsum,na.rm = T)<limiar){
    df |> head(n=1) |> mutate(SoE.99=NA)
  }else{
    df |> filter(Ucsum>=limiar) |> head(n=1) |> mutate(SoE.99=as.numeric(as.character(lado_factor)))
  }
}
#