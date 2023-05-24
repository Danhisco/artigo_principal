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
  
f_resultsMN <- function(df){
  # dados
  SADobs <- read_csv(df$SADobs.path) |> arrange(N) |> pull(N)
  m_SADrep <- read.csv(df$SADrep.path)
  df_return <- with(df,cbind(m_SADrep[,1:2],land_type))
  # subfunções
  f_KS <- function(SADrep, v_obs = SADobs, boots=1000){
    testeKS <- ks_test(a=v_obs,b=SADrep,nboots = boots)
    testeKS[2]
  }
  f_repKS <- function(X){
    v_rep <- sort(table(X))
    f_KS(SADrep = v_rep)
  }
  f_S <- function(X) length(table(X))
  # rotinas
  df_return$p.KS <- apply(as.matrix(m_SADrep[,-(1:2)]),1,f_repKS)
  df_return$S <- apply(as.matrix(m_SADrep[,-(1:2)]),1,f_S)
  # return
  return(df_return)
}
#
## congruência contrastes
# exemplo de uso:
# df_congContrastes <- ddply(df_SADrep,"SiteCode",f_congContrastes,.parallel = TRUE)
# função para se aplicar por sítio
f_congContrastes <- \(df_pSite){
  df_mSAD <- dlply(df_pSite,"land_type",\(x) f_sampRowsbyk(x$SADrep.path)) |> 
    bind_rows(.id="land_type") |> 
    relocate(land_type)
  f_outputKS <- \(df,nboots=1000){
    ## cont e os outros dois tratamentos:
    # dados
    df_cont <- df |> filter(land_type == "cont") |> data.frame()  |> select(-c(1:3)) |> as.matrix()
    df_comp <- df |> filter(land_type != "cont") |> data.frame()  |> select(-c(1:3)) |> as.matrix()
    # 2o nível: função aplicada em cada linha de df_comp
    f_KS <- \(v_ref){
      f_testeKS <- \(row_df){
        v_row <- sort(table(row_df))
        testeKS <- ks_test(a=v_ref,b=v_row,nboots=nboots)
        testeKS[2]
      }
      aaply(df_comp,1,f_testeKS)
    }
    # 1o nível: função aplicada em cada linha de df_cont
    f_appKS <- \(X){
      v_X <- sort(table(X))
      f_KS(v_ref=v_X)
    }
    df_return <- adply(df_cont,1,f_appKS) |> select(-1) |> 
      mutate(land_type = "cont") |> relocate(land_type)
    names(df_return)[-1] <- df |> filter(land_type != "cont") |> pull(land_type)
    df_return1 <- df_return |> 
      pivot_longer(-land_type,names_to = "pair",values_to = "p_value") |> 
      mutate(pair=paste0(land_type,".",pair),
             SiteCode = unique(df$SiteCode),
             k = unique(df$k)) |> select(-land_type) |> 
      relocate(p_value,.after = last_col())
    # 
    ## non_frag e ideal:
    # dados
    df_non_frag <- df |> filter(land_type == "non_frag") |> data.frame()  |> select(-c(1:3)) |> as.matrix()
    df_comp <- df |> filter(land_type == "ideal") |> data.frame()  |> select(-c(1:3)) |> as.matrix()
    # 2o nível: função aplicada em cada linha de df_comp
    f_KS <- \(v_ref){
      f_testeKS <- \(row_df){
        v_row <- sort(table(row_df))
        testeKS <- ks_test(a=v_ref,b=v_row,nboots=nboots)
        testeKS[2]
      }
      aaply(df_comp,1,f_testeKS)
    }
    # 1o nível: função aplicada em cada linha de df_cont
    f_appKS <- \(X){
      v_X <- sort(table(X))
      f_KS(v_ref=v_X)
    }
    df_return <- adply(df_non_frag,1,f_appKS) |> select(-1) |> 
      mutate(land_type = "non_frag") |> relocate(land_type)
    names(df_return)[-1] <- "ideal"
    df_return2 <- df_return |> 
      pivot_longer(-land_type,names_to = "pair",values_to = "p_value") |> 
      mutate(pair=paste0(land_type,".",pair),
             SiteCode = unique(df$SiteCode),
             k = unique(df$k)) |> select(-land_type) |> 
      relocate(p_value,.after = last_col())
    # return
    rbind.fill(df_return1,df_return2)
  }
  ddply(df_mSAD,"k",f_outputKS)
}
