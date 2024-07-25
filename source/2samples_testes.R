# Funções para comparar a SAD observada com as SADs réplicas usando as funções do pacote twosamples
#
# o input da função é uma linha de um data frame com:
## i) path para o csv com as réplicas
## ii) path para o csv com a SAD observada
#
# o output é um df com as colunas guias (SiteCode e k) e o valor p
#
# exemplo de uso:
## df_KSrep <- ddply(df_SADrep,"SiteCode",f_ksObsRep,.parallel = TRUE)
#
# df <- df_SADrep |> filter(SiteCode == "SPigua1")
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
f_congContrastes <- \(df_pSite,n_ks=500,n_replicate=500,repo="dados/csv/congruencia_contrastes/"){
  # objetos
  l_mSADrep <- dlply(df_pSite,"land_type",\(x) read_csv(x$SADrep.path))
  df_adply <- df_pSite |> select(SiteCode,land_type) |> 
    mutate(pair = case_when(land_type == "cont" ~ "cont.non_frag",
                            land_type == "non_frag" ~ "non_frag.ideal",
                            TRUE ~ "cont.ideal")) |> 
    select(-land_type) |> 
    nest() |> 
    expand_grid(k = unique(l_mSADrep[[1]]$k)) |> 
    unnest(cols = c(data))
  # rotinas
  f_ks <- \(v1,v2){
    testeKS <- ks_test(a=v1,b=v2,nboots=n_ks)
    testeKS[2]
  }
  f_aply <- \(row){
    v_names <- unlist(strsplit(row$pair,split="[.]"))
    df_1 <- l_mSADrep[[v_names[1]]] |> filter(k==row$k)
    df_2 <- l_mSADrep[[v_names[2]]] |> filter(k==row$k)
    f_replicate <- \(){
      v_sample <- sample(1:100,2)
      v1 <- df_1[v_sample[1],] |> select(starts_with("V")) |> t() |> table() |> sort()
      v2 <- df_2[v_sample[2],] |> select(starts_with("V")) |> t() |> table() |> sort()
      cbind(row) |> 
        data.frame(index.1_2 = paste(v_sample,collapse = "_"),
                   p.valor = f_ks(v1=v1,v2=v2),row.names = NULL)
    }
    rdply(n_replicate,f_replicate(),.id=NULL)
  }
  df_toWrite <- adply(df_adply,1,f_aply)
  write_csv(df_toWrite,file = paste0(repo,df_pSite$SiteCode[1],".csv"))
}
#
# f_summarise_SAD_MNEE
f_summarise_SAD_MNEE <- \(df){
#@ df: df por site, k, e  land_type
#@ e.g. ddply(.,c("SiteCode","k","land_type"))
 cbind(df[1,c("SiteCode","k","land_type")],with(df,data.frame(
   nCongKS = sum(p.KS>0.05),
   Smed = mean(S),
   Ssd = sd(S),
   Smin = min(S),
   Smax = max(S)))
  )
}
# f_logOR_land_type
f_logOR_land_type <- \(df,
                       # pairs=list(c("contemp","non_frag"),
                       #            c("contemp","ideal"))){
  pairs=list(c("contemp","non_frag"),
             c("contemp","ideal"),
             c("non_frag","ideal"))){
  df <- df %>% 
    mutate(propCong = case_when(nCongKS == 100 ~ 99.9 / 100,
                                nCongKS == 0 ~ 0.1 / 100,
                                TRUE ~ nCongKS / 100))
  f <- \(v_string){
    p1 <- df |> filter(land_type == v_string[1]) |> pull(propCong)
    p1 <- p1 / (1-p1)
    p2 <- df |> filter(land_type == v_string[2]) |> pull(propCong)
    p2 <- p2 / (1-p2)
    data.frame(logOR_pair = paste(v_string,collapse = "."),
               logOR_value = log( p1 / p2 ) )
  }
  cbind(df,ldply(pairs,f))
}
f_logOR_dfi <- \(x,y){
  x <- x/(1-x)
  y <- y/(1-y)
  log(x/y)
}

f_contrasteSAD <- \(dff,pares){
  dff <- dff %>% 
    pivot_wider(names_from=land_type,
                values_from=logitoKS)
  formals(f_logOR)$dfi <- dff
  f_par.i <- \(v_par){
    land_contrafac <- str_split_1(v_par,"-")
    names(land_contrafac) <- c("numerador","denominador")
    df_return <- dff %>% select(SiteCode,k) %>% distinct()
    df_return <- cbind(
      df_return,
      f_logOR(cpair=land_contrafac)  
    ) %>% mutate(contraste = v_par) %>% 
      relocate(contraste) %>% as.data.frame()
    return(df_return)
  }
  df_return <- lapply(pares,f_par.i) %>% rbind.fill()
  df_return %>% select(contraste:denominator)
}
f_logOR <- \(dfi,cpair){
  v_itself <- 
    dfi[dfi$taxaU==cpair["numerador"],
        cpair["numerador"]] - dfi[dfi$taxaU==cpair["denominador"],
                                  cpair["denominador"]]
  v_pristine <- 
    dfi[dfi$taxaU=="ideal",
        cpair["numerador"]] - dfi[dfi$taxaU=="ideal",
                                  cpair["denominador"]]
  v_num <- 
    dfi[dfi$taxaU==cpair["numerador"],
        cpair["numerador"]] - dfi[dfi$taxaU==cpair["numerador"],
                                  cpair["denominador"]]
  if(unname(cpair["denominador"])!="ideal"){
    v_deno <- 
      dfi[dfi$taxaU==cpair["denominador"],
          cpair["numerador"]] - dfi[dfi$taxaU==cpair["denominador"],
                                    cpair["denominador"]]
  }
  # return
  data.frame(
    itself = unname(v_itself),
    pristine = unname(v_pristine),
    numerator = unname(v_num),
    denominator = ifelse(exists("v_deno"),unname(v_deno),NA) %>% unlist()
  )
}


f_contraste_Umed <- \(dfUrep,
                      path_U="dados/csv/taxaU/df_U.csv",
                      path_land_effect="dados/csv/taxaU/df_contrastes.csv",
                      nobj_export=c("df_U","df_contrastes")){
  df_U <- dfUrep |> select(SiteCode:k)
  df_U$Umed <- apply(select(df_Urep,-c(SiteCode:d)),1,mean)
  df_U$Usd <- apply(select(df_Urep,-c(SiteCode:d)),1,sd)
  df_contrastes <- df_U |>  
    inner_join(df_p,by="SiteCode") |>
    pivot_wider(names_from = land_type, values_from = c(Umed,Usd)) |> 
    mutate(area_dif = Umed_non_frag - Umed_ideal,
           frag.perse_dif = Umed_contemp - Umed_non_frag,
           frag.total_dif = Umed_contemp - Umed_ideal,
           area_ratio = Umed_non_frag / Umed_ideal,
           frag.perse_ratio = Umed_contemp / Umed_non_frag,
           frag.total_ratio = Umed_contemp / Umed_ideal,
           area_logratio = log(area_ratio),
           frag.perse_logratio = log(frag.perse_ratio),
           frag.total_logratio = log(frag.total_ratio),
           area_logOR = f_logOR_dfi(Umed_non_frag, Umed_ideal),
           frag.perse_logOR = f_logOR_dfi(Umed_contemp, Umed_non_frag),
           frag.total_logOR = f_logOR_dfi(Umed_contemp, Umed_ideal)) |> 
    select(-starts_with("U"))
  write_csv(df_contrastes,file=path_land_effect)
  write_csv(df_U,file=path_U)
  if(!is.null(nobj_export)){
    sapply(nobj_export,\(i) assign(i,get(i),envir = .GlobalEnv))
  }
}
f_MoranTest_GAMM <- \(md){
  library(spdep)
  dfmd <- md$model
  dfmd$residuals <- residuals(md)
  if(sum(names(dfmd)%in%c("lat","long"))<1){
    df_coord <- read_csv(file = "dados/df_dados_disponiveis.csv") %>% 
      mutate(lat = ifelse(is.na(lat_correct),lat,lat_correct),
             long = ifelse(is.na(long_correct),long,long_correct),
             Sitecode = factor(SiteCode)) %>% 
      select(SiteCode,lat,long)
    dfmd <- inner_join(
      dfmd,
      df_coord
    )
  }
  dfmd_avgbySite <- dfmd %>% 
    group_by(SiteCode) %>% 
    summarise(mean_res = mean(residuals),
              lat = first(lat),
              long = first(long)) %>% 
    ungroup()
  # Prepare spatial data
  coordinates <- dfmd_avgbySite[, c("lat","long")]
  coordinates <- as.matrix(coordinates)
  nb <- knn2nb(knearneigh(coordinates, k=4))  
  listw <- nb2listw(nb)
  # Calculate Moran's I for residuals
  moran_output <- moran.test(dfmd_avgbySite$mean_res, listw)
  # return
  data.frame(
    Statistic = c(
      "Moran I statistic (res)", 
      "Expectation", 
      "Variance", 
      "Standard Deviate", 
      "p-value"
    ),
    Value = c(
      moran_output$estimate[["Moran I statistic"]], 
      moran_output$estimate[["Expectation"]], 
      moran_output$estimate[["Variance"]], 
      moran_output$statistic, 
      moran_output$p.value
    )
  ) %>% pivot_wider(names_from=Statistic,values_from = Value) %>% 
    as.data.frame()
}


# 
# 
# f_congContrastes <- \(df_pSite){
#   df_mSAD <- dlply(df_pSite,"land_type",\(x) f_sampRowsbyk(x$SADrep.path)) |> 
#     bind_rows(.id="land_type") |> 
#     relocate(land_type)
#   f_outputKS <- \(df,nboots=1000){
#     ## cont e os outros dois tratamentos:
#     # dados
#     df_cont <- df |> filter(land_type == "cont") |> data.frame()  |> select(-c(1:3)) |> as.matrix()
#     df_comp <- df |> filter(land_type != "cont") |> data.frame()  |> select(-c(1:3)) |> as.matrix()
#     # 2o nível: função aplicada em cada linha de df_comp
#     f_KS <- \(v_ref){
#       f_testeKS <- \(row_df){
#         v_row <- sort(table(row_df))
#         testeKS <- ks_test(a=v_ref,b=v_row,nboots=nboots)
#         testeKS[2]
#       }
#       aaply(df_comp,1,f_testeKS)
#     }
#     # 1o nível: função aplicada em cada linha de df_cont
#     f_appKS <- \(X){
#       v_X <- sort(table(X))
#       f_KS(v_ref=v_X)
#     }
#     df_return <- adply(df_cont,1,f_appKS) |> select(-1) |> 
#       mutate(land_type = "cont") |> relocate(land_type)
#     names(df_return)[-1] <- df |> filter(land_type != "cont") |> pull(land_type)
#     df_return1 <- df_return |> 
#       pivot_longer(-land_type,names_to = "pair",values_to = "p_value") |> 
#       mutate(pair=paste0(land_type,".",pair),
#              SiteCode = unique(df$SiteCode),
#              k = unique(df$k)) |> select(-land_type) |> 
#       relocate(p_value,.after = last_col())
#     # 
#     ## non_frag e ideal:
#     # dados
#     df_non_frag <- df |> filter(land_type == "non_frag") |> data.frame()  |> select(-c(1:3)) |> as.matrix()
#     df_comp <- df |> filter(land_type == "ideal") |> data.frame()  |> select(-c(1:3)) |> as.matrix()
#     # 2o nível: função aplicada em cada linha de df_comp
#     f_KS <- \(v_ref){
#       f_testeKS <- \(row_df){
#         v_row <- sort(table(row_df))
#         testeKS <- ks_test(a=v_ref,b=v_row,nboots=nboots)
#         testeKS[2]
#       }
#       aaply(df_comp,1,f_testeKS)
#     }
#     # 1o nível: função aplicada em cada linha de df_cont
#     f_appKS <- \(X){
#       v_X <- sort(table(X))
#       f_KS(v_ref=v_X)
#     }
#     df_return <- adply(df_non_frag,1,f_appKS) |> select(-1) |> 
#       mutate(land_type = "non_frag") |> relocate(land_type)
#     names(df_return)[-1] <- "ideal"
#     df_return2 <- df_return |> 
#       pivot_longer(-land_type,names_to = "pair",values_to = "p_value") |> 
#       mutate(pair=paste0(land_type,".",pair),
#              SiteCode = unique(df$SiteCode),
#              k = unique(df$k)) |> select(-land_type) |> 
#       relocate(p_value,.after = last_col())
#     # return
#     rbind.fill(df_return1,df_return2)
#   }
#   ddply(df_mSAD,"k",f_outputKS)
# }
##
#
# uma possível implementação do teste de permutação baseado na distância multivariada 
# autoria: chatgpt
# library(vegan)
# 
# # Example data
# vector1 <- c(10, 15, 5, 20)  # Species abundances in sample 1
# vector2 <- c(8, 12, 10, 18)  # Species abundances in sample 2
# 
# # Calculate the observed test statistic
# obs_stat <- vegdist(rbind(vector1, vector2), method = "euclidean")  # Bray-Curtis dissimilarity
# 
# # Set the number of permutations
# n_permutations <- 1000
# 
# # Create an empty vector to store permuted test statistics
# perm_stats <- numeric(n_permutations)
# 
# # Perform the permutation test
# for (i in 1:n_permutations) {
#   # Pool the data and randomly permute the labels
#   perm_data <- sample(c(vector1, vector2))
#   
#   # Calculate the test statistic for the permuted data
#   perm_stat <- vegdist(matrix(perm_data, ncol = length(vector1)), method = "euclidean")
#   
#   # Store the permuted test statistic
#   perm_stats[i] <- perm_stat
# }
# 
# # Calculate the p-value
# p_value <- sum(perm_stats >= obs_stat) / n_permutations
# 
# # Print the results
# cat("Observed test statistic:", obs_stat, "\n")
# cat("p-value:", p_value, "\n")
