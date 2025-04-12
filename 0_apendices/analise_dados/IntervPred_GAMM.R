# pacotes
library(gratia)
library(doMC)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
library(readr)
library(stringr)
library(tidyr)
library(bbmle)
library(DHARMa)
# library(lme4)
library(mgcv)
library(plyr)
library(dplyr)
## funções de ajuste e de plot
# source("source/2samples_testes.R")
# source("source/general_tools.R")
# source("source/GAMMtools.R")
# source("source/fig_tools.R")

# l_md <- readRDS(paste0(v_path,"rds/l_md_simples_apudPedersen2019_tp.rds"))
# df_tabsel <- read_csv(paste0(v_path,"rds/tabsel_simples_tp_e_cr.csv")) %>%
#   filter(dAICc==0,grepl("tp::",modelo)) %>% 
#   mutate(modelo = gsub("tp::","",modelo)) %>% 
#   rename(contraste=pair)
# l_md <- dlply(df_tabsel,"contraste",\(dfi){
#   with(dfi,{l_md[[contraste]][[modelo]]})
# })
#################################################################
# função para criar o new data fixo
f_dfmd2 <- \(dff,length_pred = 100,site.posteriori ="SPigua1"){
  #1o criar predição para todo o intervalo de dados
  f_bypert <- \(dfi){
    v_range <- range(dfi$Uefeito)
    v_k <- range(dfi$k_z)
    v_forest <- levels(dfi$forest_succession)
    expand.grid(Uefeito = seq(v_range[1],v_range[2],length.out=length_pred),
                k_z = seq(v_k[1],v_k[2],length.out=length_pred),
                forest_succession = v_forest,
                SiteCode=site.posteriori)
    }
  df_pred <- ddply(dff,"forest_succession",f_bypert)
  #2o cortar a parte extra
  ## dados para ajustar os modelos de corte
  f_dfbuffer <- \(dfi){
    v_range <- range(dfi$Uefeito)
    v_quant <- quantile(dfi$Uefeito,probs=c(0.01,0.99))
    filter(dfi,Uefeito%in%v_range) %>% 
      select(-c(forest_succession,SiteCode)) %>% 
      mutate(ext_class=ifelse(Uefeito==v_range[1],"min","max"),
             Uefeito_quant=ifelse(ext_class=="min",v_quant[1],v_quant[2])) %>% 
      relocate(ext_class,.after=last_col()) %>% 
      relocate(Uefeito_quant,.after=Uefeito)
  }
  df_buf <- ddply(dff,c("k_z","forest_succession"),f_dfbuffer)
  # df_buf %>%
  #   df_newpred %>% 
  #   pivot_longer(starts_with("Uefeito")) %>%
  #   filter(name=="Uefeito") %>%
  #   ggplot(aes(x=k_z,y=value)) + #,color=ext_class
  #   geom_point() +
  #   geom_line() +
  #   geom_smooth(method="gam",se=FALSE) +
  #   # scale_color_manual(values=c("darkred","darkgreen")) +
  #   facet_wrap(~forest_succession,ncol=1)
  ## ajuste e predição
  f_gam <- \(dfi){
    gam(Uefeito~s(k_z,bs="cr"),data=dfi)
  }
  l_md <- dlply(df_buf,
                c("ext_class","forest_succession"),
                f_gam)
  ## fazer a predição média do modelo e guardar
  l_df_ref <- lapply(levels(dff$forest_succession),\(li){
    lmd <- l_md[paste0(c("max.","min."),li)]
    #names(lmd) <- gsub(paste0(li,"."),"",names(lmd)) 
    df_ref <- lapply(names(lmd),\(i){
      md <- lmd[[i]]
      dfr <- filter(df_pred,forest_succession==li)
      dfr[[ gsub(paste0("\\.",li),"",i) ]] <- predict.gam(md,dfr)
      return(dfr)
    }) %>% Reduce("inner_join",.)
  })
  names(l_df_ref) <- levels(dff$forest_succession)
  ## filtrar para cada k os valores entre predições
  library(data.table)
  df_newpred <- lapply(names(l_df_ref),\(li){
    df_pred <- l_df_ref[[li]]
    setDT(df_pred)
    as.data.frame(df_pred[, 
                          .SD[.SD$min <= .SD$Uefeito & .SD$Uefeito <= .SD$max,
                              .(Uefeito,forest_succession,SiteCode)], 
                          by = k_z])
  }) %>% do.call("rbind",.)
  ## 
  return(df_newpred)
}
####################### versão predição mais simples
# f_pred_simples <- \(vpath_md){
#   l_md <- readRDS(vpath_md)
#   md <- l_md[["te(logU/U,k)"]];rm(l_md);gc() 
#   dfmd <- f_dfmd2(md$model) %>% 
#     mutate(lat = filter(md$model,SiteCode=="SPigua1") %>% pull(lat) %>% unique,
#            long = filter(md$model,SiteCode=="SPigua1") %>% pull(long) %>% unique)
#   dfmd$pred <- predict.gam(md,newdata = dfmd,exclude = "s(lat,long)")
#   rm(md);gc()
#   return(dfmd)
# }
# l_df_pred_simples <- lapply(l_path$te,f_pred_simples)
# names(l_df_pred_simples) <- str_extract(l_path$te,"(?<=md_)(.+?)(?=\\.rds)")
# saveRDS(l_df_pred_simples,file="./dados/csv_SoE/rds/l_df_predict_gam_simples.rds")











# f_dfmd <- \(dff,byforest,length_pred = 150,site.posteriori ="SPigua1"){
#   v_range <- range(dff$Uefeito)
#   df_newpred <- select(dff,-logOR,-Uefeito) %>% 
#     mutate(SiteCode=site.posteriori) %>% head(n=1)
#   try({
#     df_newpred <- cbind(
#       data.frame(Uefeito = seq(v_range[1],v_range[2],length.out=length_pred)),
#       df_newpred)
#   },silent = TRUE)
#   if(byforest){
#     df_return <- adply(as.character(unique(dff$forest_succession)),1,\(vstr){
#       vrange <- filter(dff,forest_succession==vstr) %>% pull(Uefeito) %>% range
#       filter(df_newpred,Uefeito>=vrange[1] & Uefeito<=vrange[2])
#     }) %>% select(-X1)
#     return(df_return)
#   }else{
#     return(df_newpred)
#   }
# }
# função que realiza a predição a posteriori dado o modelo, new data e o vetor de exclusão de componentes
f_predictions <- \(gamm,nsim,to_exclude,df_newpred,quants=c(0.05,0.5,0.95)){
  coef_samples <- MASS::mvrnorm(n=nsim, mu=coef(gamm), Sigma=vcov(gamm))
  matrix_lprediction <- predict(gamm,
                                type ="lpmatrix",
                                exclude = to_exclude,
                                newdata = df_newpred,
                                newdata.guaranteed = TRUE)
  df_pred <- matrix_lprediction %*% t(coef_samples)
  df_pred <- t(apply(df_pred,1,\(X) quantile(X,probs = quants))) %>% 
    as.data.frame()
  names(df_pred) <- paste0("Q_",quants)
  return(df_pred)
}
# criação do new data para predição fixo e aleatório
# f_dfmd_aleat <- \(dfmd,df_newpred){
#   ddply(dfmd,"SiteCode",\(dfi){
#     v_range <- range(dfi$Uefeito)
#     teste <- rbind.fill(
#       select(dfi,-logOR),
#       filter(df_newpred,
#              Uefeito >= v_range[1] & Uefeito <= v_range[2])
#     ) %>% arrange(Uefeito) %>% 
#       mutate(forest_succession = dfi$forest_succession[1],
#              data_year = dfi$data_year[1],
#              lat = dfi$lat[1],
#              long = dfi$long[1],
#              SiteCode = dfi$SiteCode[1])
#   })
# }
# função que simula a predição a posteriori
f_calcPI <- \(gamm,
              nsim = 1000,
              to_exclude,
              simple_area=FALSE,
              quants=c(0.05,0.5,0.95)){
  #### 1a parte: somente efeito fixo ####
  # create the new data 
  df_newpred <- f_dfmd2(gamm$model)
  # quais componentes serão zerados?
  # to_exclude <- to_exclude[
  #   grep(pattern = paste(names(df_newpred),collapse = "|"),
  #        to_exclude)]
  # obtain the predictions
  df_pred <- f_predictions(gamm,nsim,to_exclude,df_newpred)
  # save the data frame
  l_df <- list()
  l_df$`apenas fixo` <- cbind(df_newpred,df_pred)
  #### 2a parte: efeito fixo e aleatório ####
  # take the original data
  df_newpred <- gamm$model
  # quais componentes serão zerados?
  to_exclude <- NULL
  # obtain the predictions
  df_pred <- f_predictions(gamm,nsim,to_exclude,df_newpred)
  # save the data frame
  l_df$`fixo e aleat` <- cbind(df_newpred,df_pred)
  # return
  return(l_df)
}
##################################################
# # # # # # # # # # # ROTINA # # # # # # # # # # #
##################################################
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
l_path <- list()
l_path$te <-  paste0(v_path,"rds/l_md_",c("areaperse","fragperse","fragtotal"),".rds")
df_tabsel <- readRDS("./5_resultados/df_tabsel_tehgam_efeitos.rds") %>% 
  relocate(efeito) %>% 
  filter(dAICc==0)
# if(length(unique(df_tabsel$modelo))==1){
#   l_md <- list()
#   for(i in l_path$te){
#     vefeito <- str_extract(i,pattern = "(?<=md\\_)(.*?)(?=\\.rds)")
#     l_md0 <- readRDS(i)
#     l_md[[vefeito]] <- l_md0[[ unique( df_tabsel$modelo ) ]]
#     rm(l_md0);gc()
#   }
# }
vlog <- lapply(l_path$te,\(i){
  hgam <- readRDS(i)
  hgam <- hgam[[ unique( df_tabsel$modelo ) ]]
  l_df <- f_calcPI(hgam,to_exclude = c("s(lat,long)","t2(Uefeito,k_z,SiteCode)"))
  saveRDS(l_df,gsub("l_md_","l_dfpred_",i))
  rm(hgam,l_df);gc()
})
# md_area <- readRDS(paste0(v_path,l_path$U))
# md_area <- md_area$`Área per se`$`s(land)|Site : gs`
# l_df <- f_calcPI(md_area,to_exclude = "s(Uefeito,SiteCode)",simple_area = TRUE)
# saveRDS(l_df,paste0(v_path,"rds/l_dfpred_areaperse_Ugs.rds"))