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
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
# caminhos para os modelos
# paths <- list.files(path=paste0(v_path,"rds"),
#                     pattern="l_md_3aperg",full.names = T)
# names(paths) <- str_extract(paths,"(?<=cov\\_)(.*?)(?=\\.rds)")
#
l_path <- list()
l_path$te <-  paste0("rds/l_md_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$U <- "rds/l_md_simples_apudPedersen2019.rds"
df_tabsel <- read_csv(paste0(v_path,"rds/df_tabsel_geral.csv"))
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
f_dfmd2 <- \(dff,length_pred = 100,site.posteriori ="SPigua1",mdarea=FALSE){
  if(mdarea){
    v_range <- range(dff$Uefeito)
    data.frame(Uefeito = seq(v_range[1],v_range[2],length.out=length_pred),
               SiteCode=site.posteriori)
  }else{
    v_range <- range(dff$Uefeito)
    v_k <- range(dff$k_cont)
    expand.grid(Uefeito = seq(v_range[1],v_range[2],length.out=length_pred),
                k_cont = seq(v_k[1],v_k[2],length.out=length_pred),
                SiteCode=site.posteriori)  
  }
}
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
  df_newpred <- f_dfmd2(gamm$model,mdarea = simple_area)
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
# rotina
vlog <- lapply(l_path$te,\(i){
  hgam <- readRDS(paste0(v_path,i))
  hgam <- hgam[[grep("gs",names(hgam))]]
  l_df <- f_calcPI(hgam,to_exclude = "t2(Uefeito,k_cont,SiteCode)")
  saveRDS(l_df,paste0(v_path,gsub("l_md_","l_dfpred_",i)))
  rm(hgam,l_df);gc()
})
md_area <- readRDS(paste0(v_path,l_path$U))
md_area <- md_area$`Área per se`$`s(land)|Site : gs`
l_df <- f_calcPI(md_area,to_exclude = "s(Uefeito,SiteCode)",simple_area = TRUE)
saveRDS(l_df,paste0(v_path,"rds/l_dfpred_areaperse_Ugs.rds"))

# l_df_pred <- lapply(l_md,f_calcPI)
# saveRDS(l_df_pred,paste0(v_path,"rds/l_dfpred_simples_apudPedersen2019_tp.rds"))


