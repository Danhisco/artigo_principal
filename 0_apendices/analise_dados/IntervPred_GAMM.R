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
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
# caminhos para os modelos
paths <- list.files(path=paste0(v_path,"rds"),
                    pattern="l_md_3aperg",full.names = T)
names(paths) <- str_extract(paths,"(?<=cov\\_)(.*?)(?=\\.rds)")
#################################################################
# função para criar o new data fixo
f_dfmd <- \(dff,byforest,length_pred = 150,site.posteriori ="SPigua1"){
  fdfmd <- \(dff){
    v_range <- range(dff$Uefeito)
    df_newpred <- select(dff,-logOR,-Uefeito) %>% 
      mutate(SiteCode=site.posteriori) %>% head(n=1)
    try({
      df_newpred <- cbind(
        data.frame(Uefeito = seq(v_range[1],v_range[2],length.out=length_pred)),
        df_newpred)
    },silent = TRUE)
    return(df_newpred)
  }
  if(byforest){
    df_return <- ddply(dff,"forest_succession",fdfmd)
  }else{
    df_return <- fdfmd(dff)
  }
  return(df_return)
}
# função que realiza a predição a posteriori dado o modelo, new data e o vetor de exclusão de componentes
f_predictions <- \(gamm,nsim,to_exclude,df_newpred){
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
              quants=c(0.05,0.5,0.95)){
  #### 1a parte: somente efeito fixo ####
  # create the new data 
  df_newpred <- f_dfmd(
    gamm$model,
    grepl("by = fore",as.character(formula(gamm)[[3]])[2]))
  # quais componentes serão zerados?
  to_exclude <- to_exclude[
    grep(pattern = paste(c("(Intercept)",names(df_newpred)),collapse = "|"),
         to_exclude)]
  # obtain the predictions
  df_pred <- f_predictions(gamm,nsim,to_exclude,df_newpred)
  # save the data frame
  l_df <- list()
  l_df$`apenas fixo` <- cbind(df_newpred,df_pred)
  #### 2a parte: efeito fixo e aleatório ####
  # take the original data
  df_newpred <- gamm$model
  # quais componentes serão zerados?
  to_exclude <- to_exclude[-grep("Uefeito,SiteCode",to_exclude)]
  # obtain the predictions
  df_pred <- f_predictions(gamm,nsim,to_exclude,df_newpred)
  # save the data frame
  l_df$`fixo e aleat` <- cbind(df_newpred,df_pred)
  # return
  return(l_df)
}
f_df_pred <- \(vpath){
  l_md <- readRDS(vpath)
  # esses serão os splines desconsiderados para construir o fixo, e depois com o aleat
  formals(f_calcPI2)$to_exclude = c("s(Uefeito,SiteCode)",
                                    "(Intercept)",
                                    "forest_successionprimary/secondary",
                                    "forest_successionsecondary",
                                    "s(SiteCode)",
                                    "s(lat,long)",
                                    "s(data_year)")
  # rotina; detalhes em source/GAMMtools.R
  lapply(l_md,f_calcPI2)
}
# a predição de todos os modelos candidatos
l_df_pred <- lapply(paths,f_df_pred)
saveRDS(l_df_pred,paste0(v_path,"rds/l_dfpredictions_fromfixedrandom_landeffect.rds"))
# calculo do model averaging 
l_df_pred <- readRDS(paste0(v_path,"rds/l_dfpredictions_fromfixedrandom_landeffect.rds"))
## tabela de seleção da pergunta: "Quais as covariáveis adequadas?"
df_tabsel <- read_csv(paste0(v_path,"rds/tabsel_3aperg_quais_cov.csv")) %>% 
  relocate(contraste) %>% select(contraste:modelo,weight) %>% 
  mutate(nome = gsub("Área per se","areaperse",contraste) %>% 
           gsub("Frag. per se","fragperse",.) %>% 
           gsub("Frag. total","fragtotal",.))
l_df_avgpred <- dlply(df_tabsel,"contraste",\(dff){
  ## multiplicação pelo peso de evidência 
  l_df <- alply(dff,1,\(dfi){
    with(dfi,{
      l_fixrand <- l_df_pred[[nome]][[modelo]]
      l_fixrand <- lapply(l_fixrand,\(li){
        mutate(li,
               across(starts_with("Q_"),~.x*weight))
      })
      return(l_fixrand)})
    })
  names(l_df) <- names(l_df_pred[[dff$nome[1]]])
  ## soma para obter a média
  l_df_return <- lapply(sapply(l_df,names)[,1],\(li){ # os datasets 'apenas fixo' e 'fixo e aleatório'
    ldf <- lapply(l_df,\(i) i[[li]]) # os datasets de um único tipo (e.g. apenas fixo)
    vquantil <- select(ldf[[1]],starts_with("Q_")) %>% names # colunas que serão somadas
    names(vquantil) <- vquantil
    # soma per se e retorno
    lapply(vquantil,\(x1){
      Reduce("+",lapply(ldf,\(x) x[[x1]]))
    }) %>% do.call("cbind",.) %>% 
      cbind(select(ldf[[1]],-all_of(vquantil)),.)
  })
  names(l_df_return) <- sapply(l_df,names)[,1]
  return(l_df_return)
})
saveRDS(l_df_avgpred,file=paste0(v_path,"rds/l_df_avgpred.rds"))