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
#################################################################
# função para criar o new data fixo
f_dfmd2 <- \(dff,length_pred = 100,site.posteriori ="SPigua1"){
  #0) dados empíricos
  df_ref <- read_csv(file="dados/csv_SoE/taxaU/df_contrastes.csv") %>% 
    select(SiteCode:k, contains("_logratio")) #%>% 
  names(df_ref) <- gsub("\\.","",names(df_ref)) %>% 
    gsub("\\_logratio","",.)
  df_ref <- pivot_longer(df_ref,c("area","fragperse","fragtotal"),
                         values_to="Uefeito",names_to="efeitos") %>% 
    inner_join(.,y=select(dff,SiteCode,forest_succession) %>% distinct())
  #1o criar predição para todo o intervalo de dados
  f_bypert2 <- \(dff){
    v_range <- range(df_ref$Uefeito)
    v_k <- range(df_ref$k)
    dfr <- expand.grid(Uefeito = seq(v_range[1],v_range[2],length.out=length_pred),
                k = seq(v_k[1],v_k[2],length.out=length_pred),
                SiteCode=site.posteriori)
    dfr <- adply(levels(dff$forest_succession),1,\(i){
      mutate(dfr,forest_succession=i)
    },.id=NULL)
  }
  df_pred_comum <- f_bypert2(dff)
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
  df_buf <- ddply(df_ref,c("efeitos","k","forest_succession"),f_dfbuffer)
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
    gam(Uefeito~s(k,bs="cr"),data=dfi)
  }
  l_md <- dlply(df_buf,
                c("efeitos","ext_class","forest_succession"),
                f_gam)
  ## fazer a predição média do modelo e guardar
  l_df_ref <- alply(distinct(select(df_buf,forest_succession,efeitos)),1,\(li){
    vname <- paste0(li$efeitos[1],c(".max.",".min."),li$forest_succession[1])
    lmd <- l_md[vname]
    #names(lmd) <- gsub(paste0(li,"."),"",names(lmd)) 
    df_ref <- lapply(names(lmd),\(i){
      md <- lmd[[i]]
      dfr <- filter(df_pred_comum,forest_succession==li$forest_succession[1])
      dfr[[ gsub("\\.","",str_extract(i,"\\.(.*?)\\.")) ]] <- predict.gam(md,dfr)
      return(dfr)
    }) %>% Reduce("inner_join",.)
  })
  names(l_df_ref) <- apply(distinct(select(df_buf,forest_succession,efeitos)),1,\(dfx){
   paste0(dfx["efeitos"],"_",dfx["forest_succession"]) 
  })
  ## filtrar para cada k os valores entre predições
  library(data.table)
  df_newpred <- lapply(names(l_df_ref),\(li){
    dfpred <- l_df_ref[[li]]
    setDT(dfpred)
    dfr <- as.data.frame(dfpred[,
                         .SD[.SD$min <= .SD$Uefeito & .SD$Uefeito <= .SD$max,
                             .(Uefeito,forest_succession,SiteCode)], 
                          by = k])
    mutate(dfr,efeito=str_split_1(li,"\\_")[1])
  }) %>% do.call("rbind",.) %>% 
    rename(Uref=Uefeito) %>% 
    mutate(k_z = f_z2(k,x_ref=df_ref$k),
           Uefeito = f_z2(Uref,x_ref=df_ref$Uefeito))
  ## 
  return(df_newpred)
}
dff <- gamm$model
df_p_pred <- f_dfmd2(dff)
saveRDS(df_p_pred,file="dados/csv_SoE/df_p_pred")
df_p_pred <- readRDS(file="dados/csv_SoE/df_p_pred")

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
# função que simula a predição a posteriori
f_calcPI <- \(gamm,
              dfppred=df_p_pred,
              vefeito,
              nsim = 1000,
              to_exclude,
              simple_area=FALSE,
              quants=c(0.05,0.5,0.95)){
  #### 1a parte: somente efeito fixo ####
  # create the new data 
  df_newpred <- filter(dfppred,efeito==vefeito)
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
##### retorno dos dados para a escala padrão
f_z <- \(x) (x-mean(x))/sd(x)
f_z2 <- \(x,x_ref) (x-mean(x_ref))/sd(x_ref)
f_anti_z <- \(x,x_ref) (x * sd(x_ref)) + mean(x_ref)
f_escalaoriginal <- \(ipath){
  # quais efeitos não serão lidos?
  efeitos_rm = str_extract(ipath,"(?<=pred\\_)(.*?)(?=\\.rds)") %>% 
    gsub("areaperse","area",.) %>% 
    setdiff(c("area","fragperse","fragtotal"),.)
  # lista com os df de predição
  l_dfpred <- readRDS(ipath)
  # juntando o fixo e aleat com os valores empíricos
  df_ref <- read_csv(file="dados/csv_SoE/taxaU/df_contrastes.csv") %>% 
    select(SiteCode:k, contains("_logratio")) #%>% 
  names(df_ref) <- gsub("\\.","",names(df_ref)) %>% 
    gsub("\\_logratio","",.)
  df_ref <- pivot_longer(df_ref,c("area","fragperse","fragtotal"),
                         values_to="Uefeito",names_to="efeitos")
  vUref <- df_ref$Uefeito
  vkref <- df_ref$k
  ## apenas fixo
  l_dfpred[["apenas fixo"]] <- l_dfpred[["apenas fixo"]] %>% 
    mutate(Uref = f_anti_z(x = Uefeito, x_ref = vUref),
           k = f_anti_z(x = k_z, x_ref = vkref))
  # salvamento e limpeza
  saveRDS(l_dfpred,ipath)
  rm(l_dfpred);gc()
}
##################################################
# # # # # # # # # # # ROTINA # # # # # # # # # # #
##################################################
# caminho
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
# modelos ajustados
l_path <- list()
l_path$te <-  paste0(v_path,"rds/l_md_",c("areaperse","fragperse","fragtotal"),".rds")
# tabela de seleção
df_tabsel <- readRDS("./5_resultados/df_tabsel_tehgam_efeitos.rds") %>% 
  relocate(efeito) %>% 
  filter(dAICc==0)
# predição a posteriori do mais plausível
if(FALSE){
  vlog <- lapply(l_path$te,\(i){
    hgam <- readRDS(i)
    hgam <- hgam[[ unique( df_tabsel$modelo ) ]]
    vefeito <- str_extract(i,"(?<=md\\_)(.*?)(?=\\.rds)") %>% 
      gsub("areaperse","area",.)
    l_df <- f_calcPI(hgam,vefeito = vefeito,to_exclude = c("s(lat,long)","t2(Uefeito,k_z,SiteCode)"))
    saveRDS(l_df,gsub("l_md_","l_dfpred_",i))
    rm(hgam,l_df);gc()
  })
}
# retornar para a escala padrão
## df valores de referência e padronização
df_contrastes <- read_csv(file="dados/csv_SoE/taxaU/df_contrastes.csv") %>% 
  select(SiteCode:k, contains("_logratio")) %>% 
  mutate(k_z=f_z(k)) %>% relocate(k_z,.after="k")
names(df_contrastes) <- gsub("\\.","",names(df_contrastes)) %>% 
  gsub("\\_logratio","",.)
## caminhos para os rds e função para padronização
l_rds <- gsub("l_md_","l_dfpred_",l_path$te)
lapply(l_rds,f_escalaoriginal)