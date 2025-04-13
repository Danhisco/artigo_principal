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
##### retorno dos dados para a escala padrão
f_z <- \(x) (x-mean(x))/sd(x)
f_anti_z <- \(x,x_ref) (x * sd(x_ref)) + mean(x_ref)
f_escalaoriginal <- \(ipath){
  # quais efeitos não serão lidos?
  efeitos_rm = str_extract(ipath,"(?<=pred\\_)(.*?)(?=\\.rds)") %>% 
    gsub("areaperse","area",.) %>% 
    setdiff(c("area","fragperse","fragtotal"),.)
  # lista com os df de predição
  l_dfpred <- readRDS(ipath)
  # juntando o fixo e aleat com os valores empíricos
  ## fixo e aleatório
  l_dfpred[["fixo e aleat"]] <- 
    inner_join(l_dfpred[["fixo e aleat"]],
               select(df_contrastes,-all_of(efeitos_rm)) %>% 
                 rename(Uref = last_col()),
               by=c("SiteCode","k_z"))
  vUref <- l_dfpred[["fixo e aleat"]]$Uref
  vkref <- l_dfpred[["fixo e aleat"]]$k
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
    l_df <- f_calcPI(hgam,to_exclude = c("s(lat,long)","t2(Uefeito,k_z,SiteCode)"))
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