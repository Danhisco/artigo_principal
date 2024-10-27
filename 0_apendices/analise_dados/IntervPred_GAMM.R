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


################# remoção da extrapolação
library(mgcv)
## funções de ajuste e de plot
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
# objetos comuns
df_logOR <- read_csv(file="dados/csv/df_logOR.csv")
df_md <- df_logOR %>% 
  select(-Uefeito,-(pristine:denominator),-matches("(k_z|^p)")) %>% 
  rename(logOR = itself) %>% 
  filter(forest_succession!="capoeira") %>% 
  mutate(contraste =  gsub("non_frag.ideal","Área per se",contraste) %>% 
           gsub("contemp-non_frag","Frag. per se",.) %>% 
           gsub("contemp-ideal","Frag. total",.),
         across(c(contraste:SiteCode,forest_succession,k),factor))
df_contrastes <- read_csv("dados/csv/taxaU/df_contrastes.csv") %>% 
  select(SiteCode,k,ends_with("logratio")) %>% 
  pivot_longer(cols = ends_with("logratio"),
               names_to = "contraste",
               values_to = "Uefeito") %>% 
  mutate(contraste = gsub("area_logratio","Área per se",contraste) %>% 
           gsub("frag.perse_logratio","Frag. per se",.) %>% 
           gsub("frag.total_logratio","Frag. total",.),
         across(SiteCode:contraste,factor))
df_md <- inner_join(df_md,df_contrastes,by=c("SiteCode","contraste","k")) %>% 
  relocate(Uefeito,.after="logOR") %>% 
  mutate(k_cont = as.numeric(as.character(k)))
#
df_buf <- ddply(df_md,c("contraste","k"),\(dfi){
  v_range <- range(dfi$Uefeito)
  v_quant <- quantile(dfi$Uefeito,probs=c(0.01,0.99))
  filter(dfi,Uefeito%in%v_range) %>% 
    select(-logOR,-c(forest_succession:data_year),-SiteCode) %>% 
    mutate(ext_class=ifelse(Uefeito==v_range[1],"min","max"),
           Uefeito_quant=ifelse(ext_class=="min",v_quant[1],v_quant[2])) %>% 
    relocate(ext_class,.after=contraste) %>% 
    relocate(Uefeito_quant,.after=Uefeito)
}) %>% select(-k)
df_buf %>% 
  pivot_longer(starts_with("Uefeito")) %>% 
  filter(name=="Uefeito") %>% 
  ggplot(aes(x=k_cont,y=value,color=ext_class)) +
  geom_point() +
  geom_line() +
  geom_smooth(method="gam",se=FALSE) +
  scale_color_manual(values=c("darkred","darkgreen")) +
  facet_wrap(~contraste,ncol=3)
#
f_gam <- \(dfi){
  gam(Uefeito~s(k_cont,bs="cr"),data=dfi)
}
l_md <- dlply(df_buf,c("contraste","ext_class"),f_gam)
for(i in names(l_md)){
  p <- gratia::draw(l_md[[i]],residuals = TRUE) + labs(title=i)
  print(p)
  if(!isTRUE(askYesNo("Do you want to see the residuals?"))) break
  p <- gratia::appraise(l_md[[i]])
  print(p)
  if(!isTRUE(askYesNo("next model?"))) break
}
saveRDS(l_md,file=paste0(v_path,"rds/l_md_refU.rds"))
# i) carregar os dados de l_df_pred
l_paths <- paste0(v_path,"rds/l_dfpred_",c("fragtotal","fragperse","areaperse"),".rds")
l_df_pred <- lapply(l_paths,readRDS) %>% 
  lapply(.,"[[","apenas fixo")
names(l_df_pred) <- c("Frag. total","Frag. per se","Área per se")
# ii) filtrar os valores únicos de k em um novo data frame
l_df_ref <- lapply(l_df_pred,select,k_cont,SiteCode) %>% lapply(.,distinct)
# iii) fazer a predição média do modelo para cada ponto e guardar em um dataframe
l_df_ref <- lapply(names(l_df_ref),\(li){
  lmd <- l_md[grep(li,names(l_md))]
  names(lmd) <- gsub(paste0(li,"."),"",names(lmd)) 
  df_ref <- lapply(names(lmd),\(i){
    md <- lmd[[i]]
    dfr <- l_df_ref[[li]]  
    dfr[[i]] <- predict.gam(md,dfr)
    return(dfr)
  }) %>% Reduce("inner_join",.)
})
names(l_df_ref) <- c("Frag. total","Frag. per se","Área per se")
# iv) filtrar para cada k os valores entre as predições
library(data.table)
l_df_new <- lapply(names(l_df_ref),\(li){
  df_pred <- l_df_pred[[li]]
  df_ref <- l_df_ref[[li]]
  setDT(df_pred)
  alply(df_ref,1,\(dfi){
    maxv <- dfi[["max"]]
    minv <- dfi[["min"]]
    kv <- dfi[["k_cont"]]
    df_pred[k_cont==kv & between(Uefeito,minv,maxv),]
  }) %>% rbindlist
})
names(l_df_new) <- names(l_df_ref)
# salvamento
l_paths <- gsub("dfpred","dfnew",l_paths)
names(l_paths) <- names(l_df_ref)
lapply(names(l_df_ref),\(i){
  dfrds <- l_df_new[[i]]
  vpath <- l_paths[i]
  saveRDS(dfrds,file=vpath)
})
