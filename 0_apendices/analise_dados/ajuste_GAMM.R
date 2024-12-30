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
library(mgcv)
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
# objetos comuns
df_md <- readRDS(file="dados/csv_SoE/df_logOR.rds") %>% 
  mutate(k_cont = as.numeric(as.character(k))) %>% 
  relocate(k,.after="Uefeito")
saveRDS(df_md,
        file=paste0(v_path,"rds/df_md.rds"))
##########################################################################
### Funções usadas para ajustar os 3 blocos de modelos: te, logU/U e k ###
##########################################################################
## Seguindo as especifícações de Pedersen et al. 2019 https://peerj.com/articles/6876/
# te: tensor entre o efeito na riqueza e a capacidade de dispersão per se
f_gam_te <- \(dfi){
  l_md <- list()
  l_md$`te(land)|Site : gs` <- gam(
    logOR ~ 
      te(Uefeito,k_cont,
         bs=c("cr","cr"),m=2,id = "efeito_comum") +
      t2(Uefeito,k_cont,SiteCode,
         bs=c("cr","cr","re"),m=2,full=TRUE,id = "efeito_sitio"),
    data=dfi,method = "REML")
  l_md$`te(land) + 1|Site` <- gam(
    logOR ~ 
      te(Uefeito,k_cont,
         bs=c("cr","cr"),m=2,id = "efeito_comum") +
      s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  l_md$`1 + 1|Site` <- gam(
    logOR ~ 1 + s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  saveRDS(l_md,
          file=paste0(v_path,"rds/l_md_",
                      gsub("Frag. total","fragtotal",dfi$contraste[1]) %>% 
                        gsub("Frag. per se","fragperse",.) %>% 
                        gsub("Área per se","areaperse",.),
                      ".rds"))
  rm(l_md);gc()
}
lapply(split(df_md,df_md$contraste),f_gam_te)
# s(k)+s(k)|Site
f_gam3 <- \(dfi){
  l_md <- list()
  l_md$`s(k)+s(k)|Site` <- gam(
    logOR ~ 
      s(k_cont,bs="cr",m=2,id = "efeito_comum") +
      s(k_cont, SiteCode, bs = "fs", xt=list(bs = "cr"), m=2, id="efeito_sitio"),
    data=dfi,method = "REML")
  l_md$`s(k)+1|Site` <- gam(
    logOR ~ 
      s(k_cont,bs="cr",m=2,id = "efeito_comum") +
      s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  return(l_md)
}
l_md <- dlply(df_md,"contraste",f_gam3)
saveRDS(l_md,paste0(v_path,"rds/l_md_onlyk.rds"))
rm("l_md");gc()
# s(logU/U)+s(logU/U)|Site
f_gam <- \(dfi,bs_type="cr"){
  l_md <- list()
  l_md$`s(land)|Site : gs` <- gam(
    logOR ~ 
      s(Uefeito,bs=bs_type,m=2, id="efeito_comum") +
      s(Uefeito, SiteCode, bs = "fs", xt=list(bs = bs_type), m=2, id="efeito_sitio"),
    data=dfi,method = "REML")
  l_md$`s(land) + 1|Site` <- gam(
    logOR ~ 
      s(Uefeito,bs=bs_type,m=2,id="efeito_comum") +
      s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  return(l_md)
}
l_md_logOR <- dlply(df_md,"contraste",f_gam)
saveRDS(l_md_logOR,file=paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))
rm(l_md_logOR);gc()