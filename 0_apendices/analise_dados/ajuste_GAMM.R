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
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
# objetos comuns
df_logOR <- read_csv(file="dados/csv/df_logOR.csv")
df_md <- df_logOR %>% select(-Uefeito,-(pristine:denominator),-matches("(k_z|^p)")) %>% 
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
saveRDS(select(df_md,contraste,SiteCode,k_cont,Uefeito,logOR),
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
lapply(split(df_md,contraste,df_md$contraste),f_gam_te)
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
  l_md$`1 + 1|Site` <- gam(
    logOR ~ 1 + s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  return(l_md)
}
l_md_logOR <- dlply(df_md,"contraste",f_gam)
saveRDS(l_md_logOR,file=paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))
###############################################
############ tabela de seleção  ###############
###############################################
l_path <- list()
l_path$te <-  paste0("rds/l_md_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$U <- "rds/l_md_simples_apudPedersen2019.rds"
l_path$k <- "rds/l_md_onlyk.rds"
# veffect <- l_path$te[[1]]
f_single_lmd <- \(veffect){
  # nome 
  vname <- str_extract(veffect,"(?<=l_md_)(.*?)(?=\\.rds)") %>% 
    gsub("areaperse","Área per se",.) %>% 
    gsub("fragperse","Frag. per se",.) %>% 
    gsub("fragtotal","Frag. total",.)
  # te
  lmd <- readRDS(paste0(v_path,veffect))
  # s(U)
  lmd_U <- readRDS(paste0(v_path,l_path$U))
  lmd_U <- lmd_U[[vname]]
  lmd$`s(logU/U)|Site : gs` <- lmd_U$`s(land)|Site : gs`
  lmd$`s(logU/U) + 1|Site` <- lmd_U$`s(land) + 1|Site`
  rm(lmd_U);gc()
  # s(k)
  lmd_k <- readRDS(paste0(v_path,l_path$k))
  lmd_k <- lmd_k[[vname]]
  lmd$`s(k)|Site : gs` <- lmd_k$`s(k)+s(k)|Site`
  lmd$`s(k) + 1|Site` <- lmd_k$`s(k)+1|Site`
  rm(lmd_k);gc()
  #
  df_tabsel <- f_TabSelGAMM(lmd)
  rm(lmd);gc()
  return(df_tabsel)
}
df_tabsel_geral <- alply(l_path$te,1,f_single_lmd)
names(df_tabsel_geral) <- str_extract(l_path$te,"(?<=l_md_)(.*?)(?=\\.rds)") %>% 
  gsub("areaperse","Área per se",.) %>% 
  gsub("fragperse","Frag. per se",.) %>% 
  gsub("fragtotal","Frag. total",.)
df_write <- lapply(names(df_tabsel_geral),\(li){
  mutate(df_tabsel_geral[[li]],contraste=li)
}) %>% do.call("rbind",.)
write_csv(df_write,paste0(v_path,"rds/df_tabsel_geral.csv"))
########################################################
##################### diagnósticos #####################
########################################################
df_tabsel <- read_csv(paste0(v_path,"rds/df_tabsel_geral.csv")) %>% 
  filter(dAICc==0)
f_diag_e_plots <- \(veffect){
  vname <- str_extract(veffect,"(?<=l_md_)(.*?)(?=\\.rds)") %>% 
    gsub("areaperse","Área per se",.) %>% 
    gsub("fragperse","Frag. per se",.) %>% 
    gsub("fragtotal","Frag. total",.)
  hgam <- readRDS(paste0(v_path,veffect))
  hgam <- hgam[[grep("gs",names(hgam),value=T)]]
  vpath <- f_diag(hgam,vname)
}
vlog <- lapply(l_path$te,f_diag_e_plots)




# l_md_logOR <- readRDS(file=paste0(v_path,"rds/l_md_simples.rds"))
df_tabelaSelecao <- ldply(l_md_logOR,f_TabSelGAMM,.id="pair")
write_csv(df_tabelaSelecao,
          file=paste0(v_path,"rds/tabsel_simples.csv"))
###  tabelas e diagnósticos
l_md <- readRDS(paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))
df_tabsel <- read_csv(paste0(v_path,"rds/tabsel_simples.csv")) %>%
  filter(dAICc==0)
l_md <- lapply(split(df_tabsel,df_tabsel$contraste),\(dfi){
  with(dfi,{l_md[[contraste]][[modelo]]})
})
# vpaths <- f_diagplots(l_md)
l_k.check <- lapply(l_md,k.check)

###################
#
#





#############################
l_md <- readRDS(paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))
l_md <- lapply(l_md, \(li){
  names(li) <- paste0("cr::",names(li))
  return(li)
})
l_md_logOR <- lapply(l_md_logOR,\(li){
  names(li) <- paste0("tp::",names(li))
  return(li)
})
l_md_geral <- lapply(names(l_md),\(li){
  c(l_md[[li]],
    l_md_logOR[[li]])
})
names(l_md_geral) <- names(l_md)
df_tabelaSelecao_geral <- ldply(l_md_geral,f_TabSelGAMM,.id="pair")
write_csv(df_tabelaSelecao_geral,
          file=paste0(v_path,"rds/tabsel_simples_tp_e_cr.csv"))
#



df_tabsel <- df_tabelaSelecao_geral %>% filter(dAICc==0)
l_md <- dlply(df_tabsel,"pair",\(dfi){
  with(dfi,{l_md_geral[[pair]][[modelo]]})
})
# vpaths <- f_diagplots(l_md)
l_k.check <- lapply(l_md,k.check)
vpaths <- f_diagplots(l_md)




