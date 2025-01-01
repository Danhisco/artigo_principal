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
###############################################
############ tabela de seleção  ###############
###############################################
l_path <- list()
l_path$te <-  paste0("rds/l_md_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$U <- "rds/l_md_simples_apudPedersen2019.rds"
l_path$k <- "rds/l_md_onlyk.rds"
# veffect <- l_path$te[[1]]
f_single_lmd <- \(veffect,v_path){
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
formals(f_single_lmd)$v_path <- v_path
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
l_path <- list()
l_path$te <-  paste0("rds/l_md_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$U <- "rds/l_md_simples_apudPedersen2019.rds"
l_path$k <- "rds/l_md_onlyk.rds"
df_tabsel <- read_csv(paste0(v_path,"rds/df_tabsel_geral.csv")) %>% 
  filter(dAICc==0)
f_diag_e_plots <- \(veffect,
                    pattern_extract="(?<=l_md_)(.*?)(?=\\.rds)"){
  vname <- str_extract(veffect,pattern_extract) %>% 
    gsub("areaperse","Área per se",.) %>% 
    gsub("fragperse","Frag. per se",.) %>% 
    gsub("fragtotal","Frag. total",.)
  hgam <- readRDS(paste0(v_path,veffect))
  hgam <- hgam[[grep("gs",names(hgam),value=T)]]
  vpath <- f_diag(hgam,vname)
  rm(hgam);gc()
}
f_diag_maisplaus <- \(dfi){
  vname <- dfi$contraste
  hgam <- dfi$modelo
  if(grepl("logU",hgam)){
    lmd <- readRDS(paste0(v_path,l_path$U))
    lmd <- lmd[[vname]]
    lmd <- lmd[[grep("gs",names(lmd),value=TRUE)]]
    vpath <- f_diag(lmd,vname)
  }else if(grepl("s\\(k",hgam)){
    l_path$k
  }else if(grepl("te\\(land",hgam)){
    vp <- gsub("Área per se","areaperse",vname) %>% 
      gsub("Frag. per se","fragperse",.) %>% 
      gsub("Frag. total","fragtotal",.)
    vp <- grep(vp,l_path$te,value=TRUE)
    md <- readRDS(paste0(v_path,vp))
    md <- md[[hgam]]
    vpath <- f_diag(md,vname)
  }
}
vlog <- lapply(l_path$te,f_diag_maisplaus)
# 
# f_diag_e_plots2 <- \(vpath,
#                      efeitos_paisagem="Área per se",
#                      pattern_extract="(?<=l_md_)(.*?)(?=\\.rds)"){
#   hgam <- readRDS(vpath)[[efeitos_paisagem]]
#   hgam <- hgam[[grep("gs",names(hgam),value=T)]]
#   vpath <- f_diag(hgam,efeitos_paisagem)
#   rm(hgam);gc()
# }
# vpath <- list.files(path = paste0(v_path,"rds"),
#                     pattern = "Pedersen2019",
#                     full.names = TRUE)[4]
# 
# 
# 
# 
# # l_md_logOR <- readRDS(file=paste0(v_path,"rds/l_md_simples.rds"))
# df_tabelaSelecao <- ldply(l_md_logOR,f_TabSelGAMM,.id="pair")
# write_csv(df_tabelaSelecao,
#           file=paste0(v_path,"rds/tabsel_simples.csv"))
# ###  tabelas e diagnósticos
# l_md <- readRDS(paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))
# df_tabsel <- read_csv(paste0(v_path,"rds/tabsel_simples.csv")) %>%
#   filter(dAICc==0)
# l_md <- lapply(split(df_tabsel,df_tabsel$contraste),\(dfi){
#   with(dfi,{l_md[[contraste]][[modelo]]})
# })
# # vpaths <- f_diagplots(l_md)
# l_k.check <- lapply(l_md,k.check)
# 
# ###################
# #
# #
# 
# 
# 
# 
# 
# #############################
# l_md <- readRDS(paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))
# l_md <- lapply(l_md, \(li){
#   names(li) <- paste0("cr::",names(li))
#   return(li)
# })
# l_md_logOR <- lapply(l_md_logOR,\(li){
#   names(li) <- paste0("tp::",names(li))
#   return(li)
# })
# l_md_geral <- lapply(names(l_md),\(li){
#   c(l_md[[li]],
#     l_md_logOR[[li]])
# })
# names(l_md_geral) <- names(l_md)
# df_tabelaSelecao_geral <- ldply(l_md_geral,f_TabSelGAMM,.id="pair")
# write_csv(df_tabelaSelecao_geral,
#           file=paste0(v_path,"rds/tabsel_simples_tp_e_cr.csv"))
# #
# 
# 
# 
# df_tabsel <- df_tabelaSelecao_geral %>% filter(dAICc==0)
# l_md <- dlply(df_tabsel,"pair",\(dfi){
#   with(dfi,{l_md_geral[[pair]][[modelo]]})
# })
# # vpaths <- f_diagplots(l_md)
# l_k.check <- lapply(l_md,k.check)
# vpaths <- f_diagplots(l_md)
# 
# 
# 
# 
