# pacotes e dados
library(readr)
library(stringr)
library(tidyr)
library(insight)
library(bbmle)
library(DHARMa)
library(mgcv)
library(lme4)
library(plyr)
library(dplyr)
source("source/nameModel.R")
source("source/f_PredIntGAMM.R")
# objetos
# load("dados/Rdata/md_MNEE.Rdata")
# load("dados/Rdata/l_md_GAMM.Rdata")
# dados
load("dados/Rdata/l_md.contrastes.Rdata")
# df_p
df_p <- read_csv("dados/df_p.csv")
# df_resultMNEE
df_resultados <- read_csv("dados/csv/resultados_MN/df_resultados.csv")
# df_sim
df_sim <- read_csv("dados/df_simulacao.csv")
# df_ad
f_z <- function(x) (x-mean(x))/sd(x)
df_ad <- inner_join(x=distinct(select(df_resultados,-(Ssd:Smax))),
                    y=distinct(select(df_sim,-tif.path)),
                    by=c("SiteCode","k"),multiple="all") |>
  # escolhi remover depois do ajuste, para obter um melhor ajuste e depois excluir os valores
  # filter(k > 0.20) |> # figuras 5 e 6 apêndice Efeito de Escala
  mutate(k_cont = round(k,2),
         across(Ntotal:S_obs,log,.names="log_{.col}"),
         log_S.N = log_S_obs - log_Ntotal,
         across(c(p,k_cont,log_Ntotal:log_S.N),f_z,.names = "{.col}_z")) |>
  select(-c(effort_ha, Ntotal:S_obs))
# df_md
df_md <- df_ad |> 
  mutate(SiteCode = factor(SiteCode))
# df_newdata
df_newpred <- read_csv(file="dados/csv/df_newpred.csv") |> 
  select(-starts_with("log_")) |> 
  rename(k_z = k_cont_z)
# df_pred
# df_pred <- read_csv("dados/csv/resultados_MN/MNEE/df_pred.csv")
# df_newdat <- expand.grid(p_z = seq(min(df_ad$p_z),max(df_ad$p_z), length=150),
#                          k_cont_z = seq(min(df_ad$k_cont_z),max(df_ad$k_cont_z), length=150))
# df_newdat <- adply(df_newdat,1,.fun = \(x) cbind(x,distinct(select(df_ad,SiteCode,log_S_obs_z,log_Ntotal_z))))
# write_csv(df_newdat,"dados/csv/df_newdataSADs.csv")
df_contrastes <- read_csv(file="dados/csv/taxaU/df_contrastes.csv") |> 
  inner_join(df_sim |> select(SiteCode,Ntotal:S_obs) |> distinct(),
             by="SiteCode") |> 
  rename(N=Ntotal,S=S_obs) |> 
  mutate(across(N:S,log,.names="log.{.col}"),
         across(c(p,k,log.N:log.S),f_z,.names = "{.col}_z"),
         SiteCode = factor(SiteCode)) |> 
  select(-c(N:log.S))
#
df_intPred <- ldply(l_md.contrastes,f_PredInt.GAMM)