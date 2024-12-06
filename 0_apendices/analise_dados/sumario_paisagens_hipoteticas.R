####### Sumário da congruência por paisagem hipotética 
library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(mgcv)
library(lme4)
library(DHARMa)
library(ggplot2)
library(grid)
library(gridExtra)
library(bbmle)
source("source/GAMMtools.R")
df_coord <- read_csv(file = "dados/df_dados_disponiveis.csv") %>% 
  mutate(lat = ifelse(is.na(lat_correct),lat,lat_correct),
         long = ifelse(is.na(long_correct),long,long_correct),
         Sitecode = factor(SiteCode)) %>% 
  select(SiteCode,lat,long)
df_adSAD <- readRDS("5_resultados/df_adSAD.rds")
df_nSAD <- adply(unique(df_adSAD$land_type),1,\(i){
  df_return <- df_adSAD %>% filter(taxaU==i & land_type==i) %>% 
    select(-land_type,-c(Smed:Smax)) %>% 
    rename("land"="taxaU",
           "nSAD"="nCongKS")
},.id=NULL) %>% 
  mutate(k=round(k,2),
         across(c(land,SiteCode,k),factor)) %>% 
  inner_join(df_coord)
df_nSAD$SiteCode <- factor(df_nSAD$SiteCode)
f_gam <- \(vf,dfi){
  gam(formula=vf,
      family='binomial',
      data=dfi,
      method="REML")
}
l_f <- list()
l_f$`land * k + land|Site` <- cbind(nSAD,100-nSAD) ~ land * k + s(SiteCode,by=land,bs="re")
l_f$`land + k + land|Site` <- cbind(nSAD,100-nSAD) ~ land + k + s(SiteCode,by=land,bs="re")
l_f$`land + land|Site` <- cbind(nSAD,100-nSAD) ~ land + s(SiteCode,by=land,bs="re")
l_f$`k + 1|Site` <- cbind(nSAD,100-nSAD) ~ k + s(SiteCode,bs="re") 
l_f$`1 + 1|Site` <- cbind(nSAD,100-nSAD) ~ 1 + s(SiteCode,bs="re")
# l_f$`land * k + s(lat,long) + land|Site` <- cbind(nSAD,100-nSAD) ~ land * k + s(lat,long) + s(SiteCode,by=land,bs="re")
# l_f$`land + s(lat,long) + k + s(lat,long) + land|Site` <- cbind(nSAD,100-nSAD) ~ land + s(lat,long) + k + s(lat,long) + s(SiteCode,by=land,bs="re")
# l_f$`land + s(lat,long) + land|Site` <- cbind(nSAD,100-nSAD) ~ land + s(lat,long) + s(SiteCode,by=land,bs="re")
# l_f$`k + s(lat,long) + 1|Site` <- cbind(nSAD,100-nSAD) ~ k + s(lat,long) + s(SiteCode,bs="re") 
# l_f$`1 + s(lat,long) + 1|Site` <- cbind(nSAD,100-nSAD) ~ 1 + s(lat,long) + s(SiteCode,bs="re")
doMC::registerDoMC(3)
l_md <- llply(l_f,f_gam,dfi=df_nSAD,.parallel = TRUE)
df_tabsel <- f_TabSelGAMM(l_md)
write_csv(df_tabsel,"1_to_compile_dissertacao_EM_USO/00_Resultados/tabelas/tabselecao_sumario_paisagens.csv")
