####### Sumário da congruência por paisagem hipotética 
library(plyr)
library(dplyr)
library(mgcv)
library(lme4)
library(DHARMa)
library(ggplot2)
library(grid)
library(gridExtra)
df_adSAD <- readRDS("5_resultados/df_adSAD.rds")
df_nSAD <- adply(unique(df_adSAD$land_type),1,\(i){
  df_return <- df_adSAD %>% filter(taxaU==i & land_type==i) %>% 
    select(-land_type,-c(Smed:Smax)) %>% 
    rename("land"="taxaU",
           "nSAD"="nCongKS")
},.id=NULL) %>% 
  mutate(k=round(k,2),across(-nSAD,factor))
f_glmer <- \(vf,dfi){
  glmer(formula=vf,
        family='binomial',
        data=dfi,
        control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=2e9)))
}
l_f <- list()
l_f$`land * k + land|Site` <- cbind(nSAD,100-nSAD) ~ land * k + land|SiteCode
l_f$`land + k + land|Site` <- cbind(nSAD,100-nSAD) ~ land + k + land|SiteCode 
l_f$`land + land|Site` <- cbind(nSAD,100-nSAD) ~ land + land|SiteCode
l_f$`k + 1|Site` <- cbind(nSAD,100-nSAD) ~ k + 1|SiteCode 
l_f$`1 + 1|Site` <- cbind(nSAD,100-nSAD) ~ 1 + 1|SiteCode
l_md <- lapply(l_f,f_glmer,dfi=df_nSAD)