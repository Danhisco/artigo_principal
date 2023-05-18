library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(MuMIn)
library(AICcmodavg)
library(insight)
library(bbmle)
library(DHARMa)
library(mgcv)
library(lme4)
library(gamm4)
library(plyr)
library(dplyr)
df_md <- read_csv("dados/csv/df_md_logOR.csv") |> 
  mutate(SiteCode = factor(SiteCode))
#
f_gam_logOR <- function(df){
  l_md <- list()
  l_md[[1]] <- gam(logOR ~ 
                     te(p_z,k_z,k=6) +
                     t2(k_z,SiteCode,k=6,bs=c("cr","re"),full=TRUE),
                   data=df, method="REML")
  l_md[[2]] <- gam(logOR ~ 
                     te(cubeR_efeito_area_z,cubeR_efeito_frag_z,k=6) +
                     t2(cubeR_efeito_area_z,cubeR_efeito_frag_z,SiteCode,k=6,bs=c("cr","cr","re"),full=TRUE),
                   data=df, method="REML")
  names(l_md) <- c("~ f(p,k)","~ f(efeito area, efeito frag)")
  return(l_md)
}
l_md.logOR <- df_md |> 
  filter(logOR_pair == "cont.non_frag") |> 
  f_gam_logOR()
save(l_md.logOR,file="dados/Rdata/l_md.logOR.Rdata")