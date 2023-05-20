library(readr)
library(purrr)
library(stringr)
library(tidyr)
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
                     s(k_z,bs = "cr") + s(p_z,bs = "cr") + ti(p_z,k_z) +
                     s(k_z, SiteCode, bs="fs",xt=list(bs="cr")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[2]] <- gam(logOR ~ 
                     s(r3_area_z,bs = "cr") + s(r3_frag_z,bs = "cr") + ti(r3_area_z,r3_frag_z) +
                     s(r3_frag_z, SiteCode, bs="fs",xt=list(bs="cr")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[3]] <- gam(logOR ~ 
                     s(r3_area_z,bs = "cr") + s(r3_frag_z,bs = "cr") + ti(r3_area_z,r3_frag_z) +
                     s(r3_area_z, SiteCode, bs="fs",xt=list(bs="cr")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[4]] <- gam(logOR ~ 
                     s(r3_area_z,bs = "cr") + s(r3_frag_z,bs = "cr") + ti(r3_area_z,r3_frag_z) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[5]] <- gam(logOR ~ 
                     s(r3_frag_z,bs = "cr") +
                     s(r3_frag_z, SiteCode, bs="fs",xt=list(bs="cr")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[6]] <- gam(logOR ~ 
                     s(r3_area_z,bs = "cr") +
                     s(r3_area_z, SiteCode, bs="fs",xt=list(bs="cr")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[7]] <- gam(logOR ~ 
                     s(r3_area_z,bs = "cr") + s(r3_frag_z,bs = "cr") +
                     s(r3_frag_z, SiteCode, bs="fs",xt=list(bs="cr")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[8]] <- gam(logOR ~ 
                     s(r3_area_z,bs = "cr") + s(r3_frag_z,bs = "cr") +
                     s(r3_area_z, SiteCode, bs="fs",xt=list(bs="cr")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  names(l_md) <- c("~ f(p, k, k|Site)",
                   "~ f(area, frag, frag|Site)","~ f(area, frag, area|Site)","~ f(area, frag, 1|Site)",
                   "~ f(frag, frag|Site)","~ f(area, area|Site)",
                   "~ f(area,frag,frag|Site - ti","~ f(area,frag,area|Site - ti")
  return(l_md) 
}
l_md.logOR <- df_md |> 
  filter(logOR_pair == "cont.non_frag") |> 
  f_gam_logOR()
save(l_md.logOR,file="dados/Rdata/l_md.logOR.Rdata")
