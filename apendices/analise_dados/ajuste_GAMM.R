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
  mutate(SiteCode = factor(SiteCode)) |> 
  filter(logOR_pair == "cont.non_frag")
# não há capacidade computacional suficiente para rodar:
# md.logOR_t2full <- gam(logOR ~ 
#                          te(area_z,frag_z,k=15,bs=c("tp","tp"),m=2) +
#                          t2(area_z,frag_z,SiteCode,k=6,bs=c("tp","tp","re"),full=TRUE,m=2),
#                        data=df, method="REML")
f_gam_logOR <- function(df){
  l_md <- list()
  l_md[[1]] <- gam(logOR ~ 
                     s(k_z,bs = "tp") + s(p_z,bs = "tp") + ti(p_z,k_z) +
                     s(k_z, SiteCode, bs="fs",xt=list(bs="tp")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[2]] <- gam(logOR ~ 
                     s(area_z,bs = "tp") + s(frag_z,bs = "tp") + ti(area_z,frag_z) +
                     s(frag_z, SiteCode, bs="fs",xt=list(bs="tp")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[3]] <- gam(logOR ~ 
                     s(area_z,bs = "tp") + s(frag_z,bs = "tp") + ti(area_z,frag_z) +
                     s(area_z, SiteCode, bs="fs",xt=list(bs="tp")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[4]] <- gam(logOR ~ 
                     s(area_z,bs = "tp") + s(frag_z,bs = "tp") + ti(area_z,frag_z) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[5]] <- gam(logOR ~ 
                     s(frag_z,bs = "tp") +
                     s(frag_z, SiteCode, bs="fs",xt=list(bs="tp")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[6]] <- gam(logOR ~ 
                     s(area_z,bs = "tp") +
                     s(area_z, SiteCode, bs="fs",xt=list(bs="tp")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[7]] <- gam(logOR ~ 
                     s(area_z,bs = "tp") + s(frag_z,bs = "tp") +
                     s(frag_z, SiteCode, bs="fs",xt=list(bs="tp")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  l_md[[8]] <- gam(logOR ~ 
                     s(area_z,bs = "tp") + s(frag_z,bs = "tp") +
                     s(area_z, SiteCode, bs="fs",xt=list(bs="tp")) +
                     s(SiteCode,bs="re"),
                   data=df, method="REML")
  names(l_md) <- c("~ f(p, k, k|Site)",
                   "~ f(area, frag, frag|Site)","~ f(area, frag, area|Site)","~ f(area, frag, 1|Site)",
                   "~ f(frag, frag|Site)","~ f(area, area|Site)",
                   "~ f(area,frag,frag|Site - ti","~ f(area,frag,area|Site - ti")
  return(l_md) 
}
l_md.logOR_tp_r3removido <- df_md |> 
  f_gam_logOR()
save(l_md.logOR_tp_r3removido,file="dados/Rdata/l_md.logOR_tp_r3removido.Rdata")

#erro: modelo possui mais coeficientes do que dados:
# l_md[[9]] <- gam(logOR ~ 
#                    s(area_z,bs = "tp") + s(frag_z,bs = "tp") + ti(area_z,frag_z) +
#                    s(frag_z, SiteCode, bs="fs",xt=list(bs="tp")) +
#                    s(area_z, SiteCode, bs="fs",xt=list(bs="tp")) +
#                    s(SiteCode,bs="re"),
#                  data=df, method="REML")
