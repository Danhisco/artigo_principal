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
  select(SiteCode,lat,long, forest_succession) %>% 
  filter(forest_succession!="capoeira")
#
df_ad <- read_csv(file="dados/csv_SoE/df_congruencia_simulacao.csv") %>% 
  rename(nSAD = nCongKS,
         land = land_type) %>% 
  select(-c(Smed:Smax)) %>% 
  inner_join(df_coord) %>% 
  mutate(k=round(k,2),
         across(c(land,SiteCode,forest_succession),factor))
df_ad$SiteCode <- factor(df_ad$SiteCode)
#saveRDS(df_ad,file="0_apendices/analise_dados/rar/df_ad_sumario.rds")
#############
# df_adSAD <- readRDS("5_resultados/df_adSAD.rds")
# df_nSAD <- adply(unique(df_adSAD$land_type),1,\(i){
#   df_return <- df_adSAD %>% filter(taxaU==i & land_type==i) %>% 
#     select(-land_type,-c(Smed:Smax)) %>% 
#     rename("land"="taxaU",
#            "nSAD"="nCongKS")
# },.id=NULL) %>% 
#   mutate(k=round(k,2),
#          across(c(land,SiteCode,k),factor)) %>% 
#   inner_join(df_coord)
# df_nSAD$SiteCode <- factor(df_nSAD$SiteCode)
##############
f_gam <- \(vf,dfi){
  gam(formula=vf,
      family='binomial',
      data=dfi,
      method="REML")
}
l_f <- list()
# modelo cheio
l_f$`s(k,by=land + class_pert) + (lat,long)` <- 
  cbind(nSAD,100-nSAD) ~ 
  s(k,by=interaction(land,forest_succession),bs="cr",id="fixo") +
  s(lat,long) + 
  te(k,SiteCode,bs=c("cr","re"),by=land,id="random")
# modelo sem classe de perturbação
l_f$`s(k,by=land) + (lat,long)` <- 
  cbind(nSAD,100-nSAD) ~ 
  s(k,by=land,bs="cr",id="fixo") +
  s(lat,long) + 
  te(k,SiteCode,bs=c("cr","re"),by=land,id="random")
# modelo sem coordenadas
l_f$`s(k,by=land * class_pert)` <- 
  cbind(nSAD,100-nSAD) ~ 
  s(k,by=interaction(land,forest_succession),bs="cr",id="fixo") +
  te(k,SiteCode,bs=c("cr","re"),by=land,id="random")
# modelo sem cov
l_f$`s(k,by=land)` <- 
  cbind(nSAD,100-nSAD) ~ 
  s(k,by=land,bs="cr",id="fixo") +
  te(k,SiteCode,bs=c("cr","re"),by=land,id="random")
#




# 
# 
# 
# 
# l_f$`land * k + class_pert` <- 
#   cbind(nSAD,100-nSAD) ~ land * k + forest_succession + s(SiteCode,by=land,bs="re")
# l_f$`land * k` <- 
#   cbind(nSAD,100-nSAD) ~ land * k + s(SiteCode,by=land,bs="re")
# l_f$`land * k + s(lat,long) + land|Site` <- cbind(nSAD,100-nSAD) ~ land * k + s(lat,long) + s(SiteCode,by=land,bs="re")
# l_f$`land + s(lat,long) + k + s(lat,long) + land|Site` <- cbind(nSAD,100-nSAD) ~ land + s(lat,long) + k + s(lat,long) + s(SiteCode,by=land,bs="re")
# l_f$`land + s(lat,long) + land|Site` <- cbind(nSAD,100-nSAD) ~ land + s(lat,long) + s(SiteCode,by=land,bs="re")
# l_f$`k + s(lat,long) + 1|Site` <- cbind(nSAD,100-nSAD) ~ k + s(lat,long) + s(SiteCode,bs="re") 
# l_f$`1 + s(lat,long) + 1|Site` <- cbind(nSAD,100-nSAD) ~ 1 + s(lat,long) + s(SiteCode,bs="re")
if(!file.exists("1_to_compile_dissertacao_EM_USO/00_Resultados/tabelas/tabselecao_sumario_paisagens.csv")){
  # doMC::registerDoMC(2)
  l_md <- llply(l_f,f_gam,dfi=df_ad,.parallel = FALSE)
  saveRDS(l_md,file="dados/csv_SoE/Rdata/l_md_sumario")
  l_tabsel <- f_TabSelGAMM(l_md,test_moranK = TRUE)
  write_csv(l_tabsel$tabsel,"1_to_compile_dissertacao_EM_USO/00_Resultados/tabelas/tabselecao_sumario_paisagens.csv")  
  saveRDS(l_tabsel$l_moranK,file="1_to_compile_dissertacao_EM_USO/00_Resultados/rds/l_moranK_sumario_paisagens.rds")
}else{
  df_tabsel <- read_csv("1_to_compile_dissertacao_EM_USO/00_Resultados/tabelas/tabselecao_sumario_paisagens.csv")
  l_md <- readRDS(file="dados/csv_SoE/Rdata/l_md_sumario")
}
##########
l_f <- list()
l_f$`s(k,by=land) + land|Site` <- 
  cbind(nSAD,100-nSAD) ~ 
  s(k,by=land,bs="cr",id="fixo") +
  s(SiteCode,bs="re",by=land)
l_f$`land + land|Site` <- 
  cbind(nSAD,100-nSAD) ~ 
  land +
  s(SiteCode,bs="re",by=land)
doMC::registerDoMC(2)
l_md_minimo <- llply(l_f,f_gam,dfi=df_ad,.parallel = TRUE)
l_md_minimo$`s(k,by=land)` <- l_md$`s(k,by=land)`
l_md_final <- c(l_md,l_md_minimo)
df_tabsel_minimo <- f_TabSelGAMM(l_md_final)
write_csv(df_tabsel_minimo,"1_to_compile_dissertacao_EM_USO/00_Resultados/tabelas/tabselecao_sumario_minimo_paisagens.csv")  

######################################### ESTIMATIVAS
md_sumario <- l_md[[df_tabsel$modelo[1]]]
new_data <- expand.grid(
  k = unique(df_ad$k),
  land = unique(df_ad$land),
  forest_succession = unique(df_ad$forest_succession),
  lat = 0,
  long = 0,
  SiteCode = "SPigua1"  # Replace with a fixed value
)
#
l_dfpred <- list()
l_dfpred$fixo_e_aleat <- md_sumario$model
l_dfpred$fixo_e_aleat[,1] <- l_dfpred$fixo_e_aleat[,1][,1]
names(l_dfpred$fixo_e_aleat)[c(1,3)] <- c("nSAD","forest_succession")
lpred <- predict(md_sumario,
                 type = "response", 
                 se.fit = TRUE)
l_dfpred$fixo_e_aleat$fit <- lpred$fit
l_dfpred$fixo_e_aleat$se <- lpred$se.fit
l_dfpred$fixo_e_aleat <- mutate(
  l_dfpred$fixo_e_aleat,
  lower = fit - 1.96 * se,
  upper = fit + 1.96 * se,
  forest_succession = gsub("cont.","",forest_succession) %>% 
    gsub("ideal.","",.) %>% gsub("non_frag.","",.)
)
#
l_dfpred$fixo <- expand.grid(
  k = unique(df_ad$k),
  land = unique(df_ad$land),
  forest_succession = unique(df_ad$forest_succession),
  lat = 0,
  long = 0,
  SiteCode = "SPigua1"
)
lpred <-  predict(md_sumario,
                  newdata=l_dfpred$fixo,
                  type = "response",
                  exclude=c("s(lat,long)",
                            "te(k,SiteCode):landcont",
                            "te(k,SiteCode):landideal",
                            "te(k,SiteCode):landnon_frag"),
                  se.fit = TRUE)
l_dfpred$fixo$fit <- lpred$fit
l_dfpred$fixo$se <- lpred$se.fit
l_dfpred$fixo <- mutate(
  l_dfpred$fixo,
  lower = fit - 1.96 * se,
  upper = fit + 1.96 * se
)
#######
l_dfpred <- lapply(l_dfpred,\(dfi){
  mutate(dfi,across(fit:upper,~.x*100),
         land=factor(land,
                     levels=c("cont",
                              "non_frag",
                              "ideal"),
                     labels=c("fragmentada",
                              "aglomerada",
                              "prístina")),
         pert_class = factor(forest_succession,
                             levels=c("primary",
                                      "primary/secondary",
                                      "secondary"),
                             labels=c("baixa",
                                      "mediana",
                                      "alta"))
         ) %>% 
    select(-forest_succession,-lat,-long) %>% 
    relocate(pert_class,.after="k")
})
library(hexbin)
p <- l_dfpred$fixo_e_aleat %>% 
  ggplot(aes(x=k,y=nSAD)) +
  # geom_point(alpha=0.3) +
  geom_line(aes(y=fit,group = SiteCode),alpha=0.3) +
  # geom_boxplot(aes(group=k),alpha=0.6) +
  geom_ribbon(data=l_dfpred$fixo,
              aes(x=k,y=fit,ymin=lower,ymax=upper),
              color="blue",
              fill="lightblue") +
  geom_line(data=l_dfpred$fixo,
            aes(x=k,y=fit),
            color="darkred") +
  geom_hex(bins = 50,alpha=0.5) + 
  scale_fill_gradient("contagem",low = "gray", high = "black", na.value = NA) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x="Proporção de propágulos na vizinhança imediata (k)",
       y="número de SAD simuladas com boa congruência com o observado (nSAD)",
       title="Descrição da congruência absoluta da SAD simulada nas paisagens hipotéticas",
       subtitle="HGAM: nSAD ~ s(k,by=interaction(paisagem hipotética, classe de perturbação) + s(lat,long) + te(k,Sítio,by=paisagem hipotética)") +
  facet_grid(pert_class ~ land) +
  theme(strip.text = element_text(size=12, #margin=margin(),
                                  face="bold"))
ggsave(filename="figuras/descricao_congruencia_relativa.png",
       p,height = 7.33,width = 13.8)
library(magick)
img <- image_read("figuras/descricao_congruencia_relativa.png") %>% 
  image_resize("50%") %>% 
  image_trim()
image_write(img,path = "figuras/descricao_congruencia_relativa.png")
############################################################
# f_geom_final <- \(vs){
#   list(
#     scale_color_manual(values = c(
#       "Contemporânea"="darkred",
#       "Aglomerada"="darkblue",
#       "Prístina"="darkgreen")),
#     scale_y_continuous(expand=c(0,0)),
#     facet_wrap(~label),
#     labs(
#       x = "k", 
#       y = "logito(Proporção de SADs simuladas com boa congruência)", 
#       color = "Paisagem\nhipotética"
#     ),
#     theme_classic(),
#     theme(text=element_text(size=vs))
#   )
# } 
# f_efeitosumario <- \(dfi,
#                      vsize=5,
#                      vrange,
#                      vposleg=c(0.49,0.9)){
#   dfi %>% 
#     mutate(land = gsub("contemp","Contemporânea",land) %>% 
#              gsub("ideal","Prístina",.) %>% 
#              gsub("non_frag","Aglomerada",.),
#            land = factor(land,levels=c("Contemporânea",
#                                        "Aglomerada",
#                                        "Prístina")),
#            label="Estimativa sem o efeito particular da paisagem") %>% 
#     ggplot(aes(x = k_factor, y = fit, color = land, group = land)) +
#     geom_line(aes(x=k)) +
#     geom_point(size = 3, position = position_dodge(0.5)) + # Points for categories
#     geom_errorbar(
#       aes(ymin = lower, ymax = upper),
#       width = 0.2,
#       position = position_dodge(0.5)
#     ) + # Error bars for uncertainty
#     f_geom_final(vsize) +
#     scale_y_continuous(expand=c(0,0),limits = vrange)+
#     theme(legend.position = "inside",
#           legend.position.inside = vposleg,
#           legend.direction="horizontal")
# }
# #
# f_obssumario <- \(dfi,
#                   vsize=5){
#   dfi %>% 
#     mutate(k_factor=factor(k),
#            label="Distribuição dos valores observados de logito",
#            land = gsub("contemp","Contemporânea",land) %>% 
#              gsub("ideal","Prístina",.) %>% 
#              gsub("non_frag","Aglomerada",.),
#            land = factor(land,levels=c("Contemporânea",
#                                        "Aglomerada",
#                                        "Prístina"))) %>% 
#     ggplot(aes(x = k_factor, y = p, color = land, group = interaction(k_factor,land))) +
#     geom_jitter() +
#     geom_boxplot() +
#     f_geom_final(vsize) +
#     theme(legend.position = "none")
# }
# ve=0.9
# df_nSAD <- df_nSAD %>% 
#   mutate(p = (nSAD+ve)/(100+2*ve),
#          p = log(p/(1-p)))
# l_p <- list()
# l_p$obs <- f_obssumario(df_nSAD,vsize = 15)
# l_p$efeito <- f_efeitosumario(new_data,
#                               vsize = 15,
#                               vposleg = c(0.5,0.1),
#                               vrange = range(df_nSAD$p))
# p <- arrangeGrob(grobs=l_p,ncol=2)
# ggsave(filename = "figuras/sumario_paisagens.png",plot = p,
#        width = 16,height = 8)
# img <- image_read("figuras/sumario_paisagens.png") %>% 
#   image_trim() %>% 
#   image_resize("50%")
# image_write(img,"figuras/sumario_paisagens.jpeg")
# 
# ### diagnóstico
# source("source/GAMMtools.R")
# v_path <- "dados/csv_SoE/figuras/"
# f_diag(hgam = md_sumario,vname = "sumario",patsave = "hgam")
# 
# 
# 
# ### logitos empíricos
# df_mdsumar <- md_sumario$model
# df_mdsumar[,1] <- df_mdsumar[,1][,1]
# names(df_mdsumar)[1] <-  "nSAD"
# df_mdsumar$logito <- predict(md_sumario,se=FALSE)
# df_save <- df_mdsumar %>% 
#   filter(nSAD %in% c(0,100)) %>% 
#   group_by(nSAD,k,land) %>% 
#   summarise(logito_avg = mean(logito),
#             logito_sd = sd(logito))
# saveRDS(df_save,file="./dados/csv_SoE/df_mdsumar.rds")
# df_md <- left_join(df_mdsumar,
#                    df_save,
#                    by=c("nSAD","k","land")) %>% 
#   rename("forest_succession"=3) %>% 
#   mutate(
#     logito = ifelse(
#       is.na(logito_avg),
#       log(nSAD/(100-nSAD)),
#       logito_avg
#       ),
#     forest_succession = gsub("cont.|ideal.|non_frag.",
#                              "",forest_succession)
#     ) %>% 
#   select(-c(logito_avg,logito_sd))
# saveRDS(df_md,file="./dados/csv_SoE/df_mdsumar.rds")