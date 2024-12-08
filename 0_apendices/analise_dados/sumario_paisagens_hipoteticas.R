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
if(!file.exists("1_to_compile_dissertacao_EM_USO/00_Resultados/tabelas/tabselecao_sumario_paisagens.csv")){
  df_tabsel <- f_TabSelGAMM(l_md)
  write_csv(df_tabsel,"1_to_compile_dissertacao_EM_USO/00_Resultados/tabelas/tabselecao_sumario_paisagens.csv")  
}else{
  df_tabsel <- read_csv("1_to_compile_dissertacao_EM_USO/00_Resultados/tabelas/tabselecao_sumario_paisagens.csv")
}
######################################### ESTIMATIVAS
md_sumario <- l_md[[df_tabsel$modelo[1]]]
new_data <- expand.grid(
  k = unique(df_nSAD$k),
  land = unique(df_nSAD$land),
  SiteCode = "SPigua1"  # Replace with a fixed value
)
predictions <- predict(md_sumario, 
                       newdata = new_data, 
                       type = "link", 
                       se.fit = TRUE, 
                       exclude = "s(SiteCode)")
new_data$fit <- predictions$fit
new_data$se <- predictions$se.fit
new_data <- new_data %>% mutate(
  lower = fit - 1.96 * se,
  upper = fit + 1.96 * se
)
new_data$k_factor <- as.factor(new_data$k)
f_geom_final <- \(vs){
  list(
    scale_color_manual(values = c(
      "Contemporânea"="darkred",
      "Aglomerada"="darkblue",
      "Prístina"="darkgreen")),
    scale_y_continuous(expand=c(0,0)),
    facet_wrap(~label),
    labs(
      x = "k", 
      y = "logito(Proporção de SADs simuladas com boa congruência)", 
      color = "Paisagem\nhipotética"
    ),
    theme_classic(),
    theme(text=element_text(size=vs))
  )
} 
f_efeitosumario <- \(dfi,
                     vsize=5,
                     vrange,
                     vposleg=c(0.49,0.9)){
  dfi %>% 
    mutate(land = gsub("contemp","Contemporânea",land) %>% 
             gsub("ideal","Prístina",.) %>% 
             gsub("non_frag","Aglomerada",.),
           land = factor(land,levels=c("Contemporânea",
                                       "Aglomerada",
                                       "Prístina")),
           label="Estimativa sem o efeito particular da paisagem") %>% 
    ggplot(aes(x = k_factor, y = fit, color = land, group = land)) +
    geom_line(aes(x=k)) +
    geom_point(size = 3, position = position_dodge(0.5)) + # Points for categories
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.2,
      position = position_dodge(0.5)
    ) + # Error bars for uncertainty
    f_geom_final(vsize) +
    scale_y_continuous(expand=c(0,0),limits = vrange)+
    theme(legend.position = "inside",
          legend.position.inside = vposleg,
          legend.direction="horizontal")
}
#
f_obssumario <- \(dfi,
                  vsize=5){
  dfi %>% 
    mutate(k_factor=factor(k),
           label="Distribuição dos valores observados de logito",
           land = gsub("contemp","Contemporânea",land) %>% 
             gsub("ideal","Prístina",.) %>% 
             gsub("non_frag","Aglomerada",.),
           land = factor(land,levels=c("Contemporânea",
                                       "Aglomerada",
                                       "Prístina"))) %>% 
    ggplot(aes(x = k_factor, y = p, color = land, group = interaction(k_factor,land))) +
    geom_jitter() +
    geom_boxplot() +
    f_geom_final(vsize) +
    theme(legend.position = "none")
}
ve=0.9
df_nSAD <- df_nSAD %>% 
  mutate(p = (nSAD+ve)/(100+2*ve),
         p = log(p/(1-p)))
l_p <- list()
l_p$obs <- f_obssumario(df_nSAD,vsize = 15)
l_p$efeito <- f_efeitosumario(new_data,
                              vsize = 15,
                              vposleg = c(0.5,0.1),
                              vrange = range(df_nSAD$p))
p <- arrangeGrob(grobs=l_p,ncol=2)
ggsave(filename = "figuras/sumario_paisagens.png",plot = p,
       width = 16,height = 8)
img <- image_read("figuras/sumario_paisagens.png") %>% 
  image_trim() %>% 
  image_resize("50%")
image_write(img,"figuras/sumario_paisagens.jpeg")