##### Construção da região de buffer da predição a posteriori
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
library(readr)
library(stringr)
library(tidyr)
library(tidyselect)
library(bbmle)
library(DHARMa)
library(mgcv)
library(plyr)
library(dplyr)
library(mgcv)
## funções de ajuste e de plot
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
# objetos comuns
df_logOR <- read_csv(file="dados/csv/df_logOR.csv")
df_md <- df_logOR %>% 
  select(-Uefeito,-(pristine:denominator),-matches("(k_z|^p)")) %>% 
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
#
df_buf <- ddply(df_md,c("contraste","k"),\(dfi){
  v_range <- range(dfi$Uefeito)
  v_quant <- quantile(dfi$Uefeito,probs=c(0.01,0.99))
  filter(dfi,Uefeito%in%v_range) %>% 
    select(-logOR,-c(forest_succession:data_year),-SiteCode) %>% 
    mutate(ext_class=ifelse(Uefeito==v_range[1],"min","max"),
           Uefeito_quant=ifelse(ext_class=="min",v_quant[1],v_quant[2])) %>% 
    relocate(ext_class,.after=contraste) %>% 
    relocate(Uefeito_quant,.after=Uefeito)
}) %>% select(-k)
df_buf %>% 
  pivot_longer(starts_with("Uefeito")) %>% 
  filter(name=="Uefeito") %>% 
  ggplot(aes(x=k_cont,y=value,color=ext_class)) +
  geom_point() +
  geom_line() +
  geom_smooth(method="gam",se=FALSE) +
  scale_color_manual(values=c("darkred","darkgreen")) +
  facet_wrap(~contraste,ncol=3)
#
f_gam <- \(dfi){
  gam(Uefeito~s(k_cont,bs="cr"),data=dfi)
}
l_md <- dlply(df_buf,c("contraste","ext_class"),f_gam)
for(i in names(l_md)){
  p <- gratia::draw(l_md[[i]],residuals = TRUE) + labs(title=i)
  print(p)
  if(!isTRUE(askYesNo("Do you want to see the residuals?"))) break
  p <- gratia::appraise(l_md[[i]])
  print(p)
  if(!isTRUE(askYesNo("next model?"))) break
}
