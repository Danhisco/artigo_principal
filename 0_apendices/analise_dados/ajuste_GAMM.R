## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
library(mgcv)
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
###########################################################################
###### Como expressar o efeito da paisagem por sítio de amostragem?  ######
###########################################################################
# objetos comuns
df_logOR <- read_csv(file="dados/csv/df_logOR.csv")
df_md <- df_logOR %>% select(-Uefeito,-(pristine:denominator),-matches("(k_z|^p)")) %>% 
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
  relocate(Uefeito,.after="logOR")
## modelos usados (versão exploração apenas da estrutura hierarquica)
## Seguindo as especifícações de Pedersen et al. 2019 https://peerj.com/articles/6876/
f_gam <- \(dfi){
  l_md <- list()
  l_md$`s(land)|Site : gi` <- gam(
    logOR ~ 
      s(Uefeito,bs="cr",m=2, id="efeito_comum") +
      s(SiteCode,bs="re") + 
      s(Uefeito, by=SiteCode, bs="cr",m=1, id="efeito_sitio"),
    data=dfi,method = "REML")
  l_md$`s(land)|Site : gs` <- gam(
    logOR ~ 
      s(Uefeito,bs="cr",m=2, id="efeito_comum") +
      s(Uefeito, SiteCode, bs = "fs", xt=list(bs = "cr"), m=2, id="efeito_sitio"),
    data=dfi,method = "REML")
  l_md$`s(land) + 1|Site` <- gam(
    logOR ~ 
      s(Uefeito,bs="cr",m=2,id="efeito_comum") +
      s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  l_md$`1 + 1|Site` <- gam(
    logOR ~ 1 + s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  return(l_md)
}
l_md_logOR <- dlply(df_md,"contraste",f_gam)
saveRDS(l_md_logOR,file=paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))0
# l_md_logOR <- readRDS(file=paste0(v_path,"rds/l_md_simples.rds"))
df_tabelaSelecao <- ldply(l_md_logOR,f_TabSelGAMM,.id="pair")
write_csv(df_tabelaSelecao,
          file=paste0(v_path,"rds/tabsel_simples.csv"))
###  tabelas e diagnósticos
l_md <- readRDS(paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))
df_tabsel <- read_csv(paste0(v_path,"rds/tabsel_simples.csv")) %>%
  filter(dAICc==0)
l_md <- lapply(split(df_tabsel,df_tabsel$contraste),\(dfi){
  with(dfi,{l_md[[contraste]][[modelo]]})
})
vpaths <- f_diagplots(l_md)
