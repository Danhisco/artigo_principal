## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
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
f_gam <- \(dfi){
  l_md <- list()
  l_md$`s(land)|Site : gi` <- gam(
    logOR ~ 
      s(Uefeito,bs="cr",id="efeito_comum") +
      s(SiteCode,bs="re") + 
      s(Uefeito, by=SiteCode, bs="cr",id="efeito_sitio"),
    data=dfi,method = "REML")
  l_md$`s(land)|Site : gs` <- gam(
    logOR ~ 
      s(Uefeito,bs="cr",id="efeito_comum") +
      s(SiteCode,bs="re") + 
      s(Uefeito, SiteCode, bs = "fs", xt = list(bs = "cr"), id = "efeito_sitio"),
    data=dfi,method = "REML")
  l_md$`s(land) + 1|Site` <- gam(
    logOR ~ 
      s(Uefeito,bs="cr",id="efeito_comum") +
      s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  l_md$`1 + 1|Site` <- gam(
    logOR ~ 1 + s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  return(l_md)
}
## modelos usados
# f_gam <- \(dfi){
#   l_md <- list()
#   l_md$`s(Uefeito)|Site : gi` <- gam(
#     logOR ~
#         forest_succession +
#         s(Uefeito,by=forest_succession,bs="cr",id="effect_forest") +
#         s(data_year,bs="cr") +
#         s(lat,long,bs="tp") +
#         s(SiteCode,bs="re") + 
#         s(Uefeito, by=SiteCode, bs="cr",id="effet_site"),
#       data=dfi,method = "REML")
#   l_md$`s(Uefeito)|Site : gs` <- gam(
#     logOR ~
#       forest_succession +
#       s(Uefeito, by = forest_succession, bs = "cr", id = "effect_forest") +
#       s(data_year, bs = "cr") +
#       s(lat, long, bs = "tp") +
#       s(SiteCode, bs = "re") +
#       s(Uefeito, SiteCode, bs = "fs", xt = list(bs = "cr"), id = "effect_site"),
#     data = dfi, method = "REML")
#   l_md$`1|Site` <- gam(
#     logOR ~
#       forest_succession +
#       s(Uefeito, by = forest_succession, bs = "cr", id = "effect_forest") +
#       s(data_year,bs="cr") +
#       s(lat,long,bs="tp") +
#       s(SiteCode,bs="re"),
#     data=dfi,method = "REML")
#   l_md$`0|Site+1|Site` <- gam(
#     logOR ~
#       forest_succession +
#       s(data_year,bs="cr") +
#       s(lat,long,bs="tp") +
#       s(SiteCode,bs="re"),
#     data=dfi,method = "REML")
#   return(l_md)
# }
l_md_logOR <- dlply(df_md,"contraste",f_gam)
saveRDS(l_md_logOR,file=paste0(v_path,"rds/l_md_simples.rds"))
# l_md_logOR <- readRDS(file=paste0(v_path,"rds/l_md_simples.rds"))
df_tabelaSelecao <- ldply(l_md_logOR,f_TabSelGAMM,.id="pair")
write_csv(df_tabelaSelecao,
          file=paste0(v_path,"rds/tabsel_simples.csv"))



####
saveRDS(l_md_logOR,file=paste0(v_path,"rds/l_md_1aperg_efeito_paisagem_p_sitio.rds"))
## tabela de seleção
# l_md_logOR <- readRDS(file=paste0(v_path,"rds/l_md_1aperg_efeito_paisagem_p_sitio.rds"))
## tabela de selação dos modelos 
df_tabelaSelecao <- ldply(l_md_logOR,f_TabSelGAMM,.id="pair")
write_csv(df_tabelaSelecao,
          file=paste0(v_path,"rds/tabsel_1aperg_efeito_paisagem_p_sitio.csv"))
rm(list="l_md_logOR");gc()
##########################################################
############# Qual o modelo cheio suficiente? ############
##########################################################
# objectos comuns
# modelos mais plausíveis da seleção passada:
# f_SML.gamm <- \(df){md <- l_md_logOR[[df$contraste[1]]][[df$modelo[[1]]]]}
# l_md_1alike <- dlply(df_tabelaSelecao,"contraste",f_SML.gamm)
## modelos usados
## função de atualização do modelo cheio mais plausível
f_gam_selmodcheio <- \(dfi){ # uma vez que todos faoram selecionados para ser
  # formula 
  l_md <- list()
  l_md$`efeito fixo por preserv` <- 
    gam(
      logOR ~
        forest_succession +
        s(Uefeito, by = forest_succession, bs = "cr", id = "effect_forest") +
        s(data_year, bs = "cr") +
        s(lat, long, bs = "tp") +
        s(SiteCode, bs = "re") +
        s(Uefeito, SiteCode, bs = "fs", xt = list(bs = "cr"), id = "effect_site"),
      data = dfi, method = "REML")
  l_md$`efeito fixo comum` <- 
    gam(
      logOR ~
        forest_succession +
        s(Uefeito, bs = "cr", id = "effect_forest") +
        s(data_year, bs = "cr") +
        s(lat, long, bs = "tp") +
        s(SiteCode, bs = "re") +
        s(Uefeito, SiteCode, bs = "fs", xt = list(bs = "cr"), id = "effect_site"),
      data = dfi, method = "REML")
  return(l_md)
}
l_md_logOR <- dlply(df_md,"contraste",f_gam_selmodcheio)
saveRDS(l_md_logOR,file=paste0(v_path,"rds/l_md_2aperg_qual_md_cheio.rds"))
##########################################################
######### Qual o melhor conjunto de covariáveis? #########
##########################################################
# objetos necessários
df_tabelaSelecao <- ldply(l_md_logOR,f_TabSelGAMM,.id="pair")
write_csv(df_tabelaSelecao,
          file=paste0(v_path,"rds/tabsel_2aperg_qual_md_cheio.csv"))
f_SML.gamm <- \(df){md <- l_md_logOR[[df$contraste[1]]][[df$modelo[[1]]]]}
l_md_1alike <- dlply(df_tabelaSelecao,"contraste",f_SML.gamm)
f_gam_update <- \(mdname,
                  name_path = paste0(v_path,"rds/l_md_3aperg_quais_cov_")){
  mdgam <- l_md_1alike[[mdname]]
  f_ref <- formula(mdgam)
  original_data <- mdgam$model
  f_gam <- \(fi){
    gam(formula=fi,data=original_data,method = "REML")
  }
  l_md <- list()
  l_md$cheio <- f_ref
  l_md$`- ano` <- update.formula(f_ref, . ~ . - s(data_year, bs = "cr"))
  l_md$`- coord` <- update.formula(f_ref, . ~ . - s(lat, long, bs = "tp"))
  if(!grepl("by = forest_succession",as.character(f_ref)[3])){
    l_md$`- pert` <- update.formula(f_ref, . ~ . - forest_succession)
    l_md$`- pert e ano` <- update.formula(l_md$`- ano`, . ~ . - forest_succession)
    l_md$`- pert e coord` <- update.formula(l_md$`- coord`, . ~ . - forest_succession)
    l_md$`- ano e coord` <- update.formula(l_md$`- coord`, . ~ . - s(data_year, bs = "cr"))
    l_md$`- cov` <- update.formula(l_md$`- ano e coord`, . ~ . - forest_succession)
  }else{
    l_md$`- cov` <- update.formula(l_md$`- coord`, . ~ . - s(data_year, bs = "cr"))
  }
  l_md <- lapply(l_md,f_gam)
  mdname <- gsub("Área per se","areaperse",mdname) %>% 
    gsub("Frag. per se","fragperse",.) %>% 
    gsub("Frag. total","fragtotal",.)
  saveRDS(l_md,
          file=paste0(name_path,mdname,".rds"))
  rm(l_md);gc()
}
lapply(names(l_md_1alike),f_gam_update)
## tabelas de seleção
paths <- list.files(path=gsub("/l_md_3aperg_quais_cov_","",
                              eval(formals(f_gam_update)$name_path)),
                    pattern="3aperg",full.names = T)
df_tabelaSelecao <- lapply(paths,\(li){
  l_md <- readRDS(li)
  df_return <- f_TabSelGAMM(l_md)
  mutate(df_return,
         contraste=str_extract(li,"(?<=cov\\_)(.*?)(?=\\.rds)") %>% 
           gsub("areaperse","Área per se",.) %>% 
           gsub("fragperse","Frag. per se",.) %>% 
           gsub("fragtotal","Frag. total",.))
}) %>% do.call("rbind",.)
write_csv(df_tabelaSelecao,
          file=paste0(v_path,"rds/tabsel_3aperg_quais_cov.csv"))
############################################
#
#
############################################
###  tabelas e diagnósticos
l_md <- readRDS(paste0(v_path,"rds/l_md_2aperg_qual_md_cheio.rds"))
l_md <- lapply(split(df_tabsel,df_tabsel$contraste),\(dfi){
  with(dfi,{l_md[[contraste]][[modelo]]})
})
library(magick)
f_diagplots <- \(lmd){
  lapply(names(lmd),\(li){
    md <- lmd[[li]]
    # lista de objetos 
    ldiag <- list()
    ldiag$ft_md <- as_flextable(md) %>% 
      bg(., bg = "white", part = "all")
    ldiag$appraise <- gratia::appraise(md)
    ldiag$draw <- gratia::draw(md)
    # escreve as figuras como png e salva os caminhos
    lpaths <- lapply(names(ldiag),\(i){
      vpath <- paste0(v_path,"figuras/",i,"_lapply.png")
      if(i=="ft_md"){
        save_as_image(ldiag[[i]],path=vpath)
      }else{
        ggsave(plot = ldiag[[i]],width = 8,height = 6,filename = vpath)
      }})
    names(lpaths) <- names(ldiag)
    # le como image object
    l_png <- lapply(lpaths,\(li) image_read(li) %>% image_trim)
    # ajuste de tamanhos
    base_image <- lapply(l_png[-1],\(li){
      image_resize(li,
                   geometry_size_pixels(width = image_info(l_png[["ft_md"]])$width/2,
                                        height = image_info(l_png[["ft_md"]])$height))
    })
    base_image <- image_append(do.call("c",base_image),stack=FALSE)
    base_image <- image_resize(
      base_image, 
      geometry_size_pixels(width = image_info(l_png[["ft_md"]])$width, 
                           height = image_info(l_png[["ft_md"]])$height))
    # figura final salvamento e limpeza
    final_image <- image_append(c(l_png[[1]],base_image),stack = TRUE)
    v_name <- gsub("Área per se","areaperse",li) %>% 
      gsub("Frag. per se","fragperse",.) %>% 
      gsub("Frag. total","fratotal",.)
    image_write(final_image,paste0(v_path,"figuras/diag_",v_name,".png"), format = "png")
    file.remove(do.call("c",lpaths))
  })
}
