# pacotes
library(gratia)
library(doMC)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
library(readr)
library(stringr)
library(tidyr)
library(bbmle)
library(DHARMa)
# library(lme4)
library(mgcv)
library(plyr)
library(dplyr)
## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
library(mgcv)
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
# objetos comuns
df_md <- readRDS(file="dados/csv_SoE/df_logOR.rds") %>% 
  mutate(k_cont = as.numeric(as.character(k)),
         across(c(SiteCode,contraste,forest_succession),factor)) %>% 
  relocate(k,.after="Uefeito")
saveRDS(df_md,
        file=paste0(v_path,"rds/df_md.rds"))
# df_p_extensoes <- read_csv("dados/csv/df_p_extensoes.csv")
# v_filter <- with(df_p_extensoes,{
#   v_i <- sapply(c(4,2,1),
#                 \(x) which.min(abs(lado_km - x)))
#   return(lado_km[v_i])
#   }
# ) 
# df_p <- filter(df_p_extensoes,lado_km %in% v_filter) %>% 
#   mutate(lado_km = round(lado_km,0))
######################################################################################
#################### descrição de logU/U
df_logUU_pk <- readRDS(file="dados/csv_SoE/rds/df_logUUpk.rds") %>% 
  mutate(SiteCode=factor(SiteCode))
library(mgcv)
f_gam <- \(dfi){
  #
  l_md <- list()
  l_md[[1]] <- gam(
    Uefeito ~ 
      te(p,k,
         bs=c("cr","cr"),m=2,
         id = "fixo") +
      s(k, SiteCode, 
        bs = "fs", xt=list(bs = "cr"), m=2, 
        id="efeito_sitio"),
    data=dfi, method = "REML")
  l_md[[2]] <- gam(
    Uefeito ~ 
      s(k,bs="cr",m=2,id="fixo_k") +
      s(p,bs="cr",m=2,id="fixo_p") +
      s(k, SiteCode, 
        bs = "fs", xt=list(bs = "cr"), m=2, 
        id="efeito_sitio"),
    data=dfi, method = "REML")
  l_md[[3]] <- gam(
    Uefeito ~ 
      s(k,bs="cr",m=2,id="fixo") +
      s(k, SiteCode, 
        bs = "fs", xt=list(bs = "cr"), m=2, 
        id="efeito_sitio"),
    data=dfi, method = "REML")
  #
  names(l_md) <- c("~ te(p,k)",
                   "~ s(k) + s(p)",
                   "~ s(k)")
  return(l_md)
}
if(FALSE){
  doMC::registerDoMC(2)
  l_md <- dlply(df_logUU_pk,"contraste",f_gam,.parallel = FALSE)
  saveRDS(l_md,file="dados/csv_SoE/rds/l_md_logUUpk2.rds")
}else{
  l_md_logUUpk <- readRDS(file="dados/csv_SoE/rds/l_md_logUUpk.rds")
}
# df2ef <- df_logUU_pk %>% filter(contraste!="Frag. total")
f_gam <- \(df2ef){
  l_md <- list()
  l_md[[1]] <- gam(Uefeito ~ 
                     te(p,k,
                        by=contraste,
                        bs=c("cr","cr"),m=2,
                        id = "fixo") +
                     s(k, SiteCode, 
                       by=contraste,
                       bs = "fs", xt=list(bs = "cr"), m=2, 
                       id="efeito_sitio"),
                   data=df2ef, method = "REML")
  l_md[[2]] <- gam(Uefeito ~ 
                     te(p,k,
                        bs=c("cr","cr"),m=2,
                        id = "fixo") +
                     s(k, SiteCode, 
                       bs = "fs", xt=list(bs = "cr"), m=2, 
                       id="efeito_sitio"),
                   data=df2ef, method = "REML")
  names(l_md) <- c("por contraste","comum")
  return(l_md)
}
if(FALSE){
  l_md_2ef <- f_gam(
    df_logUU_pk %>% filter(contraste!="Frag. total")
  )
  saveRDS(l_md_2ef,file="dados/csv_SoE/rds/l_md_logUUpk_2ef.rds")
}

######################## diagnósticos 
## Os modelos ajustados para os efeitos individuais
# tabela de seleção
if(FALSE){
  l_tabsel <- lapply(l_md_logUUpk,f_TabSelGAMM,test_moranK=FALSE)
  df_tabsel <- lapply(names(l_tabsel),\(li){
    mutate(l_tabsel[[li]],efeito=li)}) %>% 
    do.call("rbind",.) %>% 
    relocate(efeito) %>% 
    select(efeito:dev.expl) %>% 
    mutate(across(where(is.numeric),~round(.x,2)))
  saveRDS(df_tabsel,file="dados/csv_SoE/rds/df_tabsel_logUU_pk.rds")
}else{
  df_tabsel <- readRDS(file="dados/csv_SoE/rds/df_tabsel_logUU_pk.rds")
}
# modelos mais plausíveis
df_tabsel <- lapply(names(l_tabsel),\(li){
  mutate(l_tabsel[[li]],efeito=li) %>% 
    select(-pvalue,-MoranI_stat_res) %>% 
    filter(dAICc==0)
}) %>% do.call("rbind",.)
# lista com os modelos mais plausíveis
l_md <- dlply(df_tabsel,"efeito",\(dfi){
  l_md_logUUpk[[dfi$efeito]][[dfi$modelo]]
})
# diagnostico dos mais plausíveis
if(FALSE){
  l_paths <- lapply(names(l_md),\(li){
    f_diag(hgam=l_md[[li]],v_path = v_path,vname = li,patsave = "te_pk")
  })
}
########################################### predito e observado
#li <- names(l_md)[[1]]
#hgam <- l_md[[li]]
l_dfpred <- readRDS(file = "dados/csv_SoE/rds/l_dfpred_md_cong_absoluta.rds")
vsites105 <- l_dfpred$fixo_e_aleat$SiteCode %>% unique
vSites <- read_csv(file = "dados/df_dados_disponiveis.csv") %>% 
  filter(forest_succession!="capoeira") %>% 
  pull(SiteCode) %>% unique
vSites <- intersect(vSites,vsites105)
# criação de df geral
dfpred100 <- ldply(l_md,\(li){
  li$model %>% 
    filter(p==1) %>% 
    summarise(maxU = max(Uefeito), minU = min(Uefeito))
}) %>% summarise(maxU = max(maxU), minU = min(minU))
## fixo e aleatorio: variações por sítio

df_obs_pred_plot <- lapply(names(l_md),\(li){
  hgam <- l_md[[li]]
  dfobs <- hgam$model %>% 
    filter(k>=0.49999,SiteCode %in% vSites) %>% 
    mutate(efeito=li) %>% relocate(efeito)
  dfpred <- as.data.frame(predict.gam(hgam,newdata = dfobs,se.fit = TRUE))
  dfpred <- cbind(dfobs,dfpred)
})
saveRDS(df_obs_pred_plot,"dados/csv_SoE/rds/df_obs_pred_plot_logUU_pk.rds")


f_obs_predito_bysite <- \(lmd){
  
  
  
  dfpred <- dfpred %>% 
    mutate(
      p_class = case_when(
        p==1 ~ "%CF = 100",
        p<1 & p>=0.80 ~ "80 ≤ %CF < 100",
        p<0.80 & p>=0.60 ~ "60 ≤ %CF < 80",
        p<0.60 & p>=0.30 ~ "30 ≤ %CF < 60",
        p<0.30 ~ "%CF < 30"),
      p_class = factor(p_class,levels=c("%CF < 30",
                                        "30 ≤ %CF < 60",
                                        "60 ≤ %CF < 80",
                                        "80 ≤ %CF < 100",
                                        "%CF = 100")),
      lower = fit - 1.96 * se.fit,
      upper = fit + 1.96 * se.fit
      )
  df_ij <- dfpred %>% 
    select(SiteCode,p_class) %>% 
    distinct() %>% 
    group_by(p_class) %>% 
    tally() %>% 
    mutate(label = paste0(p_class,", n Sítios=",n)) %>% 
    select(-n)
  df_ij$label <- factor(df_ij$label,
                        levels=unique(df_ij$label))
  dfp <- inner_join(dfpred,df_ij)
  # dfp <- filter(dfp100,p_class!="%CF = 100")
  # dfp100 <- filter(dfp100,p_class=="%CF = 100")
  dfp %>% 
    ggplot(aes(x=k,y=Uefeito)) +
      geom_boxplot(aes(group=k),alpha=0.6) +
      geom_hline(yintercept = 0,color="black") +
      geom_ribbon(aes(x=k,y=fit,ymin=lower,ymax=upper,group = SiteCode),
                  fill="lightblue",
                  alpha=0.4) + 
      geom_line(aes(y=fit,group = SiteCode),color="darkred",alpha=0.5) +
      geom_point(aes(color=p),alpha=0.7) +
      scale_colour_gradient2("% CF",midpoint=0.5,
                             low="red",
                             mid = "yellow",
                             high = "darkgreen") +
      labs(x="k",y="logU/U",title=vname) +
      facet_wrap(~label,ncol=3) +
      theme_classic() +
      theme(plot.margin=unit(c(0,0.2,0,0), "cm"),
            legend.position = "inside",
            legend.position.inside = c(0.49,0.92),
            legend.direction="horizontal") 
}
l_p_fixo_e_aleat <- lapply(names(l_md),\(li) f_obs_predito_bysite(hgam=l_md[[li]],vname=li))
names(l_p_fixo_e_aleat) <- names(l_md)
saveRDS(l_p_fixo_e_aleat,file="./figuras/logUU_construcao/l_p_fixo_e_aleat.rds")

## fixo:o efeito médio desconsiderando a variabilidade entre sítios
f_calcPI <- \(hgam,vname,
              ctextsize=5){
  dfobs <- hgam$model %>% filter(k>=0.49999,SiteCode %in% vSites)
  dfrange <- data.frame(
    X = c("p","k"),
    max = sapply(dfobs[,c("p","k")],max),
    min = sapply(dfobs[,c("p","k")],min)
  )
  dfpred0 <- expand.grid(
    p = seq(dfrange[dfrange$X=="p","min"],
            dfrange[dfrange$X=="p","max"],
            length.out=100),
    k = seq(dfrange[dfrange$X=="k","min"],
            dfrange[dfrange$X=="k","max"],
            length.out=100)
  ) %>% mutate(SiteCode=factor("BAjuss"))
  dfpred <- cbind(
    dfpred0,
    as.data.frame(
      predict.gam(hgam,newdata = dfpred0,se.fit=TRUE,
                  exclude = c("s(k,SiteCode)"))
      )
    ) %>% 
    mutate(efeito = vname)
  return(dfpred)
}
dfpred <- lapply(names(l_md),\(li) f_calcPI(hgam=l_md[[li]],vname=li)) %>% 
  do.call("rbind",.)
saveRDS(dfpred,file="dados/csv_SoE/rds/dfpred_efeitofixo_logUUpk.rds")
vrangelogUU <- range(dfpred$fit)
library(metR)
f_ggplot <- \(dfp,
              xlab="p",
              ylab="k",
              zlab="fit",
              scales_free=TRUE,
              lw=1,ctextsize=5,stripspace=0.5,textsize=15){
  f_perfit <- 
  
  
  if(scales_free){
    vbreaks <- quantile(dfp[[zlab]],c(0.05,0.25,0.50,0.75,0.95)) %>%
      round(.,digits=3)
  }else{
    vbreaks <- c(-0.10,0,0.10,0.20,0.30)
  }
  p <- dfp %>% 
    pivot_longer(c(fit,se.fit)) %>% 
    ggplot(aes(x=.data[[xlab]],y=.data[[ylab]],z=value)) +
    geom_tile(aes(fill=value),width = 0.0090, height = 0.0090) +
    geom_contour(color = "black",linewidth=lw,
                 breaks=vbreaks) +
    geom_text_contour(size=ctextsize,
                      breaks=vbreaks,
                      color="black",
                      check_overlap = TRUE, 
                      stroke=0.1) +
    scale_fill_gradient2("logU/U",
                         midpoint=0,
                         low="#440154",
                         mid = "#21908CFF",
                         high = "#FDE725") +
    labs(y=ifelse(ylab=="p","%CF","grau de limitação de dispersão (k)"),
         x=ifelse(xlab=="p","%CF","grau de limitação de dispersão (k)")) +
    facet_wrap(efeito~name,ncol=1) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position="top",
          # legend.direction = "horizontal",
          # legend.background = element_rect(
          #   fill = "transparent",          
          #   color = NA                     
          # ),
          legend.key = element_rect(       
            fill = "transparent",          
            color = NA                     
          ),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          legend.key.size = unit(0.5, "cm"),     
          legend.spacing = unit(0.2, "cm"),      
          legend.margin = margin(2, 2, 2, 2),
          strip.text = element_text(size=textsize,
                                    margin=margin(t=stripspace,
                                                  b=stripspace)),
#          aspect.ratio=1,#
          axis.text = element_text(size=textsize),
          axis.title = element_text(size=textsize))
}
l_p_apenas_fixo <- dlply(dfpred,"efeito",f_ggplot)
saveRDS(l_p_apenas_fixo,file="./figuras/logUU_construcao/l_p_apenas_fixo.rds")

#grid.arrange(grobs=l_p_apenas_fixo,ncol=3)

#############################
############################################################################################
####################### ANTIGO -> FUTURO ARTIGO: estudo do logOR ###########################
############################################################################################

##########################################################################
### Funções usadas para ajustar os 3 blocos de modelos: te, logU/U e k ###
##########################################################################
## Seguindo as especifícações de Pedersen et al. 2019 https://peerj.com/articles/6876/
# te: tensor entre o efeito na riqueza e a capacidade de dispersão per se
f_gam_te <- \(dfi){
  l_md <- list()
  l_md$`te(logU/U,k)` <- gam(
    logOR ~ 
      te(Uefeito,k_z,
         bs=c("cr","cr"),m=2,
         by=forest_succession,
         id = "fixo") +
      s(lat,long) + 
      t2(Uefeito,k_z,SiteCode,
         bs=c("cr","cr","re"),
         m=2,full=TRUE,
         id = "random"),
    data=dfi,method = "REML")
  l_md$`s(k)` <- gam(
    logOR ~ 
      s(k_z,bs="cr",
        by=forest_succession,id = "fixo") +
      s(lat,long) + 
      te(k_z,SiteCode,bs=c("cr","re"),id = "random"),
    data=dfi,method = "REML")
  saveRDS(l_md,
          file=paste0(v_path,"rds/l_md_",
                      gsub("cont-ideal","fragtotal",dfi$contraste[1]) %>% 
                        gsub("cont-non_frag","fragperse",.) %>% 
                        gsub("non_frag-ideal","areaperse",.),
                      ".rds"))
  rm(l_md);gc()
}
lapply(split(df_md,df_md$contraste),f_gam_te)
#
f_gam_te2 <- \(dfi){
  l_md <- list()
  l_md$`te(logU/U,k) - te(logU/U` <- gam(
    logOR ~ 
      te(Uefeito,k_z,
         bs=c("cr","cr"),m=2,
         by=forest_succession,
         id = "fixo") +
      s(lat,long) + 
      s(SiteCode,bs="re"),
    data=dfi,method = "REML")
  saveRDS(l_md,
          file=paste0(v_path,"rds/l_md_semT2sitio_",
                      gsub("cont-ideal","fragtotal",dfi$contraste[1]) %>% 
                        gsub("cont-non_frag","fragperse",.) %>% 
                        gsub("non_frag-ideal","areaperse",.),
                      ".rds"))
  rm(l_md);gc()
}
lapply(split(df_md,df_md$contraste),f_gam_te2)
#



# s(k)+s(k)|Site
# f_gam3 <- \(dfi){
#   l_md <- list()
#   l_md$`s(k)+s(k)|Site` <- gam(
#     logOR ~
#       s(k_cont,bs="cr",m=2,id = "efeito_comum") +
#       s(k_cont, SiteCode, bs = "fs", xt=list(bs = "cr"), m=2, id="efeito_sitio"),
#     data=dfi,method = "REML")
#   l_md$`s(k)+1|Site` <- gam(
#     logOR ~
#       s(k_cont,bs="cr",m=2,id = "efeito_comum") +
#       s(SiteCode,bs="re"),
#     data=dfi,method = "REML")
#   return(l_md)
# }
#   l_md <- dlply(df_md,"contraste",f_gam3)
#   saveRDS(l_md,paste0(v_path,"rds/l_md_onlyk.rds"))
#   rm("l_md");gc()
#   # s(logU/U)+s(logU/U)|Site
#   f_gam <- \(dfi,bs_type="cr"){
#     l_md <- list()
#     l_md$`s(land)|Site : gs` <- gam(
#       logOR ~
#         s(Uefeito,bs=bs_type,m=2, id="efeito_comum") +
#         s(Uefeito, SiteCode, bs = "fs", xt=list(bs = bs_type), m=2, id="efeito_sitio"),
#       data=dfi,method = "REML")
#     l_md$`s(land) + 1|Site` <- gam(
#       logOR ~
#         s(Uefeito,bs=bs_type,m=2,id="efeito_comum") +
#         s(SiteCode,bs="re"),
#       data=dfi,method = "REML")
#     return(l_md)
#   }
#   l_md_logOR <- dlply(df_md,"contraste",f_gam)
#   saveRDS(l_md_logOR,file=paste0(v_path,"rds/l_md_simples_apudPedersen2019.rds"))
#   rm(l_md_logOR);gc()