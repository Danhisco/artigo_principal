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
library(magick)
library(mgcv)
library(plyr)
library(dplyr)
## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
library(mgcv)
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
l_path <- paste0(v_path,"rds/l_dfpred_",
                 c("areaperse",
                   "fragperse",
                   "fragtotal"),
                 ".rds")
l_df <- lapply(l_path,readRDS)
names(l_df) <- gsub("areaperse","Área per se",
                    str_extract(l_path,"(?<=dfpred_)(.*?)(?=\\.rds)")) %>% 
  gsub("fragperse","Frag. per se",.) %>% 
  gsub("fragtotal","Frag. total",.)
l_df <- lapply(l_df,"[[","fixo e aleat")
######### gráficos diagnósticos para o artigo principal
library(hexbin)
f_plot <- \(dfp){
  dfp <- dfp %>% 
    pivot_longer(starts_with("Q_"),
                 names_to="quantile",
                 values_to = "preditor") %>% 
    mutate(quantile=gsub("Q_","",quantile) %>% 
             as.numeric(),
           quantile=as.character(quantile*100))
  fggplot <- \(dfi){
    dfi %>% 
      ggplot(aes(x=logOR,y=preditor)) +
      geom_hex(bins = 50) + 
      scale_fill_viridis_c(option = "plasma") +
      geom_abline(slope=1,intercept=0,color="green") +
      geom_smooth(aes(group=SiteCode),
                  method="lm",
                  se=FALSE,
                  color="gray",
                  alpha=0.2) +
      labs(title=paste("quantile:",dfi$quantile[1])) +
      theme(legend.position = c(0.90,0.25))
 }
  lp <- dlply(dfp,"quantile",fggplot)
  lp <- lapply(lp,f_imagefunc)
  limg <- lapply(lp,\(li) image_trim(image_read(li)))
  base_img <- image_append(do.call("c",limg[c("5","95")]),stack=FALSE)
  base_img <- f_resize_2rectangle(ref_img = limg[["50"]],
                                  toresize_img = base_img,
                                  ref_side = "width",
                                  tore_side = "width")
  img_final <- image_append(c(limg[["50"]],base_img),stack = TRUE) %>% 
    image_resize("50%")
  return(img_final)
}
l_img_diag <- lapply(l_df,f_plot)
names(l_img_diag) <- gsub("rds/l_dfpred","figuras/diagfinal",l_path) %>% 
  gsub("\\.rds","\\.png",.)
lapply(names(l_img_diag),\(li) image_write(l_img_diag[[li]],path=li))


f_dfplot <- \(dff,dataset){
  vcols <- c("logOR","k_cont")
  dff %>% 
    select(-any_of(vcols)) %>% 
    pivot_longer(starts_with("Q_"),names_to="quantile",
                 values_to=paste0("logOR_",dataset)) %>% 
    mutate(quantile=gsub("Q_","",quantile) %>% 
             as.numeric(),
           quantile=as.character(quantile*100))
}
df_te <- lapply(l_df,"[[","fixo e aleat")[["Área per se"]] %>% f_dfplot(dataset="te")
df_s <- readRDS(paste0(v_path,"rds/l_dfpred_areaperse_Ugs.rds"))[["fixo e aleat"]] %>% 
  f_dfplot(dataset="s")
df_plot <- inner_join(df_te,df_s) %>% 
  inner_join(lapply(l_df,"[[","fixo e aleat")[["Área per se"]] %>% 
               select(logOR,Uefeito,SiteCode))
f_ggplot <- \(dfi){
  l_p <- list()
  l_p$predpred <- dfi %>% 
    ggplot(aes(x=logOR_te,y=logOR_s,color=quantile)) +
    geom_abline(intercept = 0,slope=1) +
    geom_point() +
    geom_smooth(method="lm",se=FALSE) +
    scale_color_manual("%",
                       values=c("50"="darkgreen",
                                "5"="darkred",
                                "95"="darkred")) +
    labs(title="predito X predito") +
    theme(legend.position = "top",
          aspect.ratio = 1) +
    facet_wrap(~SiteCode) 
  l_p$obs_uefeito_predicoes <- dfi %>% 
    # preparação dos dados
    pivot_longer(matches("_te|_s"),
                 names_to = "hgam",
                 values_to = "logOR_pred") %>% 
    mutate(hgam=gsub("logOR_","",hgam)) %>% 
    ggplot(aes(x=logOR,y=logOR_pred,color=quantile)) +
    geom_abline(intercept = 0,slope=1) +
    geom_point() +
    geom_smooth(method="lm",se=FALSE) +
    scale_color_manual("%",
                       values=c("50"="darkgreen",
                                "5"="darkred",
                                "95"="darkred")) +
    labs(x="logOR observado",y="logOR predito",
         title=dfi$SiteCode[1]) +
    theme(legend.position = "top",
          aspect.ratio = 1) +
    facet_wrap(~hgam) 
  l_p <- lapply(l_p,\(li) ggsave(tempfile(fileext = ".png"),plot = li))
  l_p <- lapply(l_p,\(li) image_read(li) %>% image_trim)
  img_final <- image_append(do.call("c",l_p),stack = FALSE)
  image_write(
    img_final,
    path=paste0(v_path,
                "figuras/PIeOBS_sites/compara_areaperse/",
                dfi$SiteCode[1],
                ".png")
  )
}
lapply(levels(df_s$SiteCode),\(i){
  filter(df_plot,SiteCode==i) %>% f_ggplot()
})


