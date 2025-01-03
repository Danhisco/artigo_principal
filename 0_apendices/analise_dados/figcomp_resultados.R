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
# dados
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
l_df <- list()
l_df$contraste <- read_csv(file="dados/csv_SoE/df_logOR.csv") %>% 
  select(SiteCode:contraste,logOR:Uefeito) %>%
  rename("logU/U" = "Uefeito") %>% 
  pivot_longer(cols=c("logOR","logU/U"))
df_U <- read_csv(file="dados/csv_SoE/taxaU/df_U.csv")
df_nSAD <- read_csv(file="dados/csv_SoE/df_congruencia_simulacao.csv")
l_df$paisagem <- inner_join(
  select(df_nSAD,SiteCode:nCongKS),
  select(df_U,-Usd)
) %>% 
  filter(SiteCode %in% unique(l_df$contraste$SiteCode)) %>% 
  rename("paisagem"="land_type","n SAD"="nCongKS","U avg"="Umed") %>% 
  mutate(paisagem = gsub("non_frag","aglomer",paisagem)) %>% 
  pivot_longer(.,cols=c("n SAD","U avg"))
v_sites <- lapply(l_df, "[[","SiteCode") %>% 
  sapply(.,unique)
v_sites <- v_sites[,1]
vc <- "SPigua1"
f_todasasRespostas_bySite <- \(vc){
  ldf <- lapply(l_df,filter,SiteCode==vc)
  f_geom_final <- \(vdf){
    if(vdf=="paisagem"){
      list(scale_color_manual("paisagem",
                              values=c("cont"="darkred",
                                       "ideal"="darkgreen",
                                       "aglomer"="darkblue")
                              )
           )
    }else{
      list(scale_color_manual("contraste",
                              values=c("Frag. total"="darkorange",
                                       "Área per se"="#301934",
                                       "Frag. per se"="#158F6E")
                              )
      )
    }
  }
  f_ggplot <- \(dfi,ncolor){
    dfi %>% 
      ggplot(aes(x=k,
                 y=value,
                 color=.data[[ncolor]],
                 group=.data[[ncolor]])) +
      geom_point() +
      geom_line() +
      f_geom_final(ncolor) +
      labs(y="") +
      facet_wrap(as.formula(paste0("~","name")),
                 scales="free") +
      theme(legend.position="top")
  }
  lp <- lapply(names(ldf),\(li){
    f_ggplot(dfi=ldf[[li]],ncolor=li)
  })
  names(lp) <- names(ldf)
  grid.arrange(grobs=lp,ncol=1)
}


## funções de ajuste e de plot





##### figuras complementares dos resultados
lapply(l_path$te,f_plot_te2)
f_juntapreditos <- \(vsite){
  lpng <- list.files(path="figuras/predito_te_sites",
                     pattern=vsite,
                     full.names = T)
  # if(length(lpng)%in%c(0:1)) return("sítio feito")
  img_return <- lapply(lpng,image_read)
  img_return <- image_append(do.call("c",img_return),stack = FALSE) %>% 
    image_resize("50%")
  image_write(img_return,
              path = paste0("figuras/predito_te_sites/",
                            vsite,".png"))
  file.remove(lpng)
  rm(img_return);gc()
}
#### criação da terceira linha: se o hgam +lik é te então fazer dois plots por efeito:
## a) apenas fixo: filtrar k mais próximo dos 20 simulados
## b) fixo + por sítio: plotar em função de k e logU/U para mostrar a variabilidade entre sítios
#
# objetos comuns
l_paths <- paste0("rds/l_dfnew_",c("areaperse","fragperse","fragtotal"),".rds")
l_dfnew <- lapply(l_paths,readRDS)
names(l_dfnew) <- c("areaperse","fragperse","fragtotal")
l_dfpred <- lapply(gsub("dfnew","dfpred",l_paths),readRDS)
names(l_dfpred) <- c("areaperse","fragperse","fragtotal")
# funções
f_imgpng <- \(list_ggplot,vw=5,vh=8){
  l_png <- lapply(list_ggplot,\(li) ggsave(tempfile(fileext = ".png"),plot = li,
                                           width = vw,height=vh))
  lapply(l_png,\(li) image_read(li) %>% image_trim %>% image_resize("50%"))  
}
# f_plot1: logOR ~ Xi (~cut(Xj))
f_plot1 <- \(veffect,ldfref=l_dfnew){
  dff <- l_dfnew[[veffect]]
  k_pred <- dff$k_cont %>% unique
  k_sim <- sapply(c(0.05,0.25,0.50,0.75,0.95),\(x){
    k_pred[which.min(abs(k_pred - x))]
  })
  quantU <- quantile(dff$Uefeito,probs=c(0.05,0.25,0.50,0.75,0.95))
  lp <- list()
  fggplot <- \(dfp,vx,vcolor,vposition){
    dfp %>% 
      mutate(label=gsub("Uefeito","logU/U",vx) %>% gsub("k_cont","k",.)) %>% 
      ggplot(aes(x=.data[[vx]],y=Q_0.5,
                 color=.data[[vcolor]],group=.data[[vcolor]])) +
      geom_hline(yintercept = 0,color="darkgray") +
      geom_vline(xintercept = 0,color="darkgray") +
      geom_ribbon(aes(ymax=Q_0.95,ymin=Q_0.05,
                      fill=.data[[vcolor]],color=.data[[vcolor]]),alpha=0.2) +
      geom_line() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      labs(y="logOR",
           x=gsub("Uefeito","logU/U",vx) %>% 
             gsub("k_cont","k",.), 
           color = gsub("Uefeito","logU/U",vcolor) %>% 
             gsub("k_cont","k",.)) +
      guides(fill = "none") +
      theme_classic() +
      theme(legend.position = "inside",
            legend.position.inside = vposition,
            aspect.ratio = 1) +
      facet_wrap(~label)
  }
  lp$Uefeito <- fggplot(dfp = dff[k_cont%in%k_sim,], #dff[k_cont%in%k_sim,]
                        vx = "Uefeito",
                        vcolor = "k_cont",
                        vposition=c(0.15,0.30))
  lp$k_cont <- fggplot(dfp = dff[Uefeito%in%quantU,],#dff[Uefeito%in%quantU,]
                       vx = "k_cont",
                       vcolor = "Uefeito",
                       vposition=c(0.5,0.30))
  arrangeGrob(grobs=lp,ncol=1,
              top="Apenas Fixo")
}
l_apenasfixo <- lapply(names(l_dfnew),f_plot1)
names(l_apenasfixo) <- names(l_dfnew)
# grid.arrange(grobs=l_apenasfixo,ncol=3,
#              top=NULL,bottom = NULL, left = NULL, right = NULL)
l_png_fixo <- f_imgpng(l_apenasfixo)
names(l_png_fixo) <- names(l_apenasfixo)
# img_final <- image_append(do.call("c",l_png),stack = FALSE)
# f_plot2: logOR ~ Xi e ~ Xj (por sítio)
f_plot2 <- \(veffect,ldfref=l_dfpred){
  df_pred <- ldfref[[veffect]][["fixo e aleat"]]
  fggplot <- \(xvar,dfi=df_pred){
    mutate(dfi,
           label=gsub("k_cont","k",xvar) %>% 
             gsub("Uefeito","logU/U",.)) %>% 
      ggplot(aes(x=.data[[xvar]],y=Q_0.5,group=SiteCode)) +
      geom_hline(yintercept = 0,color="darkgray") +
      geom_vline(xintercept = 0,color="darkgray") +
      geom_ribbon(aes(ymax=Q_0.95,ymin=Q_0.05),alpha=0.2,color="#986868",fill="#986868") +
      geom_line() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      labs(x=gsub("k_cont","k",xvar) %>% 
             gsub("Uefeito","logU/U",.),
           y="logOR") +
      theme_classic() +
      facet_wrap(~label)
  }
  lp <- lapply(c("Uefeito","k_cont"),fggplot)
  names(lp) <- c("Uefeito","k_cont")
  arrangeGrob(grobs=lp,ncol=1,
              top="fixo e por sítio")
}
l_fixoEsitio <- lapply(names(l_dfpred),f_plot2)
names(l_fixoEsitio) <- names(l_dfpred)
# grid.arrange(grobs=l_apenasfixo,ncol=3,
#              top=NULL,bottom = NULL, left = NULL, right = NULL)
l_png_fixoEsitio <- f_imgpng(l_fixoEsitio)
names(l_png_fixoEsitio) <- names(l_fixoEsitio)
# f_figura final
l_png <- lapply(names(l_png_fixo),\(li){
  lp <- list()
  lp$fixo <- l_png_fixo[[li]]
  lp$fixoEsitio <- l_png_fixoEsitio[[li]]
  img_final <- image_append(do.call("c",lp),stack = FALSE)
  img_final <- image_title(img_final,
                           vtitle = gsub("areaperse","Área per se",li) %>% 
                             gsub("fragperse","Frag. per se",.) %>% 
                             gsub("fragtotal","Frag. total",.),
                           vheight = 50,
                           vsize = 40) %>% 
    image_resize("50%")
  return(img_final)
})
names(l_png) <- names(l_png_fixo)
lapply(names(l_png),\(li){
  image_write(l_png[[li]],
              path=paste0("figuras/3alinha_",li,".png"))
})








################ antigo ? ################
#
# lpng <- list.files(path="figuras/predito_te_sites",pattern=".png") %>% 
#   str_extract(pattern="(?<=_)(.*?)(?=.png)")
# table(lpng) %>% table()
#
df_site <- readRDS(file="rds/l_dfpred_areaperse_Ugs.rds")
v_site <- df_site$`fixo e aleat` %>% 
  # filter(SiteCode!="SPigua1") %>% 
  pull(SiteCode) %>% levels
lapply(v_site,f_juntapreditos)
# f_juntapreditos(vsite)
# for(i in v_site){
#   f_juntapreditos(i)
# }
# # ficaram com problemas!...
v_sitesaud <- c("SCside","SCilho","MGvico1")
l_path$te
f_aud <- \(vsite){
  f_porefeito <- \(veffect){
    #
    vname <- str_extract(veffect,pattern_extract) %>% 
      gsub("areaperse","Área per se",.) %>% 
      gsub("fragperse","Frag. per se",.) %>% 
      gsub("fragtotal","Frag. total",.)
    #
    dfpred <- readRDS(paste0(v_path,veffect))$`fixo e aleat` %>% 
      mutate(SiteCode=factor(SiteCode)) %>% select(-any_of("logOR"))
    dfpred <- pivot_longer(dfpred,starts_with("Q_"),
                           names_to="quantiles",
                           values_to="logOR") %>% 
      mutate(quantiles=factor(100 * as.numeric(gsub("Q_","",quantiles)),
                              levels=c(5,50,95)) ) %>% 
      rename(k = k_cont) %>% 
      filter(SiteCode==vsite)
    #
    l_p <- dlply(dfpred,"quantiles",\(dfi){
      dfi %>% 
        ggplot(aes(x=k,y=Uefeito,fill=logOR)) +
        geom_point(shape=21,size=15) +
        scale_fill_viridis_c(name = "logOR",
                             option = "magma") +
        labs(y="logU/U") +
        facet_wrap(~quantiles) +
        theme_classic() +
        theme(legend.position="top",
              aspect.ratio=1,
              legend.text = element_text(angle=90,size = 15),
              legend.title = element_text(size = 16,vjust = 0.75),
              legend.key.size = unit(1.0, 'cm'),
        )
    })
    l_p <- lapply(l_p,\(li) ggsave(tempfile(fileext = ".png"),plot = li,
                                   width = 5, height = 7))
    l_p <- lapply(l_p,\(li) image_read(li) %>% image_trim)
    # base
    img_final <- image_append(do.call("c",l_p[c("5","95")]))
    # padronizacao da mediana
    rect_info <- image_info(img_final)
    rect_width <- rect_info$width
    rect_height <- rect_info$height
    longest_side <- max(rect_width, rect_height)
    l_p[["50"]] <- image_resize(l_p[["50"]], geometry_size_pixels(longest_side))
    img_final <- image_append(c(l_p[["50"]],img_final),stack = TRUE)
    img_final <- image_annotate(
      img_final,
      text = vsite,
      location = "-1000+300",
      gravity = "north",
      size=180,
      boxcolor="lightblue",
      color = "black")
    img_final <- image_resize(img_final,"50%")
    img_final <- image_annotate(
      img_final,
      text = vname,
      location = "-525+15",
      gravity = "north",
      size=80,
      boxcolor="lightgreen",
      color = "black")
    vpath <- image_write(img_final,path = tempfile(fileext = ".png"))
    rm(img_final);gc()
    return(vpath)
  }
  l_img <- lapply(l_path$te,f_porefeito)
  l_img <- lapply(l_img, image_read) %>% 
    lapply(.,image_trim)
  img_return <- image_append(do.call("c",l_img),stack = FALSE)
  image_write(img_return,
              path=paste0("figuras/predito_te_sites/",vsite,".png"))
}
lapply(v_sitesaud,f_aud)



##