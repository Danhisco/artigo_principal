# Pacotes
library(gratia)
library(doMC)
library(gridExtra)
library(ggplot2)
library(magick)
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
library(magick)
image_title <- \(imgobj,
                 vtitle,
                 vheight=150,
                 bgcolor="white",
                 vsize=80,
                 vgrav="north",
                 vcolor="black",
                 vloc="+0+20"){
  image_info <- image_info(imgobj)
  image_width <- image_info$width
  white_canvas <- image_blank(
    width = image_width, 
    height = vheight, 
    color = bgcolor)
  imgobj <- image_append(c(white_canvas, 
                           imgobj), 
                         stack = TRUE)
  image_annotate(
    imgobj,
    text = vtitle,
    size = vsize,
    gravity = vgrav,
    color = vcolor,
    location = vloc
  ) %>% image_trim
}
# objetos
## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
df_p <- read_csv("dados/df_p.csv")
setwd(v_path)
source("figuras_e_tabelas.R")
# source("figuras_e_tabelas.R")
#
####################### 1a linha
l_figfinal <- list()
#
l_path <- paste0(v_path,"rds/l_dfpred_",c("fragtotal","fragperse","areaperse"),".rds")
l_df_pred <- lapply(l_path,readRDS) %>% 
  lapply(.,"[[","apenas fixo")
names(l_df_pred) <- c("Frag. total","Frag. per se","Área per se")
# ii) filtrar os valores únicos de k em um novo data frame
df_ref <- lapply(l_df_pred,select,k_cont,SiteCode) %>% lapply(.,distinct)
df_ref <- df_ref[[1]]
l_md <- readRDS(file=paste0(v_path,"rds/l_md_refU.rds"))
l_df_ref <- lapply(names(l_df_pred),\(li){
  lmd <- l_md[grep(li,names(l_md))]
  names(lmd) <- gsub(paste0(li,"."),"",names(lmd)) 
  df_ref <- lapply(names(lmd),\(i){
    md <- lmd[[i]]
    dfr <- l_df_ref[[li]]  
    dfr[[i]] <- predict.gam(md,dfr)
    return(dfr)
  }) %>% Reduce("inner_join",.)
})
names(l_df_ref) <- names(l_df_pred)
df_ref <- lapply(names(l_df_ref),\(li){
  mutate(l_df_ref[[li]],
         name=gsub("Área per se","area",li) %>% 
           gsub("Frag. total","frag.total",.) %>% 
           gsub("Frag. per se","frag.perse",.))
}) %>% rbindlist() %>% rename(k=k_cont) %>% 
  mutate(label=case_when(
    grepl("area",name) ~ "Área per se",
    grepl("total",name) ~ "Frag. Total",
    grepl("perse",name) ~ "Frag. per se"),
    label=factor(label,levels=c("Frag. Total","Frag. per se","Área per se"))
  )
l_figfinal$`1alinha` <- df_md %>% 
  mutate(label=case_when(
    grepl("area",name) ~ "Área per se",
    grepl("total",name) ~ "Frag. Total",
    grepl("perse",name) ~ "Frag. per se"),
    label=factor(label,levels=c("Frag. Total","Frag. per se","Área per se"))
  ) %>% 
  ggplot(aes(x=k,y=value,group=SiteCode,color=p)) +
  geom_boxplot(inherit.aes = FALSE,
               aes(x=k,y=value,group=k)) +
  geom_hline(yintercept = 0,color="black",linetype=3) +
  geom_hline(yintercept = v_hline,color="darkred") + 
  geom_line(alpha=0.75) +
  geom_point(alpha=0.75) +
  geom_line(data=df_ref,aes(y=max,x=k),color="black") +
  geom_line(data=df_ref,aes(y=min,x=k),color="black") +
  scale_colour_gradient2("% CF",midpoint=0.5,
                         low="red",
                         mid = "yellow",
                         high = "darkgreen") +
  labs(x="k (prop. de propágulos até a vizinhança imediata)",
       y="log(U/U)") +
  scale_y_continuous(expand=c(0.01,0.01)) +
  theme_classic() +
  theme(plot.margin=unit(c(0,0.2,0,0), "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.49,0.9),
        legend.direction="horizontal") +
  # guides(color=guide_legend(direction='horizontal')) +
  facet_wrap(~label,ncol=3)
ggsave(paste0(v_path,"figuras/pedacofigfinal_1alinha.png"),
       l_figfinal$`1alinha`,
       width = 12,
       height = 5)
img_obj <- image_read(paste0(v_path,"figuras/pedacofigfinal_1alinha.png")) %>% 
  image_trim() %>% 
  image_resize("50%")
image_write(img_obj,paste0(v_path,"figuras/pedacofigfinal_1alinha.png"))










## predição a posteriori 
#
### te for every model
l_path <- list()
# l_path$te <-  paste0("rds/l_dfpred_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$te <-  paste0("rds/l_dfnew_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$U <- "rds/l_dfpred_areaperse_Ugs.rds"
df_tabsel <- read_csv("rds/df_tabsel_geral.csv") %>% 
  filter(dAICc==0)
##############
### 1a parte: figfinal_te_EFEITO PAISAGEM
### @descrição: heatmap x=k, y=logU/U, z=
#v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
df_ref <- read_csv(paste0(v_path,
                           "rds/df_limites_predicao_aposteriori_logOR.csv"))
vrange <-  df_ref %>% 
  summarise(max=max(max),min=min(min)) %>% 
  unlist
setwd(v_path)
f_plot_te <- \(veffect,
               vrange=v_range,
               pattern_extract="(?<=l_dfnew_)(.*?)(?=\\.rds)"){
  # veffect <- l_path$te[3]
  f_ggplot_main <- \(dfiQ50){
    dfref <- rename(l_df_ref[[vname]],k=k_cont)
    ggplot(dfi,aes(x=k,y=Uefeito,z=logOR)) +
      geom_raster(aes(fill=logOR)) +
      geom_contour(color = "black") +
      geom_text_contour() +
      geom_line(data=dfref,
                aes(y=max,x=k),color="black",
                inherit.aes = FALSE,
                linewidth=1) +
      geom_line(data=dfref,
                aes(y=min,x=k),color="black",
                inherit.aes = FALSE,
                linewidth=1) +
      scale_fill_viridis_c(name = "logOR",
                           option = "magma") +
      labs(y="logU/U") +
      facet_wrap(~quantiles) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0),
                         limits=vrange[2:1]) +
      theme(legend.position="none",
            strip.text = element_text(size=15,margin=margin()),
            aspect.ratio=1,
            axis.text = element_text(size=15),
            axis.title = element_text(size=15))
  }
  f_ggplot_4quantiles <- \(dfpred){
    k_pred <- dfpred$k %>% unique
    k_sim <- sapply(c(0.25,0.50,0.75,0.99),\(x){
      k_pred[which.min(abs(k_pred - x))]
    })
    dfpred %>% 
      filter(k%in%k_sim) %>% 
      mutate(kf=factor(paste0("k = ",round(k,2)),
                       levels=paste0("k = ",c(0.25,0.50,0.75,0.99))
                       )
             ) %>% 
      pivot_wider(names_from=quantiles,values_from = logOR) %>% 
      ggplot(aes(x=Uefeito,y=`50`)) +
      geom_ribbon(aes(ymin=`5`,ymax=`95`),
                  fill="darkgreen",
                  alpha=0.3) +
      geom_line(color="black",alpha=0.7) +
      ylab("logOR") +
      xlab("logU/U") +
      facet_wrap(~kf,ncol=2) +
      theme(strip.text = element_text(size=15,margin=margin()),
            axis.title = element_text(size=15),
            axis.text = element_text(size=15))
  }
  
  
  library(metR)
  vname <- str_extract(veffect,pattern_extract) %>% 
    gsub("areaperse","Área per se",.) %>% 
    gsub("fragperse","Frag. per se",.) %>% 
    gsub("fragtotal","Frag. total",.)
  # dfpred <- readRDS(veffect)$`apenas fixo`
  dfpred <- readRDS(veffect)
  dfpred <- pivot_longer(dfpred,starts_with("Q_"),
                         names_to="quantiles",
                         values_to="logOR") %>% 
    mutate(quantiles=factor(100 * as.numeric(gsub("Q_","",quantiles)),
                            levels=c(5,50,95)) ) %>% 
    rename(k = k_cont)
  #
  p50 <- f_ggplot_main(filter(dfpred,quantiles=="50"))
  
  
  
  
  
  l_p <- dlply(dfpred,"quantiles",\(dfi){
    
  })
  l_p <- lapply(l_p,\(li) ggsave(tempfile(fileext = ".png"),plot = li))
  l_p <- lapply(l_p,\(li) image_read(li) %>% image_trim)
  # base
  img_final <- image_append(do.call("c",l_p[c("5","95")]))
  # padronizacao da mediana
  rect_info <- image_info(img_final)
  rect_width <- rect_info$width
  rect_height <- rect_info$height
  longest_side <- max(rect_width, rect_height)
  l_p[["50"]] <- image_resize(l_p[["50"]], geometry_size_pixels(longest_side))
  gc()
  # 3 quadros: mediana e os quantis extremos
  img_final <- image_append(c(l_p[["50"]],img_final),stack = TRUE)
  rm(l_p);gc()
  # título
  img_final <- image_title(img_final,vtitle=vname,vsize = 100)
  gc()
  # salvamento
  img_final <- image_resize(img_final, "50%")
  image_write(img_final,
              path=paste0(v_path,"figuras/figfinal_te_",
                          str_extract(veffect,pattern_extract),
                          ".png")
                )
  # 
  image_destroy(img_final)
  rm(img_final);gc()
}
lapply(l_path$te,f_plot_te)
#
#### 
path_ldf = "rds/l_dfpred_simples_apudPedersen2019.rds"
l_df <- readRDS(path_ldf)
# criação dos gráficos de PI
p <- f_plotPI_shgam("Área per se")
ggsave("figuras/figfinal_areaperse_shgam.png",p,
       width = 6,height = 8)
img_p <- image_trim(image_read("figuras/figfinal_areaperse_shgam.png"))
img_p <- image_title(img_p,vtitle="Área per se, ~ logU/U",vsize = 60)
image_write(img_p,path = "figuras/figfinal_areaperse_shgam.png")
#
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