# Pacotes
library(gratia)
library(doMC)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(magick)
theme_set(theme_bw())
library(readr)
library(stringr)
library(tidyr)
library(bbmle)
library(DHARMa)
library(data.table)
library(mgcv)
library(plyr)
library(dplyr)
library(magick)
f_z <- function(x) (x-mean(x))/sd(x)
f_imagefunc <- \(ggobj){
  tempfile_path <- tempfile(fileext = ".png")
  ggsave(tempfile_path,ggobj,width = 7.17,height = 7.20)
  img <- image_trim( image_read(tempfile_path))
  image_write(img,tempfile_path)
}
f_resize_2rectangle <- \(ref_img,toresize_img,ref_side,tore_side){
  tore_info <- image_info(toresize_img)
  ref_info <- image_info(ref_img)
  refside_info <- ref_info[[ref_side]]
  v_command <- ifelse(tore_side=="height",
                      paste0("x",refside_info),
                      paste0(refside_info,"x"))
  image_resize(toresize_img,v_command)
}
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
df_p <- read_csv("dados/df_p.csv")
df_sim <- read_csv("dados/df_simulacao.csv") |> 
  inner_join(x=df_p,by="SiteCode")
df_md <- read_csv("dados/csv_SoE/df_logOR.csv")
## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
# setwd(v_path)
# source(paste0(v_path,"figuras_e_tabelas.R"))
# source("figuras_e_tabelas.R")
#
####################### 1a linha
#
l_path <- paste0(v_path,"rds/l_dfpred_",c("fragtotal","fragperse","areaperse"),".rds")
l_df_pred <- lapply(l_path,readRDS) %>% 
  lapply(.,"[[","apenas fixo") %>% 
  lapply(.,rename,k=k_cont)
names(l_df_pred) <- case_when(
  grepl("fragtotal",l_path) ~ "Frag. total",
  grepl("fragperse",l_path) ~ "Frag. per se",
  grepl("areaperse",l_path) ~ "Área per se"
)
l_df_newpred <- lapply(
  gsub("dfpred","dfnew",l_path),
  readRDS
)
names(l_df_newpred) <- case_when(
  grepl("fragtotal",l_path) ~ "Frag. total",
  grepl("fragperse",l_path) ~ "Frag. per se",
  grepl("areaperse",l_path) ~ "Área per se"
)


# ii) filtrar os valores únicos de k em um novo data frame
l_md <- readRDS(file=paste0(v_path,"rds/l_md_refU.rds"))
l_df_ref <- lapply(names(l_df_pred),\(li){
  lmd <- l_md[grep(li,names(l_md))]
  names(lmd) <- gsub(paste0(li,"."),"",names(lmd)) 
  df_ref <- lapply(names(lmd),\(i){
    md <- lmd[[i]]
    dfr <- l_df_pred[[li]]  
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
}) %>% rbindlist() %>% 
  mutate(label=case_when(
    grepl("area",name) ~ "Área per se",
    grepl("total",name) ~ "Frag. total",
    grepl("perse",name) ~ "Frag. per se"),
    label=factor(label,levels=c("Frag. total","Frag. per se","Área per se"))
  ) %>% select(-name) %>% 
  as.data.frame()
v_sites_RefNulo <- df_md %>% filter(p>=0.975) %>% pull(SiteCode) %>% unique
v_range_RefNulo <- df_md %>% filter(p>=0.975) %>% pull(Uefeito) %>% range
v_hline <- filter(df_md,p>=0.95) %>% 
  pull(Uefeito) %>% quantile(.,c(0.25,0.50,0.75))
# figuras per se
l_figfinal <- list()


l_figfinal$`1alinha` <- df_md %>%
  mutate(across(c(SiteCode,contraste),factor),
         label = contraste) %>% 
  ggplot(aes(x=k,y=Uefeito,group=SiteCode,color=p)) +
  geom_boxplot(inherit.aes = FALSE,
               aes(x=k,y=Uefeito,group=k)) +
  geom_hline(yintercept = 0,color="black",linetype=3) +
  geom_hline(yintercept = v_hline,color="darkred") + 
  geom_line(alpha=0.75) +
  geom_point(alpha=0.75) +
  geom_line(inherit.aes = FALSE,
            data=df_ref,aes(y=max,x=k),color="black") +
  geom_line(inherit.aes = FALSE,
            data=df_ref,aes(y=min,x=k),color="black") +
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
l_path$te <-  paste0(v_path,
                     "rds/l_dfnew_",
                     c("areaperse","fragperse","fragtotal"),
                     ".rds")
l_path$U <- paste0(v_path,
                   "rds/l_dfpred_areaperse_Ugs.rds")
df_tabsel <- read_csv(paste0(v_path,"rds/df_tabsel_geral.csv")) %>% 
  filter(dAICc==0)
##############
### 1a parte: figfinal_te_EFEITO PAISAGEM
### @descrição: heatmap x=k, y=logU/U, z=
#v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
df_ref <- read_csv(paste0(v_path,
                           "rds/df_limites_predicao_aposteriori_logOR.csv"))
v_range <-  df_ref %>% 
  summarise(max=max(max),min=min(min)) %>% 
  unlist
setwd(v_path)
f_plot_te <- \(veffect,
               vrange=v_range,
               pattern_extract="(?<=l_dfnew_)(.*?)(?=\\.rds)",
               legendfillposition=c(0.5,0.1),
               stripspace=0.5,
               stripspace_fa=1,
               textsize=15,
               textsize_4q=10,
               ctextsize=10,
               facetspace_4q=0.5,
               propred_fixo = 0.9,xfixo=0.447, yfixo=0.68,
               propred_fa = 0.7, xfa=0.589, yfa=0.68,
               vw=7,vh=7){
  # veffect <- l_path$te[3]
  f_ggplot_main <- \(dfiQ50,
                     linew=1
                     ){
    dfref <- rename(l_df_ref[[vname]],k=k_cont)
    dfiQ50 %>% 
      mutate(quantiles=paste0(vname,": mediana")) %>% 
      ggplot(aes(x=k,y=Uefeito,z=logOR)) +
      geom_raster(aes(fill=logOR)) +
      geom_contour(color = "black",linewidth=linew) +
      geom_text_contour(size=ctextsize,color="black") +
      geom_line(data=dfref,
                aes(y=max,x=k),color="black",
                inherit.aes = FALSE,
                linewidth=linew) +
      geom_line(data=dfref,
                aes(y=min,x=k),color="black",
                inherit.aes = FALSE,
                linewidth=linew) +
      scale_fill_viridis_c(name = "logOR",
                           option = "magma") +
      labs(y="logU/U") +
      facet_wrap(~quantiles) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0),
                         limits=vrange[2:1]) +
      theme(legend.position=legendfillposition,
            legend.direction = "horizontal",
            legend.background = element_rect(
              fill = "transparent",          
              color = NA                     
            ),
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
            aspect.ratio=1,
            axis.text = element_text(size=textsize),
            axis.title = element_text(size=textsize))
  }
  f_ggplot_4quantiles <- \(dfpred,
                           tsize=textsize_4q){
    k_pred <- dfpred$k %>% unique
    k_sim <- sapply(c(0.25,0.50,0.75,0.99),\(x){
      k_pred[which.min(abs(k_pred - x))]
    })
    dfpred %>% 
      filter(k%in%k_sim) %>% 
      mutate(kf=factor(paste0("k = ",round(k,2)),
                       levels=paste0("k = ",c(0.25,0.50,0.75,0.99)))
             ) %>% 
      pivot_wider(names_from=quantiles,values_from = logOR) %>% 
      ggplot(aes(x=Uefeito,y=`50`)) +
      geom_hline(yintercept = 0,color="darkgray",alpha=0.3) +
      geom_vline(xintercept = 0,color="darkgray",alpha=0.3) +
      geom_ribbon(aes(ymin=`5`,ymax=`95`),
                  fill="darkgreen",
                  alpha=0.3) +
      geom_line(color="black",alpha=0.7) +
      ylab("logOR") +
      xlab("logU/U") +
      facet_wrap(~kf,ncol=2) +
      theme(aspect.ratio = 1,
            strip.text = element_text(size=tsize,
                                      margin=margin(t=stripspace,
                                                    b=stripspace)),
            axis.title = element_text(size=tsize),
            axis.text = element_text(size=tsize),
            panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            plot.margin = margin(),
            panel.spacing = unit(facetspace_4q, "cm", data = NULL))
  }
  f_fixo_e_aleatorio <- \(dfi,
                          vspacestrip=stripspace_fa,
                          vtextsize=textsize,
                          striptextsize=15){
    dfi %>% 
      mutate(label="fixo e aleatório") %>% 
      ggplot(aes(x = Uefeito, group = SiteCode)) +
      # quantiile interval
      geom_ribbon(aes(ymin = Q_0.05, ymax = Q_0.95, 
                    fill = "Quantile range", group = SiteCode), 
                alpha = 0.2) +
      geom_line(aes(y = Q_0.5, color = "Median", group = SiteCode), 
                alpha = 0.5) +
      scale_fill_manual(values = c("Quantile range" = "#986868")) +
      scale_color_manual(values = c("Median" = "#C04000")) +
      # 0 x 0 
      geom_hline(yintercept = 0,color="darkgray",alpha=0.75) +
      geom_vline(xintercept = 0,color="darkgray",alpha=0.75) +
      # ajustes
      scale_x_continuous(expand = expansion(add = c(0,0))) +
      scale_y_continuous(expand = expansion(add = c(0,0))) +
      ylab("logOR") +
      xlab("logU/U") +
      facet_wrap(~label) +
      theme_classic() +
      theme(
        axis.title.x = element_text(hjust = 0.5, vjust = 0.5,margin = margin(t = -2.5)),
        axis.title.y = element_text(hjust = 0.5, vjust = 0.5,margin = margin(t = -30)),
        legend.position = "none",
        aspect.ratio = 1,
        axis.text = element_text(size=vtextsize),
        axis.title = element_text(size=vtextsize),
        strip.text = element_text(size=striptextsize,
                                  margin=margin(t=vspacestrip,
                                                b=vspacestrip)),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)
      )
  }
  library(metR)
  theme_set(theme_classic())
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
  # apenas fixo
  p50 <- f_ggplot_main(dfiQ50 = filter(dfpred,quantiles=="50"))
  p4quan <- f_ggplot_4quantiles(dfpred = dfpred)
  p_final0 <- ggdraw() +
    draw_plot(p50) +
    draw_plot(p4quan,
              height=0.4*propred_fixo,
              width=0.25*propred_fixo,
              x=xfixo,y=yfixo) 
  # incluir o fixo com aleatório
  df_fa <- readRDS(gsub("dfnew","dfpred",veffect))
  df_fa <- df_fa[["fixo e aleat"]]
  p_fa <- f_fixo_e_aleatorio(df_fa,vtextsize = 10,striptextsize = 10)
  p_final <- ggdraw() +
    draw_plot(p_final0) +
    draw_plot(p_fa,
              height=0.4*propred_fa,
              width=0.25*propred_fa,
              x=xfa,y=yfa)
  vpath <- ggsave(plot=p_final,
                  filename=tempfile(fileext=".png"),
                  width = vw,height = vh)
  return(vpath)
}
####### salvamento da amostra
l_img <- lapply(l_path$te,\(li){
  f_plot_te(veffect = li,
            textsize_4q = 10,
            ctextsize=5,
            legendfillposition=c(0.5,0.05),
            xfixo = 0.45, yfixo = 0.610, propred_fixo = 1,
            facetspace_4q = 0.1,
            propred_fa = 1.175, xfa = 0.70, yfa=0.58)
  }
)
names(l_img) <- gsub("/rds/","/figuras/",l_path$te) %>% 
  gsub("l_dfnew_","figfinal_",.) %>% 
  gsub("rds$","jpeg",.)
lapply(names(l_img),\(li){
  img <- image_read(l_img[[li]]) %>% 
    image_trim() %>% 
    image_resize("20%")
  image_write(img,path = li)
})
#
#
############## plot do modelo mais plausível para área per se:
f_plotPI_shgam <- \(nefeito,
                    stripspace=0.5,
                    striptextsize=10,
                    textsize=15){
  # objeto para o gráfico
  v_range_x <- sapply(l_df,\(li){
    range(li[["apenas fixo"]]$Uefeito)
  }) %>% range()
  v_range_y <- sapply(l_df,\(li){
    range(select(li[["apenas fixo"]],starts_with("Q_")))
  }) %>% range()
  f_geom_legend <- list(
    # guide name
    annotate("text", x = 0.01, y =-5.5, 
             label = "Posterior Prediction Interval", hjust = 0) ,
    # mediana
    annotate("segment", x = 0.01, xend = 0.06, y = -6, yend = -6, 
             color = "black", linewidth = 1),
    annotate("text", x = 0.08, y =-6, 
             label = "median", hjust = 0) ,
    # quantile range
    annotate("rect", xmin = 0.01, xmax = 0.06, ymin = -6.75, ymax = -6.25,
             fill = "gray", alpha = 0.7),
    annotate("text", x = 0.08, y = -6.5,
             label = "90% quant. range", hjust = 0),
    # Site
    annotate("rect", xmin = 0.01, xmax = 0.06, ymin = -7.5, ymax = -7,
             fill = "#986868", alpha = 0.7),
    annotate("text", x = 0.08, y = -7.25,
             label = "s(log(U/U)) + s(log(U/U))|Site", hjust = 0),
    # overall
    annotate("rect", xmin = 0.01, xmax = 0.06, ymin = -8.25, ymax = -7.75,
             fill = "darkgreen", alpha = 0.7),
    annotate("text", x = 0.08, y = -8,
             label = "s(log(U/U)) + 0|Site", hjust = 0)
  )
  # padronização dos dados
  l_df[[nefeito]][["apenas fixo"]]$SiteCode <- "apenas fixo"
  ldfi <- lapply(l_df[[nefeito]],mutate,
                 contraste = nefeito,
                 SiteCode = factor(SiteCode)) %>% 
    lapply(.,select,-any_of("logOR"))
  # gráfico
  p <- ldfi[["fixo e aleat"]] %>% 
    mutate(label=nefeito) %>% 
    ggplot(aes(x = Uefeito, group = SiteCode)) +
    # fixo e aleat
    geom_ribbon(aes(ymin = Q_0.05, ymax = Q_0.95, 
                    fill = "Quantile range", group = SiteCode), 
                alpha = 0.2) +
    geom_line(aes(y = Q_0.5, color = "Median", group = SiteCode), 
              alpha = 0.5) +
    scale_fill_manual(values = c("Quantile range" = "#986868")) +
    scale_color_manual(values = c("Median" = "#C04000")) +
    # apenas fixo
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    geom_ribbon(data=ldfi[["apenas fixo"]],
                aes(ymin = Q_0.05, ymax = Q_0.95, 
                    fill = "Quantile range", group = SiteCode), 
                alpha = 0.3) +
    geom_line(data=ldfi[["apenas fixo"]],
              aes(y = Q_0.5, color = "Median", group = SiteCode), 
              alpha = 0.7) +
    scale_fill_manual(values = c("Quantile range" = "darkgreen")) +
    scale_color_manual(values = c("Median" = "black")) +
    # 0 x 0 
    geom_hline(yintercept = 0,color="darkgray",alpha=0.75) +
    geom_vline(xintercept = 0,color="darkgray",alpha=0.75) +
    # ajustes
    scale_x_continuous(expand = expansion(add = c(0,0))) +
    scale_y_continuous(expand = expansion(add = c(0,0))) +
    labs(x="logU/U",y="log OR") +
    facet_wrap(~label) +
    theme_classic() +
    theme(
      axis.title.x = element_text(hjust = 0.5, vjust = 0.5,margin = margin(t = -2.5)),
      axis.title.y = element_text(hjust = 0.5, vjust = 0.5,margin = margin(t = -30)),
      legend.position = "none",
      aspect.ratio = 1,
      axis.text = element_text(size=textsize),
      axis.title = element_text(size=textsize),
      strip.text = element_text(size=striptextsize,
                                margin=margin(t=stripspace,
                                              b=stripspace))
    )
  if(nefeito=="Área per se"){
    p <- p + f_geom_legend
  }
  return(p)
}
f_plot_shgam <- \(efeito="Área per se",
                  vpath="figuras/figfinal_areaperse_shgam.png",
                  vw=6,vh=8,
                  vstripspace=0.5,
                  vstriptextsize=10,
                  vtextsize=15){
  # bases com a predição a posteriori para os modelos com spline simples
  path_ldf = "rds/l_dfpred_simples_apudPedersen2019.rds"
  l_df <- readRDS(path_ldf)
  # criação dos gráficos de PI
  p <- f_plotPI_shgam("Área per se",
                      stripspace = vstripspace,
                      striptextsize=vstriptextsize,
                      textsize=vtextsize)
  ggsave(vpath, p, width = vw, height = vh)
  # padronização
  img_p <- image_trim(image_read(vpath))
  image_write(img_p,path = vpath)
  return(vpath)
}
########## leitura e salvamento da amostra
image_read(f_plot_shgam(vtextsize = 13))
#
#
#
###### combinação das 3 figuras
f_figfinal <- \(df_ajuste,
                vlegpos=c(0.5,0.04),
                vimgrs="75%"){
  # frag. per se e frag total
  l_img <- dlply(filter(df_ajuste,grepl("Frag." ,efeito)),"efeito",\(dfi){
    with(dfi,{
      v_effect <- ifelse(efeito=="Frag. total","fragtotal","fragperse")
      f_plot_te(veffect = grep(v_effect,l_path$te,value=TRUE),
                legendfillposition=vlegpos,
                textsize_4q = tsize4q,
                ctextsize=ctextsize,
                xfixo = xfixo, yfixo = yfixo, propred_fixo = propred_fixo,
                facetspace_4q = facetspace_4q,
                propred_fa = propred_fa, xfa = xfa, yfa = yfa)  
    })
  }) %>% lapply(.,image_read)
  # área em branco
  # l_img$blank <- 
  # área per se
  l_img$`Área per se` <- with(filter(df_ajuste,efeito=="Área per se"),{
      image_read(f_plot_shgam(vtextsize = vtextsize,
                              vw = vw, vh = vh,
                              vstripspace = vstripspace,
                              vstriptextsize = vstriptextsize))
    })
  ## padronização
  l_img <- lapply(l_img,image_trim) %>% 
    lapply(.,image_resize,vimgrs)
  l_img$`Área per se` <- f_resize_2rectangle(ref_img = l_img$`Frag. per se`,
                                             toresize_img = l_img$`Área per se`,
                                             ref_side = "height",
                                             tore_side = "height")
  l_img$`Área per se` <- f_resize_2rectangle(ref_img = l_img$`Frag. per se`,
                                             toresize_img = l_img$`Área per se`,
                                             ref_side = "width",
                                             tore_side = "width")
  # imagem combinada
  image_append(do.call("c",l_img),stack = FALSE)
}
dfajuste <- data.frame(
  efeito = factor(c("Área per se","Frag. per se", "Frag. total"),
                  levels=c("Frag. total","Frag. per se","Área per se")),
  tsize4q = c(NA,10,10),
  ctextsize = c(NA,5,5),
  xfixo = c(NA,0.45,0.45),
  yfixo = c(NA,0.610,0.610), 
  propred_fixo = c(NA,1,1),
  facetspace_4q = c(NA,0.1,0.1),
  propred_fa = c(NA,1.175,1.175), 
  xfa = c(NA,0.70,0.70), 
  yfa=c(NA,0.58,0.58),
  vw=c(7,NA,NA),
  vh=c(9,NA,NA),
  vstripspace=c(0.5,NA,NA),
  vstriptextsize=c(10,NA,NA),
  vtextsize=c(15,NA,NA),
  qrange_fixo="NA",
  median_fixo="NA",
  qrange_aleat="NA",
  median_aleat="NA"
)
img_final <- f_figfinal(dfajuste)
image_write(
  img_final,
  path="/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/figuras/comparacao_de_efeitos.jpeg"
  )
