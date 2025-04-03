library(gridExtra)
library(ggplot2)
library(magick)
library(flextable)
library(officer)
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
f_resize_2rectangle <- \(ref_img,toresize_img,ref_side,tore_side){
  tore_info <- image_info(toresize_img)
  ref_info <- image_info(ref_img)
  refside_info <- ref_info[[ref_side]]
  v_command <- ifelse(tore_side=="height",
                      paste0("x",refside_info),
                      paste0(refside_info,"x"))
  image_resize(toresize_img,v_command)
}
# dados
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
l_df <- list()
df_info <- read_csv("dados/df_dados_disponiveis.csv") %>% 
  select(SiteCode,effort_ha,Ntotal,S_obs)
l_df$site_info <- readRDS(file="dados/csv_SoE/df_logOR.rds") %>% 
  select(SiteCode,forest_succession:data_year) %>%
  distinct() %>% 
  inner_join(.,df_info)
l_df$contraste <- readRDS(file="dados/csv_SoE/df_logOR.rds") %>% 
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
# vc <- "MGipia1"
f_todasasRespostas_bySite <- \(vc){
  ldf <- lapply(l_df,filter,SiteCode==vc)
  # gráfico dos dados
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
                              ),
           geom_hline(yintercept = 0,color="darkgray",alpha=0.35)
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
      theme(legend.position = "top",
            legend.title = element_text(size=8),
            plot.margin = margin(0,0.75,0,0))
  }
  f_p_extensoes <- \(mat){
    lp <- c()
    lp[1] <- round(sum(mat[mat!=0]) / length(mat),2)
    mat_size <- nrow(mat) 
    # coord centr
    center_x <- mat_size / 2 
    center_y <- mat_size / 2
    # coord 2km x 2km
    square_size_2km <- mat_size / 4 
    x_min_2km <- center_x - square_size_2km 
    x_max_2km <- center_x + square_size_2km
    y_min_2km <- center_y - square_size_2km
    y_max_2km <- center_y + square_size_2km
    mat2 <- mat[x_min_2km:x_max_2km,
                y_min_2km:y_max_2km]
    lp[2] <- round(sum(mat2[mat2!=0]) / length(mat2),2)
    # coord 1km x 1km
    square_size_1km <- mat_size / 8
    x_min_1km <- center_x - square_size_1km
    x_max_1km <- center_x + square_size_1km
    y_min_1km <- center_y - square_size_1km
    y_max_1km <- center_y + square_size_1km
    mat3 <- mat[x_min_1km:x_max_1km,
                y_min_1km:y_max_1km]
    lp[3] <- round(sum(mat3[mat3!=0]) / length(mat3),2)
    lp[lp>1] <- 1
    #
    names(lp) <- 4/(c(1,2,4))
    return(lp)
  }
  f_plot <- \(mat, filename = "matrix_plot.png") {
    colors <- c("gray", "darkgreen", "darkorange")
    mat_size <- 1 
    # coord centr
    center_x <- mat_size / 2 
    center_y <- mat_size / 2
    # coord 2km x 2km
    square_size_2km <- mat_size / 4 
    x_min_2km <- center_x - square_size_2km 
    x_max_2km <- center_x + square_size_2km
    y_min_2km <- center_y - square_size_2km
    y_max_2km <- center_y + square_size_2km
    # coord 1km x 1km
    square_size_1km <- mat_size / 8
    x_min_1km <- center_x - square_size_1km
    x_max_1km <- center_x + square_size_1km
    y_min_1km <- center_y - square_size_1km
    y_max_1km <- center_y + square_size_1km
    #
    png(filename = filename) 
    #
    image(mat, col = colors, axes = FALSE)
    #
    rect(xleft = x_min_2km, xright = x_max_2km, ybottom = y_min_2km, ytop = y_max_2km, 
         border = "black", lwd = 2)
    #
    rect(xleft = x_min_1km, xright = x_max_1km, ybottom = y_min_1km, ytop = y_max_1km, 
         border = "black", lwd = 2)
    #
    lp <- f_p_extensoes(mat)
    vprop <- 0.1
    text(x = mat_size - 0.05 * vprop, 
         y = mat_size - 0.05 * vprop, 
         labels = paste0("p=", lp["4"]), 
         adj = c(1, 1))
    text(x = x_max_2km - 0.05 * vprop, 
         y = y_max_2km - 0.05 * vprop, 
         labels = paste0("p=", lp["2"]), 
         adj = c(1, 1)) 
    text(x = x_max_1km - 0.05 * vprop, 
         y = y_max_1km - 0.05 * vprop, 
         labels = paste0("p=", lp["1"]), 
         adj = c(1, 1)) 
    dev.off()
  }
  #
  lp <- lapply(names(ldf)[2:3],\(li){
    f_ggplot(dfi=ldf[[li]],ncolor=li)
  })
  names(lp) <- names(ldf)[2:3]
  p <- arrangeGrob(grobs=lp,ncol=1)
  vpath1 <- tempfile(fileext = ".png")
  ggsave(vpath1,p,
         width = 7.5,height = 7)
  img1 <- image_read(vpath1) %>% image_trim
  # mapa de cobertura florestal usado nas simulações
  vpath <- paste0(
    "dados/simulacao/",
    ldf$site_info$SiteCode,
    ".txt"
  ) 
  vtxt <- read.table(vpath) %>% as.matrix()
  vpath <- tempfile(fileext = ".jpeg")
  f_plot(mat = vtxt,
         filename = vpath)
  img2 <- image_read(vpath) %>% 
    image_trim()
  
  img_width_inches <- 
    image_info(img2)$width / as.numeric(strsplit(image_info(img2)$density, 
                                                 "x")[[1]])[1]
  # tabela:
  ft <- flextable(ldf$site_info) %>% 
    autofit() %>% 
    # width(width = unit(img_width_inches, "in")) %>% 
    bg(bg = "lightgrey", part = "header") %>% 
    bg(bg = "white", part = "body") %>% 
    set_header_labels(
      SiteCode = "Código do Sítio",
      forest_succession = "Classe de perturbação",
      lat = "Lat",
      long = "Long",
      data_year = "ano mais próximo",
      effort_ha = "área amostrada (ha)",
      Ntotal = "N",
      S_obs = "S"
    ) %>% 
    fontsize(size = 12, part = "all")
  vpath3 <- tempfile(fileext = ".png")
  save_as_image(ft,vpath3)
  img3 <- image_read(vpath3) %>% image_trim
  img3 <- f_resize_2rectangle(ref_img = img2,
                              toresize_img = img3,
                              ref_side = "width",
                              tore_side  = "width")
  #
  limg <- list()
  limg[[1]] <- image_append(c(img3,img2),stack = TRUE)
  limg[[2]] <- f_resize_2rectangle(ref_img = limg[[1]],
                                   toresize_img = img1,
                                   ref_side = "height",
                                   tore_side = "height")
  img_final <- image_append(
    do.call("c",limg),stack = FALSE
  )
  image_write(img_final,
              path = paste0("dados/csv_SoE/figuras/por_sitio/",
                            vc,".jpeg"))
}
lapply(v_sites,f_todasasRespostas_bySite)
#
#
############################ antigo ############################
## funções de ajuste e de plot
##### figuras complementares dos resultados
# lapply(l_path$te,f_plot_te2)
# f_juntapreditos <- \(vsite){
#   lpng <- list.files(path="figuras/predito_te_sites",
#                      pattern=vsite,
#                      full.names = T)
#   # if(length(lpng)%in%c(0:1)) return("sítio feito")
#   img_return <- lapply(lpng,image_read)
#   img_return <- image_append(do.call("c",img_return),stack = FALSE) %>% 
#     image_resize("50%")
#   image_write(img_return,
#               path = paste0("figuras/predito_te_sites/",
#                             vsite,".png"))
#   file.remove(lpng)
#   rm(img_return);gc()
# }
# #### criação da terceira linha: se o hgam +lik é te então fazer dois plots por efeito:
# ## a) apenas fixo: filtrar k mais próximo dos 20 simulados
# ## b) fixo + por sítio: plotar em função de k e logU/U para mostrar a variabilidade entre sítios
# #
# # objetos comuns
# l_paths <- paste0("rds/l_dfnew_",c("areaperse","fragperse","fragtotal"),".rds")
# l_dfnew <- lapply(l_paths,readRDS)
# names(l_dfnew) <- c("areaperse","fragperse","fragtotal")
# l_dfpred <- lapply(gsub("dfnew","dfpred",l_paths),readRDS)
# names(l_dfpred) <- c("areaperse","fragperse","fragtotal")
# # funções
# f_imgpng <- \(list_ggplot,vw=5,vh=8){
#   l_png <- lapply(list_ggplot,\(li) ggsave(tempfile(fileext = ".png"),plot = li,
#                                            width = vw,height=vh))
#   lapply(l_png,\(li) image_read(li) %>% image_trim %>% image_resize("50%"))  
# }
# # f_plot1: logOR ~ Xi (~cut(Xj))
# f_plot1 <- \(veffect,ldfref=l_dfnew){
#   dff <- l_dfnew[[veffect]]
#   k_pred <- dff$k_cont %>% unique
#   k_sim <- sapply(c(0.05,0.25,0.50,0.75,0.95),\(x){
#     k_pred[which.min(abs(k_pred - x))]
#   })
#   quantU <- quantile(dff$Uefeito,probs=c(0.05,0.25,0.50,0.75,0.95))
#   lp <- list()
#   fggplot <- \(dfp,vx,vcolor,vposition){
#     dfp %>% 
#       mutate(label=gsub("Uefeito","logU/U",vx) %>% gsub("k_cont","k",.)) %>% 
#       ggplot(aes(x=.data[[vx]],y=Q_0.5,
#                  color=.data[[vcolor]],group=.data[[vcolor]])) +
#       geom_hline(yintercept = 0,color="darkgray") +
#       geom_vline(xintercept = 0,color="darkgray") +
#       geom_ribbon(aes(ymax=Q_0.95,ymin=Q_0.05,
#                       fill=.data[[vcolor]],color=.data[[vcolor]]),alpha=0.2) +
#       geom_line() +
#       scale_x_continuous(expand = c(0,0)) +
#       scale_y_continuous(expand = c(0,0)) +
#       labs(y="logOR",
#            x=gsub("Uefeito","logU/U",vx) %>% 
#              gsub("k_cont","k",.), 
#            color = gsub("Uefeito","logU/U",vcolor) %>% 
#              gsub("k_cont","k",.)) +
#       guides(fill = "none") +
#       theme_classic() +
#       theme(legend.position = "inside",
#             legend.position.inside = vposition,
#             aspect.ratio = 1) +
#       facet_wrap(~label)
#   }
#   lp$Uefeito <- fggplot(dfp = dff[k_cont%in%k_sim,], #dff[k_cont%in%k_sim,]
#                         vx = "Uefeito",
#                         vcolor = "k_cont",
#                         vposition=c(0.15,0.30))
#   lp$k_cont <- fggplot(dfp = dff[Uefeito%in%quantU,],#dff[Uefeito%in%quantU,]
#                        vx = "k_cont",
#                        vcolor = "Uefeito",
#                        vposition=c(0.5,0.30))
#   arrangeGrob(grobs=lp,ncol=1,
#               top="Apenas Fixo")
# }
# l_apenasfixo <- lapply(names(l_dfnew),f_plot1)
# names(l_apenasfixo) <- names(l_dfnew)
# # grid.arrange(grobs=l_apenasfixo,ncol=3,
# #              top=NULL,bottom = NULL, left = NULL, right = NULL)
# l_png_fixo <- f_imgpng(l_apenasfixo)
# names(l_png_fixo) <- names(l_apenasfixo)
# # img_final <- image_append(do.call("c",l_png),stack = FALSE)
# # f_plot2: logOR ~ Xi e ~ Xj (por sítio)
# f_plot2 <- \(veffect,ldfref=l_dfpred){
#   df_pred <- ldfref[[veffect]][["fixo e aleat"]]
#   fggplot <- \(xvar,dfi=df_pred){
#     mutate(dfi,
#            label=gsub("k_cont","k",xvar) %>% 
#              gsub("Uefeito","logU/U",.)) %>% 
#       ggplot(aes(x=.data[[xvar]],y=Q_0.5,group=SiteCode)) +
#       geom_hline(yintercept = 0,color="darkgray") +
#       geom_vline(xintercept = 0,color="darkgray") +
#       geom_ribbon(aes(ymax=Q_0.95,ymin=Q_0.05),alpha=0.2,color="#986868",fill="#986868") +
#       geom_line() +
#       scale_x_continuous(expand = c(0,0)) +
#       scale_y_continuous(expand = c(0,0)) +
#       labs(x=gsub("k_cont","k",xvar) %>% 
#              gsub("Uefeito","logU/U",.),
#            y="logOR") +
#       theme_classic() +
#       facet_wrap(~label)
#   }
#   lp <- lapply(c("Uefeito","k_cont"),fggplot)
#   names(lp) <- c("Uefeito","k_cont")
#   arrangeGrob(grobs=lp,ncol=1,
#               top="fixo e por sítio")
# }
# l_fixoEsitio <- lapply(names(l_dfpred),f_plot2)
# names(l_fixoEsitio) <- names(l_dfpred)
# # grid.arrange(grobs=l_apenasfixo,ncol=3,
# #              top=NULL,bottom = NULL, left = NULL, right = NULL)
# l_png_fixoEsitio <- f_imgpng(l_fixoEsitio)
# names(l_png_fixoEsitio) <- names(l_fixoEsitio)
# # f_figura final
# l_png <- lapply(names(l_png_fixo),\(li){
#   lp <- list()
#   lp$fixo <- l_png_fixo[[li]]
#   lp$fixoEsitio <- l_png_fixoEsitio[[li]]
#   img_final <- image_append(do.call("c",lp),stack = FALSE)
#   img_final <- image_title(img_final,
#                            vtitle = gsub("areaperse","Área per se",li) %>% 
#                              gsub("fragperse","Frag. per se",.) %>% 
#                              gsub("fragtotal","Frag. total",.),
#                            vheight = 50,
#                            vsize = 40) %>% 
#     image_resize("50%")
#   return(img_final)
# })
# names(l_png) <- names(l_png_fixo)
# lapply(names(l_png),\(li){
#   image_write(l_png[[li]],
#               path=paste0("figuras/3alinha_",li,".png"))
# })
# 
# 
# 
# 
# 
# 
# 
# 
# ################ antigo ? ################
# #
# # lpng <- list.files(path="figuras/predito_te_sites",pattern=".png") %>% 
# #   str_extract(pattern="(?<=_)(.*?)(?=.png)")
# # table(lpng) %>% table()
# #
# df_site <- readRDS(file="rds/l_dfpred_areaperse_Ugs.rds")
# v_site <- df_site$`fixo e aleat` %>% 
#   # filter(SiteCode!="SPigua1") %>% 
#   pull(SiteCode) %>% levels
# lapply(v_site,f_juntapreditos)
# # f_juntapreditos(vsite)
# # for(i in v_site){
# #   f_juntapreditos(i)
# # }
# # # ficaram com problemas!...
# v_sitesaud <- c("SCside","SCilho","MGvico1")
# l_path$te
# f_aud <- \(vsite){
#   f_porefeito <- \(veffect){
#     #
#     vname <- str_extract(veffect,pattern_extract) %>% 
#       gsub("areaperse","Área per se",.) %>% 
#       gsub("fragperse","Frag. per se",.) %>% 
#       gsub("fragtotal","Frag. total",.)
#     #
#     dfpred <- readRDS(paste0(v_path,veffect))$`fixo e aleat` %>% 
#       mutate(SiteCode=factor(SiteCode)) %>% select(-any_of("logOR"))
#     dfpred <- pivot_longer(dfpred,starts_with("Q_"),
#                            names_to="quantiles",
#                            values_to="logOR") %>% 
#       mutate(quantiles=factor(100 * as.numeric(gsub("Q_","",quantiles)),
#                               levels=c(5,50,95)) ) %>% 
#       rename(k = k_cont) %>% 
#       filter(SiteCode==vsite)
#     #
#     l_p <- dlply(dfpred,"quantiles",\(dfi){
#       dfi %>% 
#         ggplot(aes(x=k,y=Uefeito,fill=logOR)) +
#         geom_point(shape=21,size=15) +
#         scale_fill_viridis_c(name = "logOR",
#                              option = "magma") +
#         labs(y="logU/U") +
#         facet_wrap(~quantiles) +
#         theme_classic() +
#         theme(legend.position="top",
#               aspect.ratio=1,
#               legend.text = element_text(angle=90,size = 15),
#               legend.title = element_text(size = 16,vjust = 0.75),
#               legend.key.size = unit(1.0, 'cm'),
#         )
#     })
#     l_p <- lapply(l_p,\(li) ggsave(tempfile(fileext = ".png"),plot = li,
#                                    width = 5, height = 7))
#     l_p <- lapply(l_p,\(li) image_read(li) %>% image_trim)
#     # base
#     img_final <- image_append(do.call("c",l_p[c("5","95")]))
#     # padronizacao da mediana
#     rect_info <- image_info(img_final)
#     rect_width <- rect_info$width
#     rect_height <- rect_info$height
#     longest_side <- max(rect_width, rect_height)
#     l_p[["50"]] <- image_resize(l_p[["50"]], geometry_size_pixels(longest_side))
#     img_final <- image_append(c(l_p[["50"]],img_final),stack = TRUE)
#     img_final <- image_annotate(
#       img_final,
#       text = vsite,
#       location = "-1000+300",
#       gravity = "north",
#       size=180,
#       boxcolor="lightblue",
#       color = "black")
#     img_final <- image_resize(img_final,"50%")
#     img_final <- image_annotate(
#       img_final,
#       text = vname,
#       location = "-525+15",
#       gravity = "north",
#       size=80,
#       boxcolor="lightgreen",
#       color = "black")
#     vpath <- image_write(img_final,path = tempfile(fileext = ".png"))
#     rm(img_final);gc()
#     return(vpath)
#   }
#   l_img <- lapply(l_path$te,f_porefeito)
#   l_img <- lapply(l_img, image_read) %>% 
#     lapply(.,image_trim)
#   img_return <- image_append(do.call("c",l_img),stack = FALSE)
#   image_write(img_return,
#               path=paste0("figuras/predito_te_sites/",vsite,".png"))
# }
# lapply(v_sitesaud,f_aud)
# 
# 
# 
# ##