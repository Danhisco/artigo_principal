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
## predição a posteriori 
#
### te for every model
l_path <- list()
l_path$te <-  paste0("rds/l_dfpred_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$U <- "rds/l_dfpred_areaperse_Ugs.rds"
df_tabsel <- read_csv("rds/df_tabsel_geral.csv") %>% 
  filter(dAICc==0)
# veffect <- l_path$te[3]
f_plot_te <- \(veffect,
               pattern_extract="(?<=l_dfpred_)(.*?)(?=\\.rds)"){
  library(metR)
  vname <- str_extract(veffect,pattern_extract) %>% 
    gsub("areaperse","Área per se",.) %>% 
    gsub("fragperse","Frag. per se",.) %>% 
    gsub("fragtotal","Frag. total",.)
  dfpred <- readRDS(paste0(v_path,veffect))$`apenas fixo`
  dfpred <- pivot_longer(dfpred,starts_with("Q_"),
                         names_to="quantiles",
                         values_to="logOR") %>% 
    mutate(quantiles=factor(100 * as.numeric(gsub("Q_","",quantiles)),
                            levels=c(5,50,95)) ) %>% 
    rename(k = k_cont)
  #
  l_p <- dlply(dfpred,"quantiles",\(dfi){
    ggplot(dfi,aes(x=k,y=Uefeito,z=logOR)) +
      geom_contour_fill() +
      geom_contour(color = "black") +
      geom_text_contour() +
      scale_fill_viridis_c(name = "logOR",
                           option = "magma") +
      labs(y="logU/U") +
      facet_wrap(~quantiles) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_classic() +
      theme(legend.position="none",
            aspect.ratio=1)
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
# lapply(l_path$te,f_plot_te)
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