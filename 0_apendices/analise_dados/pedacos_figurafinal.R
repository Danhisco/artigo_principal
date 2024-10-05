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
# objetos
## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
## predição a posteriori 
#
### te for every model
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
l_path <- list()
l_path$te <-  paste0("rds/l_dfpred_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$U <- "rds/l_dfpred_areaperse_Ugs.rds"
df_tabsel <- read_csv(paste0(v_path,"rds/df_tabsel_geral.csv")) %>% 
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
  image_info <- image_info(img_final)
  image_width <- image_info$width
  white_canvas <- image_blank(width = image_width, height = 150, color = "white")
  img_final <- image_append(c(white_canvas, img_final), stack = TRUE)
  img_final <- image_annotate(
    img_final,
    text = vname,
    size = 100,
    gravity = "north",
    color = "black",
    location = "+0+20"
  ) %>% image_trim
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
