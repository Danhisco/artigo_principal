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
## predição a posteriori 
#
### te for every model
l_path <- list()
l_path$te <-  paste0("rds/l_dfpred_",c("areaperse","fragperse","fragtotal"),".rds")
l_path$U <- "rds/l_dfpred_areaperse_Ugs.rds"
df_tabsel <- read_csv(paste0(v_path,"rds/df_tabsel_geral.csv")) %>% 
  filter(dAICc==0)
# veffect <- l_path$te[1]
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
  
    
}
