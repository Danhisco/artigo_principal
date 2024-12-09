library(gratia)
library(doMC)
library(gridExtra)
library(ggplot2)
library(cowplot)
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
  library(magick)
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
  )# %>% image_trim
}
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
    df_slopes <- ddply(dfi,"SiteCode",\(dfii){
      mdl <- as.data.frame(lm(logOR~preditor,data=dfii)$coefficients)
      mdl$coef = rownames(mdl) %>% gsub("\\(","",.) %>% gsub("\\)","",.)
      rownames(mdl) <- NULL
      names(mdl)[1] <- "valor"
      pivot_wider(mdl,names_from="coef",values_from="valor")
    }) %>% pivot_longer(-SiteCode) %>% 
      mutate(vline = ifelse(name=="preditor",1,0),
             name = gsub("Intercept","intercepto",name) %>% 
               gsub("preditor","inclinação",.),
             name=factor(name,levels=c("intercepto","inclinação")))
    l_p <- list()
    l_p$coef <- ggplot(df_slopes,aes(x=value)) +
      geom_histogram(aes(y=after_stat(density)),bins=15,fill="skyblue",color="black") +
      geom_density(color="darkgreen",linewidth=1) +
      geom_vline(aes(xintercept=vline),color="red",linetype=2,linewidth=1) +
      scale_y_continuous("densidade",expand=c(0,0)) +
      xlab("estimativa") +
      facet_wrap(~name,ncol=1,scales="free") +
      theme(strip.text = element_text(size=13,margin = margin()),
            plot.margin = margin(0,0,0,0))
    l_p$predobs <- dfi %>% 
      mutate(quantile=paste0("Quantil: ",quantile,"%")) %>% 
      ggplot(aes(x=logOR,y=preditor)) +
      geom_smooth(aes(group=SiteCode),
                  method="lm",
                  se=FALSE,
                  color="#4682B4",
                  alpha=0.5) +
      geom_hex(bins = 50,alpha=0.5) + 
      geom_abline(slope=1,intercept=0,color="black",linewidth=1.2,linetype=2) +
      scale_fill_gradient("contagem",low = "yellow", high = "red", na.value = NA) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      labs(x="logOR observado",y="logOR predito") +
      facet_wrap(~quantile) +
      theme(legend.position = c(0.90,0.25),
            strip.text = element_text(size=13,margin = margin()),
            plot.title = element_text(size=10,
                                      hjust=0.5,
                                      face="bold"),
            axis.title.x = element_text(margin=margin(t=0,b=0,l=0,r=0)),
            axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=0))
            )
    if(dfi$quantile[1]=="50"){
      p <- ggdraw() +
        draw_plot(l_p[["predobs"]]) +
        draw_plot(l_p[["coef"]],
                  height = 0.4, width = 0.25,
                  x=0.08, y=0.55)
      vpath <- f_imagefunc(p)
      return(vpath)
    }else{
      p <- l_p$predobs +
        theme(plot.tile=element_text(size=))
      vpath <- f_imagefunc(l_p$predobs)
      return(vpath)
    }
  }
  vpath <- fggplot(filter(dfp,quantile=="50"))
  # lp <- dlply(dfp,"quantile",fggplot)
  # limg <- lapply(lp,\(li) image_trim(image_read(li)))
  # base_img <- image_append(do.call("c",limg[c("5","95")]),stack=FALSE)
  # base_img <- f_resize_2rectangle(ref_img = limg[["50"]],
  #                                 toresize_img = base_img,
  #                                 ref_side = "width",
  #                                 tore_side = "width")
  # img_final <- image_append(c(limg[["50"]],base_img),stack = TRUE) %>% 
  #   image_resize("50%")
  img <- image_read(vpath) %>% 
    image_trim() %>% 
    image_resize("50%")
  return(img)
}
l_img_diag <- lapply(l_df,f_plot)
l_img_diag <- lapply(names(l_img_diag),\(li){
  img <- l_img_diag[[li]]
  image_title(imgobj = img, 
              vtitle = li,
              vsize = 50,
              vheight = 90) %>% 
    image_trim
})
names(l_img_diag) <- paste0(
  "figuras/diagfinal_",
  case_when(
    grepl("Área",names(l_df)) ~ "areaperse",
    grepl("Frag. per se",names(l_df)) ~ "fragperse",
    grepl("total",names(l_df)) ~ "fragtotal",
  ),
  ".jpeg"
)
lapply(names(l_img_diag),\(li) image_write(l_img_diag[[li]],path=li))
#
# ljpg <- list.files("figuras",pattern="jpg",full.names = TRUE)
# file.remove(ljpg)
#

# f_dfplot <- \(dff,dataset){
#   vcols <- c("logOR","k_cont")
#   dff %>% 
#     select(-any_of(vcols)) %>% 
#     pivot_longer(starts_with("Q_"),names_to="quantile",
#                  values_to=paste0("logOR_",dataset)) %>% 
#     mutate(quantile=gsub("Q_","",quantile) %>% 
#              as.numeric(),
#            quantile=as.character(quantile*100))
# }
# df_te <- lapply(l_df,"[[","fixo e aleat")[["Área per se"]] %>% f_dfplot(dataset="te")
# df_s <- readRDS(paste0(v_path,"rds/l_dfpred_areaperse_Ugs.rds"))[["fixo e aleat"]] %>% 
#   f_dfplot(dataset="s")
# df_plot <- inner_join(df_te,df_s) %>% 
#   inner_join(lapply(l_df,"[[","fixo e aleat")[["Área per se"]] %>% 
#                select(logOR,Uefeito,SiteCode))
# f_ggplot <- \(dfi){
#   l_p <- list()
#   l_p$predpred <- dfi %>% 
#     ggplot(aes(x=logOR_te,y=logOR_s,color=quantile)) +
#     geom_abline(intercept = 0,slope=1) +
#     geom_point() +
#     geom_smooth(method="lm",se=FALSE) +
#     scale_color_manual("%",
#                        values=c("50"="darkgreen",
#                                 "5"="darkred",
#                                 "95"="darkred")) +
#     labs(title="predito X predito") +
#     theme(legend.position = "top",
#           aspect.ratio = 1) +
#     facet_wrap(~SiteCode) 
#   l_p$obs_uefeito_predicoes <- dfi %>% 
#     # preparação dos dados
#     pivot_longer(matches("_te|_s"),
#                  names_to = "hgam",
#                  values_to = "logOR_pred") %>% 
#     mutate(hgam=gsub("logOR_","",hgam)) %>% 
#     ggplot(aes(x=logOR,y=logOR_pred,color=quantile)) +
#     geom_abline(intercept = 0,slope=1) +
#     geom_point() +
#     geom_smooth(method="lm",se=FALSE) +
#     scale_color_manual("%",
#                        values=c("50"="darkgreen",
#                                 "5"="darkred",
#                                 "95"="darkred")) +
#     labs(x="logOR observado",y="logOR predito",
#          title=dfi$SiteCode[1]) +
#     theme(legend.position = "top",
#           aspect.ratio = 1) +
#     facet_wrap(~hgam) 
#   l_p <- lapply(l_p,\(li) ggsave(tempfile(fileext = ".png"),plot = li))
#   l_p <- lapply(l_p,\(li) image_read(li) %>% image_trim)
#   img_final <- image_append(do.call("c",l_p),stack = FALSE)
#   image_write(
#     img_final,
#     path=paste0(v_path,
#                 "figuras/PIeOBS_sites/compara_areaperse/",
#                 dfi$SiteCode[1],
#                 ".png")
#   )
# }
# lapply(levels(df_s$SiteCode),\(i){
#   filter(df_plot,SiteCode==i) %>% f_ggplot()
# })

##### combinação das figuras de diagnóstico
l_img <- list.files("figuras",pattern = "diagfinal",full.names = TRUE)
l_img <- lapply(l_img,image_read)
names(l_img) <- list.files("figuras",pattern = "diagfinal",full.names = TRUE)
img_final <- image_append(do.call("c",l_img[3:1]),stack = FALSE) %>% 
  image_resize("50%")
image_write(img_final,"figuras/diagfinal_3efeitos.jpeg",format="jpeg")