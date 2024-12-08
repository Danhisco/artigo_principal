v_figuras <- c("figuras/tabela_box1.png",
               "figuras/GE_dados_disponiveis.png",
               "figuras/figSI_S_N_sizeDA_fillSN.png",
               "figuras/fig2a_tendencia_media.png",
               "figuras/fig2_U16_k.png",
               "figuras/fig3_Uk_ladoKM.png",
               "figuras/fig5_SoE.png",
               "figuras/pedacofigfinal_1alinha.png",
               "figuras/tabsel_fragtotal.png",
               "figuras/diagfinal_fragtotal.png",
               "figuras/figfinal_te_fragtotal.png",
               "figuras/3alinha_fragtotal.png",
               "figuras/tabsel_fragperse.png",
               "figuras/diagfinal_fragperse.png",
               "figuras/figfinal_te_fragperse.png",
               "figuras/3alinha_fragperse.png",
               "figuras/tabsel_areaperse.png",
               "figuras/diagfinal_areaperse.png",
               "figuras/figfinal_te_areaperse.png",
               "figuras/3alinha_areaperse.png",
               "figuras/beta0_figura_final_logOR.png")
library(magick)
#
l_img <- lapply(v_figuras,image_read) %>% 
  lapply(.,image_trim) %>% 
  lapply(.,image_resize,geometry="50%")
names(l_img) <- v_figuras
lapply(names(l_img),\(li){
  image_write(l_img[[li]],path=gsub(".png",".jpeg",li))
})
#
l_jpeg <- list.files(path="figuras",pattern=".jpeg",full.names = TRUE)
l_png <- list.files(path="figuras",pattern=".png",full.names = TRUE)
lapply(l_png,\(li){
  file.rename(from=li,
              to=gsub("figuras","figuras_extras",li))
})
lapply(l_jpeg,\(li){
  img <- image_read(li)
  image_write(img,gsub("jpeg","jpg",li))
})
#
#
#
###############################
lf <- list.files(path="./figuras",pattern="jpeg",full.names = TRUE)
library(magick)
lapply(lf,\(li){
  img <- image_read(li)
  image_write(img, li, format = "jpeg")
})

vpaths <- c("GE_dados_disponiveis",
            "figSI_S_N_sizeDA_fillSN",
            "fig2a_tendencia_media",
            "fig2_U16_k",
            "fig3_Uk_ladoKM",
            "fig5_SoE",
            "pedacofigfinal_1alinha")
vpaths <- paste0("figuras/",vpaths,".jpeg")
lapply(vpaths, \(li){
  img <- image_read(li)
  image_write(img,path=gsub("jpeg","png",li),format="png")
})
