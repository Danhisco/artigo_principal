library(raster)
library(png)
# função para ajustar resolução de um arquivo png
func_png.ajust <- function(file, densidade){ # atualizar para o pacote 'magick'
  system(paste(
    "convert ",file, " -resize ", densidade*A_landscape,"@ ", file,  
    sep = ""
  ))
}
# função para desenhar a área amostral quadrada com N indivíduos
f_area.simulada <- function(matriz, N, coordenada_central=TRUE,index_janela=8/4){
  # determinação da matriz para aplicçaão da rotina: usar a coordenada central do tif ou uma janela com nova coordenada central?
  if(coordenada_central){
    l <- ceiling(dim(matriz)[1]/2)
    c <- ceiling(dim(matriz)[2]/2)
    d <- ceiling(sqrt(N)*index_janela)
    m_temp <- matriz[(l-d):(l+d),(c-d):(c+d)]
  } else{
    m_temp <- matriz
  }
  # Janela de preparo é suficiente para criar a área da parcela?
  if(length(m_temp[m_temp==1]) < N) { 
    stop("habitat insuficiente na janela de observação")
  } else if (length(m_temp[m_temp==1]) == N) { 
    stop("area amostral igual janela de observacao")
  } else { # caso seja adequado, use esse código antigo e misterioso:
    col_cresc <- which(m_temp==m_temp, arr.ind = T)
    col_decre <- col_cresc[dim(col_cresc)[1]:1,]
    row_cresc <- col_cresc[order(col_cresc[,1],decreasing = FALSE),] 
    row_decre <- row_cresc[dim(col_cresc)[1]:1,] 
    ciclo <- (dim(m_temp)[1]-1)/2
    l_mat_index <- list()
    dim_temp <- dim(m_temp)[1]
    for(i in 1:ciclo){
      a1 <- col_cresc[col_cresc[,"col"]==i,]
      b1 <- row_cresc[row_cresc[,"row"]==dim_temp+1-i,]
      c1 <- col_decre[col_decre[,"col"]==dim_temp+1-i,]
      d1 <- row_decre[row_decre[,"row"]==i,]
      l_mat_index[[i]] <- do.call(rbind,list(a1,b1,c1,d1))
    }
    l_mat_index[[(dim_temp+1)/2]] <- col_cresc[col_cresc[,"col"]==(dim_temp+1)/2,]
    mat_ref <- unique(do.call(rbind, l_mat_index))
    length_ref <- length(m_temp[mat_ref][m_temp[mat_ref]==1])
    m_temp[mat_ref][m_temp[mat_ref]==1][(1+length_ref-N):length_ref] <- 2
  }
  # se tiver usado a coordenada central então retorne a paisagem completa
  if(coordenada_central){
    matriz[(l-d):(l+d),(c-d):(c+d)] <- m_temp
    return(matriz)
  } else{ # se tiver usado a janela centrada em outra coordenada, apenas a matrix de trabalho
    return(m_temp)
  }
}
# função que aplica as duas funções anteriores e salva o resultado em .txt
f_landscape_for_MNEE <- function(df,txt_repo = "dados/simulacao//"){
  # df precisa de:
  # @ tif.path: caminho para o arquivo .tif com o mapa completo da paisagem ao redor
  # @ lado_km: númerico que indica a dimensão da extensão esspacial, por exemplo, 16.02 (km) 
  # @ DA: densidade óbservada na parcela 
  # @ Ntotal: número de indivíduos na parcela
  # @ SiteCode: código único de cada parcela
  #
  # descrição:
  # preparo os arquivos .tif para os mapas adequados para rodar MNEE:
  # 1) reduz a paisagem ao redor para 16.02 x 16.02 km2
  # 2) ajuste a resolução do .tif tal que a densidade de pixels é igual a densidade de indivíduos na parcela
  # 3) cria a parcela quadrada no centro da paisagem, trocando 1s (unidades de habitat) por 2s (habitat da parcela)
  # 4) salva em .txt
  #
  # leitura e preparo da paisagem .tif
  m_raster <- as.matrix(raster(df$tif.path))
  i_centro <- nrow(m_raster)/2
  i_ladoKm <- df$lado_km * 1000/60
  m_i <- m_raster[(i_centro+1-i_ladoKm):(i_centro+i_ladoKm),
                  (i_centro+1-i_ladoKm):(i_centro+i_ladoKm)]
  # escrita da paisagem no temporario .png
  png.path <- tempfile(fileext = ".png")
  writePNG(image = m_i,target = png.path)
  # processamento para MNEE
  A_landscape <- (df$lado_km ^ 2) * 100
  func_png.ajust(file=png.path,densidade = df$DA) # função que ajusta a resolução
  m_png <-  png::readPNG(png.path) |> as.matrix()
  m_png[m_png>=0.7] <- 1 # o ajuste da resolução cria, 
  m_png[m_png<0.7] <- 0  # valores entre 0 e 1
  mat_tri <- try(f_area.simulada(matriz = m_png, N = df$Ntotal)) # função que desenha a parcela
  # salvamento no formato .txt
  if(is.matrix(mat_tri)){ 
    try(write.table(x = mat_tri, 
                    file = paste0(txt_repo,df$SiteCode,".txt"),
                    sep = " ", row.names = FALSE, col.names = FALSE))
  } else{print("f_area.simulada: not a matrix")}
  unlink(png.path)
}
