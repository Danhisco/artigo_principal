library(doMC)
library(raster)
library(readr)
library(purrr)
library(stringr) 
library(tidyr)
library(plyr)
library(dplyr)
setwd("Documentos/mestrado_Ecologia/artigo_principal")
# funções
source("source/dinamica_coalescente.R")
source("source/SoE_MNEE.R")
#
## dados gerais
df_sim <- read_csv("dados/df_simulacao.csv")
# df_sorteio_sitios
df_sorteio_sitios <- read_csv(file="dados/csv/df_sorteio_sitios.csv") |> 
  select(-d) |> 
  inner_join(x=select(df_sim,SiteCode,d,k),by="SiteCode")
df_sorteio_sitios$txt.file <- gsub(".txt","_null.txt",df_sorteio_sitios$txt.file)
#
## rotina 
setwd("dados/simulacao")
df_sitios_2lote <- data.frame(path.file = list.files("../csv/SoE",".csv",full.names = T)) |>
  mutate(ctime = as.Date(file.info(path.file)$ctime),
         SiteCode = str_extract(path.file,"(?<=E\\/).*?(?=\\_k)"),
         k = as.numeric(str_extract(path.file,"(?<=\\_k).*?(?=\\.csv)")),
         label = paste(SiteCode, k/100)) |>
  filter(ctime=="2023-02-01")
df_se2 <- df_sorteio_sitios |> 
  mutate(label = paste(SiteCode, round(k,2))) |>
  inner_join(select(df_sitios_2lote,label),by="label") |> 
  select(-label)
registerDoMC(2)
a_ply(df_se2,1,f_SoE_MNEE_null,.parallel = FALSE)
