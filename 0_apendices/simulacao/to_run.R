# funções
source("source/dinamica_coalescente.R")
source("source/mapa_p_MNEE.R")
# pacotes
library(doMC)
library(gridExtra)
library(ggplot2)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
## dados comuns
select <- dplyr::select
# df_p
df_p <- read_csv("dados/df_p.csv")
#df_sítios simuláveis em MNEE contemporâneo
v_site <- read_csv("dados/csv/resultados_MN/MNEE/cont/df_resultados.csv") |> 
  pull(SiteCode) |> unique()
# df_sim
df_sim <- read_csv("dados/df_simulacao.csv") |> 
  dplyr::select(-tif.path) |> 
  inner_join(y=mutate(
    data.frame(txt.path = list.files(path="dados/simulacao",pattern=".txt")),
    SiteCode = str_extract(txt.path,".*?(?=.txt)")),
    by="SiteCode") |> 
  filter(SiteCode %in% v_site) |>
  left_join(df_p,by="SiteCode")
df_sim <- cbind(df_sim,land_type = rep(x=c("ideal","non_frag"),each=nrow(df_sim)))
df_simSAD <- inner_join(
  df_sim %>% filter(land_type=="ideal") %>% select(-land_type),
  read_csv("dados/csv/taxaU/df_U.csv") %>% filter(land_type!="ideal")
) %>% select(-Usd)
setwd("dados/simulacao")
# df <- df_sim |> filter(SiteCode == "SPigua1" & land_type == "non_frag")
registerDoMC(3)
df_filter <- data.frame(paths = list.files(path="../csv/SADs_neutras/MNEE_taxaU_contemp",
                                           pattern = ".csv",full.names = TRUE))
df_filter$SiteCode <- aaply(df_filter$paths,1,
                            \(x) str_extract(x,"(?<=\\_\\_)(.*?)(?=.csv)"))
df_filter$land_type <- aaply(df_filter$paths,1,
                            \(x) str_extract(x,"(?<=zado\\/)(.*?)(?=\\_\\_)"))
v_filter1 <- df_filter %>% group_by(SiteCode) %>%
  summarise(nland = n()) %>%
  filter(nland==2) %>%
  pull(SiteCode)
# df_filter %>% count(SiteCode)
df_filter <- data.frame(paths = list.files(path="../csv/SADs_neutras/MNEE_taxaU_non_frag",
                                           pattern = ".csv",full.names = TRUE))
df_filter$SiteCode <- aaply(df_filter$paths,1,
                            \(x) str_extract(x,"(?<=\\_\\_)(.*?)(?=.csv)"))
df_filter$land_type <- aaply(df_filter$paths,1,
                             \(x) str_extract(x,"(?<=zado\\/)(.*?)(?=\\_\\_)"))
df_filter %>% count(SiteCode)
v_filter2 <- df_filter %>% group_by(SiteCode) %>%
  summarise(nland = n()) %>%
  filter(nland==2) %>%
  pull(SiteCode)
setdiff(v_filter2,v_filter1)
v_filter <- unique(c(v_filter1,v_filter2))
# df_bySiteLand <- df_simSAD %>% filter(SiteCode=="SPigua1",land_type=="contemp")
d_ply(filter(df_simSAD,!SiteCode%in%v_filter),
      c("SiteCode","land_type"),f_simMNEE3,.parallel = TRUE)
