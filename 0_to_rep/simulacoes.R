source("source/dinamica_coalescente.R")
source("source/mapa_p_MNEE.R")
source("source/2samples_testes.R")
# pacotes #
library(twosamples)
library(doMC)
library(gridExtra)
library(ggplot2)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
# objetos comums #
# df_p
df_p <- read_csv("dados/df_p.csv")
#df_sítios simuláveis em MNEE contemporâneo
df_coord <- read_csv(file = "dados/df_dados_disponiveis.csv") %>% 
  mutate(lat = ifelse(is.na(lat_correct),lat,lat_correct),
         long = ifelse(is.na(long_correct),long,long_correct),
         Sitecode = factor(SiteCode)) %>% 
  select(SiteCode,lat,long, forest_succession) %>% 
  filter(forest_succession!="capoeira")

v_site <- read_csv(file="dados/csv_SoE/df_congruencia_simulacao.csv") %>% 
  pull(SiteCode) |> unique()
v_site <- intersect(df_coord$SiteCode,v_site)
# df_sim
df_sim <- read_csv("dados/df_simulacao.csv") |> 
  dplyr::select(-tif.path) |> 
  inner_join(y=mutate(
    data.frame(txt.path = list.files(path="dados/simulacao",pattern=".txt")),
    SiteCode = str_extract(txt.path,".*?(?=.txt)")),
    by="SiteCode") |> 
  filter(SiteCode %in% v_site) |>
  left_join(df_p,by="SiteCode") %>% 
  mutate(k=round(k,2))
# df com o efeito de escalar
vk_sim <- 0.33333
vk_sim <- 0.001 
df_tabsoe <- readRDS(file="1_to_compile_dissertacao_EM_USO/09_SI/RDS/df_tabsoe.rds") %>%
  filter(limiar==0.75) %>% 
  mutate(Li=as.numeric(as.character(Li)),
         SoE=ifelse(k<=0.30,4,2))
# df com as informações para simulação
df_sim <- inner_join(
  df_sim,
  select(df_tabsoe,-c(k_factor,EoS))
) %>% 
  filter(k>=vk_sim)

#-----------------------------------------------------------------------
#-------------------------------  SIMULAÇÃO ----------------------------
#-----------------------------------------------------------------------

setwd("dados/simulacao")
# dfi <- df_sim %>% filter(SiteCode=="SPigua1")
formals(f_simMNEE)$general_path <- ".."
formals(f_simMNEE)$pathU <- "/csv_SoE/taxaU/p_consolidar/"
formals(f_simMNEE)$pathSAD <- "/csv_SoE/SADs_neutras/p_consolidar/"
registerDoMC(12)
d_ply(df_sim,"SiteCode",f_simMNEE,.parallel = TRUE)


#-----------------------------------------------------------------------
#-----------------------------  COMPARAÇÃO SADS-------------------------
#-----------------------------------------------------------------------


setwd("../..")
df_SADrep <- mutate(
  data.frame(
    SADrep.path = list.files(path = "dados/csv_SoE/SADs_neutras/MNEE",
                             pattern = ".csv",
                             full.names = T,
                             recursive = T)
  ),
  SiteCode = str_extract(SADrep.path,
                         "(?<=[cont|frag|ideal]\\/).*?(?=_So)"),
  land_type = str_extract(SADrep.path,
                          "(?<=MNEE\\/).*?(?=\\/)"),
  SoE = str_extract(gsub("dados/csv_SoE","",SADrep.path),
                    "(?<=SoE).*?(?=km)"),
  SADobs.path = paste0("dados/csv_SoE/SADs_observadas/",SiteCode,".csv")
)
#
# função que une os data frames de diferentes escaslas em um só e aplica o teste KS
f_bySite_e_Land <- \(dfi,fname = "f_resultsMN",colname="SADrep.path"){
  #
  l_csv <- lapply(dfi[[colname]],read_csv)
  names(l_csv) <- dfi$SoE
  #
  df_sad <- do.call("rbind",l_csv)
  temppath <- tempfile(fileext = ".csv")
  write_csv(df_sad,temppath)
  #
  dfi <- select(dfi,-SoE) %>% 
    mutate(SADrep.path=temppath) %>% 
    distinct()
  f_dfi <- get(fname)
  dfr <- f_dfi(dfi)
  #
  file.remove(temppath)
  #
  return(dfr)
}
#registerDoMC(2)
df_KSrep <- ddply(df_SADrep,
                  c("SiteCode","land_type"),
                  f_bySite_e_Land,
                  .parallel = FALSE)
write_csv(df_KSrep,file="dados/csv_SoE/df_KSrep.csv")

#--------------------------------------------------------------
#------------- Sumário dos resultados do teste KS -------------
#--------------------------------------------------------------

f_summarise_SAD_MNEE <- \(df){
  #@ df: df por site, k, e  land_type
  #@ e.g. ddply(.,c("SiteCode","k","land_type"))
  cbind(df[1,c("SiteCode","k","land_type")],with(df,data.frame(
    nCongKS = sum(p.KS>0.05),
    Smed = mean(S),
    Ssd = sd(S),
    Smin = min(S),
    Smax = max(S)))
  )
}
df_KSrep <- read_csv(file="dados/csv_SoE/df_KSrep.csv")
registerDoMC(3)
df_ad <- ddply(df_KSrep,c("SiteCode","k","land_type"),
               f_summarise_SAD_MNEE,
               .parallel = TRUE)
write_csv(df_ad,file="dados/csv_SoE/df_congruencia_simulacao.csv")


##----------------------------------------------------------
#------------------- Contraste logU/U ----------------------
##----------------------------------------------------------

f_contraste_Umed <- \(dfUrep,
                      path_U="dados/csv/taxaU/df_U.csv",
                      path_land_effect="dados/csv/taxaU/df_contrastes.csv"){
  df_U <- dfUrep |> select(SiteCode:k)
  df_U$Umed <- apply(select(df_Urep,-c(SiteCode:d)),1,mean)
  df_U$Usd <- apply(select(df_Urep,-c(SiteCode:d)),1,sd)
  df_contrastes <- df_U |>
    distinct() %>% 
    pivot_wider(names_from = "land_type", 
                values_from = c("Umed","Usd")) |> 
    mutate(FLF = log( Umed_ideal / Umed_cont),
           FPS = log( Umed_non_frag / Umed_cont),
           LFC = log( Umed_ideal / Umed_non_frag)) |> 
    select(-starts_with("U"))
  write_csv(df_contrastes,file=path_land_effect)
  write_csv(df_U,file=path_U)
}
# dados
df_Urep <- mutate(
  data.frame(
    Urep.path = list.files(path = "dados/csv_SoE/taxaU/MNEE",
                           pattern = ".csv",
                           full.names = T,
                           recursive = T)
  ),
  SiteCode = str_extract(Urep.path,
                         "(?<=[cont|frag|ideal]\\/).*?(?=_So)"),
  land_type = str_extract(Urep.path,
                          "(?<=MNEE\\/).*?(?=\\/)"),
  SoE = str_extract(gsub("dados/csv_SoE","",Urep.path),
                    "(?<=SoE).*?(?=km)")
)
# função para unir as bases de dados de diferentes SoE
f_bySite_e_Land <- \(dfi,
                     colname="Urep.path"){
  l_csv <- lapply(dfi[[colname]],read_csv)
  do.call("rbind",l_csv)
}
registerDoMC(2)
df_Urep <- ddply(df_Urep,c("SiteCode","land_type"),
                 f_bySite_e_Land,
                 .parallel = TRUE)
df_Urep <- relocate(df_Urep,SiteCode)
# sumarização
f_contraste_Umed(df_Urep,
                 path_U = "dados/csv_SoE/taxaU/df_U.csv",
                 path_land_effect <- "dados/csv_SoE/taxaU/df_contrastes.csv")

#--------------------------------------------------------------------
#------------------- Sumário paisagens hipotéticas ------------------
#--------------------------------------------------------------------

source("0_apendices/analise_dados/sumario_paisagens_hipoteticas.R")