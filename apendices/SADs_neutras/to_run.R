# funções
library(twosamples)
library(doMC)
library(sads)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
source("source/2samples_testes.R")
# rotinas
if(!file.exists("dados/csv/resultados_MN/MNEE/df_repMNEE2.csv")){
  df_SADrep2 <- mutate(data.frame(
    SADrep.path = list.files(path = "dados/csv/SADs_neutras/MNEE_taxaU_idealizado",
                             pattern = ".csv",full.names = T,recursive = T)),
    SiteCode = str_extract(SADrep.path,"(?<=\\_\\_).*?(?=.csv)"),
    land_type = str_extract(SADrep.path,"(?<=zado\\/).*?(?=\\_\\_)"),
    SADobs.path = paste0("dados/csv/SADs_observadas/",SiteCode,".csv"))
  registerDoMC(3)
  df_repMNEE <- adply(df_SADrep2,1,f_resultsMN,.parallel = TRUE)
  write_csv(df_repMNEE, file = "dados/csv/resultados_MN/MNEE/df_repMNEE2.csv")  
}else if(!file.exists("./dados/csv/resultados_MN/MNEE/df_repMNEE2.csv")){
  registerDoMC(3)
  df_repMNEE <- read_csv("./dados/csv/resultados_MN/MNEE/df_repMNEE2.csv")
  df_repRes <- ddply(df_repMNEE,c("SiteCode","k","land_type"),f_summarise_SAD_MNEE,.parallel = TRUE)
  write_csv(df_repRes, file = "dados/csv/resultados_MN/MNEE/df_repRes2.csv") 
}else if(file.exists("./dados/csv/resultados_MN/MNEE/df_repMNEE2.csv") &
         !file.exists("./dados/csv/resultados_MN/MNEE/df_resultados2.csv")){
  df_resultados <- read_csv("dados/csv/resultados_MN/df_resultados.csv")
  registerDoMC(3)
  df_MNEE_Uidealizado <- read.csv("./dados/csv/resultados_MN/MNEE/df_repRes2.csv") %>% 
    rbind.fill(df_resultados %>% filter(land_type=="ideal") %>% select(SiteCode,k:Smax)) %>% 
    arrange(SiteCode,k) %>%
    ddply(.,c("SiteCode","k"),f_logOR_land_type,.parallel = TRUE)
  write_csv(df_MNEE_Uidealizado,
            file="./dados/csv/resultados_MN/MNEE/df_resultados2.csv")
}else{
  warning("Tem trabalho para fazer? Se sim, remover os .csv")
}
