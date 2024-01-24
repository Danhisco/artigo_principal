# funções
library(twosamples)
library(sads)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
source("source/2samples_testes.R")
# dados
df_SADrep2 <- mutate(data.frame(
  SADrep.path = list.files(path = "dados/csv/SADs_neutras/MNEE_taxaU_idealizado",
                           pattern = ".csv",full.names = T,recursive = T)),
  SiteCode = str_extract(SADrep.path,"(?<=\\_\\_).*?(?=.csv)"),
  land_type = str_extract(SADrep.path,"(?<=zado\\/).*?(?=\\_\\_)"),
  SADobs.path = paste0("dados/csv/SADs_observadas/",SiteCode,".csv"))
# rotina
registerDoMC(3)
df_repMNEE <- adply(df_SADrep2,1,f_resultsMN,.parallel = TRUE)
write_csv(df_repMNEE, file = "dados/csv/resultados_MN/MNEE/df_repMNEE2.csv")