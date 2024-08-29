# pacotes
# library(grid)
# library(gtable)
# library(ggpubr)
# library(gt)
# library(flextable)
# library(dagitty)
# library(ggdag)
library(gratia)
# library(sads)
library(doMC)
# library(metR)
# library(kableExtra)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
library(readr)
# library(purrr)
library(stringr)
library(tidyr)
# library(MuMIn)
# library(AICcmodavg)
# library(insight)
library(bbmle)
library(DHARMa)
library(mgcv)
library(lme4)
# library(data.table)
library(plyr)
library(dplyr)
## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
paste0(v_path,"rds/l_md_3aperg_quais_cov_")
paths <- list.files(path=paste0(v_path,"rds"),
                    pattern="l_md_3aperg",full.names = T)
# 
f_df_pred <- \(vpath){
  l_md <- readRDS(vpath)
  gamm <- l_md[[1]]
  
}


f_lapply <- \(ldf){
  f_PI <- \(gamm){
    lpi <- list(
      with_random = list("apenas fixo" = FALSE,
                         "aleatório e fixo" = TRUE),
      to_exclude = list("apenas fixo"= c("s(Uefeito,SiteCode)",
                                         "s(SiteCode)",
                                         "s(lat,long)",
                                         "s(data_year)"),
                        "aleatório e fixo" = c("s(lat,long)",
                                               "s(data_year)"))
    )
    adply(c("apenas fixo","aleatório e fixo"),1,\(i){
      f_calcPI2(gamm=gamm,
                with_random = lpi$with_random[[i]],
                to_exclude = lpi$to_exclude[[i]])
    },.id = NULL)
  }
  registerDoMC(3)
  ldply(ldf,f_PI,.id="gamm",.parallel = TRUE)
}