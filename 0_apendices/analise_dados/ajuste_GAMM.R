library(readr)
library(stringr)
library(tidyr)
library(mgcv)
library(plyr)
library(dplyr)
source("source/nameModel.R")
source("source/GAMMtools.R")
load("./dados/Rdata/md_prcong.Rdata")
df_tow <- f_calcPI(gamm=md_prcong,
                   ad = "ad4",
                   newdata_path="./dados/csv/df_newdata_ad4.csv",
                   link_scale = TRUE,
                   n.posteriori_samples = 1000)
write_csv(df_tow,file="./dados/csv/df_newpred_ad4.csv")
