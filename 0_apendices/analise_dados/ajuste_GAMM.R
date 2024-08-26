



## modelos com ?????
l_md_logOR <- readRDS(file="./dados/Rdata/l_md_logOR_cgi_c_perturbacao.rds")
## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
## tabela de selação dos modelos 
df_tabsel_logOR <- ldply(l_md_logOR,f_TabSelGAMM,.id="pair")
write_csv(df_tabsel_logOR,file="dados/csv/df_tabelaSelecao_logOR_cgi_cperturbacao.csv")
df_tabelaSelecao <- read_csv(file="dados/csv/df_tabelaSelecao_logOR_cgi_cperturbacao.csv")
# write_csv(df_tabelaSelecao,
#           file=paste0(v_path,"tabelas/df_tabelaSelecaologOR_cgi.csv"))

f_SML.gamm <- \(df){md <- l_md_logOR[[df$pair[1]]][[df$modelo[[1]]]]}
l_md_1alike <- dlply(df_tabelaSelecao,"pair",f_SML.gamm)
saveRDS(l_md_1alike,file="dados/Rdata/l_md_logOR_itself_1alike_cperturbacao.rds")


##########################################################
############# Qual o modelo cheio suficiente? ############
##########################################################
# objectos comuns
## dados
df_md <- df_logOR %>% filter(forest_succession!="capoeira")
## modelos usados
l_md_1alike <- readRDS(file="dados/Rdata/l_md_logOR_itself_1alike_cperturbacao.rds")
## função de atualização do modelo cheio mais plausível
f_gam_update <- \(gamm){
  # setup
  f_ref <- formula(gamm)
  original_data <- gamm$model
  f_gam <- \(fi){
    gam(formula=fi,data=original_data,method = "REML")
  }
  # formula 
  l_md <- list()
  l_md$`-Uefeito*perturb` <- update.formula(
    f_ref, . ~ . - s(Uefeito, by = forest_succession, bs = "cr") + s(Uefeito, bs = "cr"))
  # ajuste das formulas
  l_md <- lapply(l_md,f_gam)
  # inclusão do modelo cheio
  l_md$cheio <- gamm
  return(l_md)
}
l_md_cheios <-lapply(l_md_1alike,f_gam_update) 
saveRDS(l_md_cheios,"dados/Rdata/l_md_cheios_cperturbacao.rds")
##########################################################
######### Qual o melhor conjunto de covariáveis? #########
##########################################################
# objetos necessários
## modelos cheios com perturbação
l_md_logOR_aud <- readRDS(file="dados/Rdata/l_md_cheios_cperturbacao.rds")
## tabela de seleção
if(!file.exists("dados/csv/df_tabelaSelecao_logOR_cpert.csv")){
  df_tabsel_logOR_aud <- ldply(l_md_logOR_aud,f_TabSelGAMM,.id="pair")
  write_csv(df_tabsel_logOR_aud,
            file="dados/csv/df_tabelaSelecao_logOR_cpert.csv")
}else{
  df_tabsel_logOR_aud <- 
    read_csv(file="dados/csv/df_tabelaSelecao_logOR_cpert.csv")
}

