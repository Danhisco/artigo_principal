# tabela de seleção com deviance explained
# input: uma lista com os modelos candidatos nomeados
# output um data frame com o output de bbmle::AICctab + mgcv::summary$dev.exp
f_TabSelGAMM <- function(l_md){
  df_aicctab <- AICctab(l_md,weights=TRUE) |> as.data.frame()
  df_aicctab$modelo <- row.names(df_aicctab)
  row.names(df_aicctab) <- NULL
  df_dev.exp <- ldply(l_md,.fun = \(x) summary(x)$dev.expl)
  names(df_dev.exp) <- c("modelo","dev.expl")
  df_dev.exp |> 
    inner_join(y=df_aicctab,"modelo") |> 
    arrange(dAICc) |> 
    select(modelo,dAICc:weight,dev.expl)
}
#
# 

