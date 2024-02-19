# taela de selação com R2mc p glmm
#@ aceita-se lista em q cada elemento é um modelo
#
## tabela de seleção com R2
f_tselR2 <- \(l_md){
  # tab sel
  df_return <- AICctab(l_md,weights=TRUE) %>% as.data.frame 
  df_return$formula <- df_return %>% row.names()
  row.names(df_return) <- NULL
  # R2 m e c
  df_R2 <- ldply(l_md,\(li){
    suppressWarnings(MuMIn::r.squaredGLMM(li)) %>% 
      as.data.frame() %>% apply(.,2,mean)
  },.id = "formula")
  # return
  inner_join(df_return,df_R2,"formula") %>% 
    relocate(formula)
}
# residuals DHARMa
f_resDharma <- \(md,nrep=1000,png_fullpath,...){
  resDharma <- simulateResiduals(md,n=nrep)
  if(!missing(...)){
    png(filename = png_fullpath,...)
    plot(resDharma)
    dev.off()
  }else{
    print("quais os par. gráf. para o png?")
  }
}