source("source/nameModel.R")
#
f_refit <- \(md){
  v_mes <- md@optinfo$conv$lme4$messages
  if(is.null(v_mes)){
    return(md)
  }else{
    ss <- getME(md,c("theta","fixef"))
    md <- update(md,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
    l_ms_ss <- list("modelo" = md, "starts" = ss)
    # md <- update(md,control=glmerControl(optCtrl=list(maxfun=2e6)))  
    return(l_ms_ss)
  }
}
f_ts_R2 <- \(l_md){
  df_data <- df_plot
  f_R2 <- \(namemd){
    ss <- l_starts[[namemd]]
    md <- l_mds[[namemd]]
    MuMIn::r.squaredGLMM(md) %>% as.data.frame() %>% apply(.,2,mean)
  }
  v_starts <- sapply(l_md,is.list)
  l_starts <- sapply(l_md, \(l){
    if(is.list(l)){
      return(l$starts)
    }else{
      return(NA)
    }
  })
  l_mds <- lapply(l_md,\(l){
    if(is.list(l)){
      return(l$modelo)
    }else{
      return(l)
    }
  })
  df_ts <- AICctab(l_mds,mnames = names(l_mds),weights=TRUE) %>% as.data.frame() %>% 
    mutate(modelo = names(l_mds))
  row.names(df_ts) <- NULL
  df_r2 <- adply(names(l_md),1,f_R2,.id = "md")
  df_r2$md <- names(l_md)
  inner_join(df_ts,df_r2,by=c("modelo"="md")) %>% relocate(modelo)
}
f_predicaoBootStrap <- \(lme,df_newdat){
  ## Passo 2: crie as função que devem ser calculadas dos modelos a cada simulação
  ## Previstos por efeitos fixos e aleatórios
  f1 <- function(.) predict(., newdata=df_newdat)
  ## Previstos por efeitos fixos (argumento re.form=~0)
  f2 <- function(.) predict(., newdata=df_newdat, re.form=~0)
  ## Os dois bootstraps. Ajuste o argumento ncpus para o numero de cores de seu computador
  b3 <- bootMer(lme, FUN = f1, nsim=1000, parallel="multicore", ncpus=2)
  b4 <- bootMer(lme, FUN = f2, nsim=1000, parallel="multicore", ncpus=2)
  ## calcula as médias e intervalos de confiança quantílicos para cada combinação de preditoras
  ## no novo conjunto de dados
  df_newdat$mean <- apply(b3$t,2,mean)
  df_newdat$IC.low <- apply(b3$t,2,quantile, 0.025)
  df_newdat$IC.upp <- apply(b3$t,2,quantile, 0.975)
  df_newdat$mean.fixed <- apply(b4$t,2,mean)
  df_newdat$IC.low.fixed <- apply(b4$t,2,quantile, 0.025)
  df_newdat$IC.upp.fixed <- apply(b4$t,2,quantile, 0.975)
  return(df_newdat)
}
f_ts_converg.md <- \(l_md_i){
  f_warning <- function(x) x@optinfo$conv$lme4$messages %>% 
    paste(.,collapse = " , ") %>% sub("$","NA",.)
  f_pseudoR2 <- \(md){
    MuMIn::r.squaredGLMM(md) %>% as.data.frame() %>% apply(.,2,mean)
  }
  df_warnings <- ldply(l_md_i,f_warning,.id="glmer")
  v_convmd <- df_warnings %>% filter(V1=="NA") %>% pull("glmer") %>% as.character()
  l_md_i <- l_md_i[v_convmd]
  df_ts <- AICctab(l_md_i,weights=TRUE,mnames = names(l_md_i)) %>% 
    {if(length(v_convmd)==1) as.data.frame(.,row.names = v_convmd) else as.data.frame(.)} %>% 
    mutate(.,md=rownames(.))
  row.names(df_ts) <- NULL
  df_r2 <- ldply(l_md_i,f_pseudoR2,.id = "md")
  inner_join(df_ts,df_r2,by="md") %>% relocate(md)
}
f_loadll <- \(l_path){
  v_name <- c()
  l_md <- list()
  for(i in 1:length(l_path)){
    v_name[i] <- load(l_path[[i]],verbose = TRUE)
    l_md[[i]] <- get(v_name[i])
  }
  names(l_md) <- v_name
  return(l_md)
}
f_glmm <- \(lf,df_data=df_md,Rdata_path="./dados/Rdata/glmm_"){
  df_code <- data.frame(
    nome=c("nCong"),
    funcao=c("function(f) glmer(as.formula(f),family='binomial',data=df_data,control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=2e9)))"),
    f_resp = c("cbind(nCong,100-nCong) ~"))
  f_a_ply <- \(df){
    env <- environment()
    f_md <- eval(parse(text = df$funcao),envir=env)
    f <- gsub("~",df$f_resp,lf)
    l_md <- lapply(f,f_md)
    names(l_md) <- names(lf)
    name_l <- df$nome
    name_l <- paste0("l_md_",name_l)
    assign(name_l,l_md)
    save(list=name_l,file=paste0(Rdata_path,name_l,".Rdata"))
  }
  a_ply(df_code,1,f_a_ply)
}
f_evalhere <- \(s) eval(parse(text = s),envir=env)
ll_ggpng <- \(l_paths){
  if(!is.list(l_paths)) l_paths <- as.list(l_paths)
  lapply(l_paths,png::readPNG) %>% 
    lapply(.,\(p) ggplot() +
             annotation_custom(grid::rasterGrob(p), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
             theme_void())
}
f_gt <- \(l_d,vtitle=NULL,v_footnote=NULL){
  # vtitle <- c("Códigos atuais para FluA ('PCR_FLUASU')",
  #             "Códigos atuais para FluB ('PCR_FLUBLI')",
  #             "Códigos dicionário 09-12 ('*_ETIOL')",
  #             "Códigos dicionário 13-18 ('RES_FLUASU')")
  names(vtitle) <- names(l_d)
  if(!is.null(v_footnote)) names(v_footnote) <- names(l_d)
  l_d <- lapply(names(l_d),\(s) l_d[[s]] %>% 
                  gt::gt() %>% 
                  gt::tab_header(title=vtitle[s]) %>% 
                  {if(s=="preditoras")  gt::tab_footnote(.,
                    footnote = v_footnote[s],
                    locations = cells_body(
                      columns = nome,
                      rows = 1)
                  ) else .} %>% 
                  {if(s=="contrastes")  gt::tab_footnote(.,
                                                     footnote = v_footnote[s],
                                                     locations = cells_body(
                                                       columns = `interpretação`,
                                                       rows = 3)
                  ) else .})
  names(l_d) <- names(vtitle)
  return(l_d)
} # f especifica para esse chunk: correspondência entre os nomes dos df e os títulos das tabelas
f_gt_to_ggplot <- \(l_d){
  l_temp <- replicate(length(l_d),tempfile(fileext = '.png'))
  a_ply(1:length(l_d),1,\(i) gtsave(data = l_d[[i]],filename = l_temp[i]))
  l_p <- ll_ggpng(l_temp)
  names(l_p) <- names(l_d)
  l_ply(l_temp,file.remove)
  return(l_p)
}
