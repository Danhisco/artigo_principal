source("source/nameModel.R")
f_evalhere <- \(s) eval(parse(text = s),envir=env)
ll_ggpng <- \(l_paths){
  if(!is.list(l_paths)) l_paths <- as.list(l_paths)
  lapply(l_paths,png::readPNG) %>% 
    lapply(.,\(p) ggplot() +
             annotation_custom(grid::rasterGrob(p), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
             theme_void())
}
f_gt <- \(l_d,vtitle,v_footnote){
  # vtitle <- c("Códigos atuais para FluA ('PCR_FLUASU')",
  #             "Códigos atuais para FluB ('PCR_FLUBLI')",
  #             "Códigos dicionário 09-12 ('*_ETIOL')",
  #             "Códigos dicionário 13-18 ('RES_FLUASU')")
  names(vtitle) <- names(l_d)
  names(v_footnote) <- names(l_d)
  l_d <- lapply(names(l_d),\(s) l_d[[s]] %>% 
                  gt() %>% 
                  tab_header(title=vtitle[s]) %>% 
                  {if(s=="preditoras")  tab_footnote(.,
                    footnote = v_footnote[s],
                    locations = cells_body(
                      columns = nome,
                      rows = 1)
                  ) else .} %>% 
                  {if(s=="contrastes")  tab_footnote(.,
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
