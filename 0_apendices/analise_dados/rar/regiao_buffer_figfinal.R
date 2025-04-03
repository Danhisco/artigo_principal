##### Construção da região de buffer da predição a posteriori
library(gridExtra)
library(ggplot2)
theme_set(theme_bw())
library(readr)
library(stringr)
library(tidyr)
library(tidyselect)
library(bbmle)
library(DHARMa)
library(mgcv)
library(plyr)
library(dplyr)
library(mgcv)
# objetos comuns
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
df_md <- read_csv(file="dados/csv_SoE/df_logOR.csv")
#
#############################################################
###### máximos e mínimos dos contrastes em função de k ######
#############################################################
df_buf <- ddply(df_md,c("contraste","k"),\(dfi){
  v_range <- range(dfi$Uefeito)
  v_quant <- quantile(dfi$Uefeito,probs=c(0.01,0.99))
  filter(dfi,Uefeito%in%v_range) %>% 
    select(-logOR,-c(forest_succession:data_year),-SiteCode) %>% 
    mutate(ext_class=ifelse(Uefeito==v_range[1],"min","max"),
           Uefeito_quant=ifelse(ext_class=="min",v_quant[1],v_quant[2])) %>% 
    relocate(ext_class,.after=contraste) %>% 
    relocate(Uefeito_quant,.after=Uefeito)
}) %>% select(-starts_with("p"))
df_buf %>% 
  pivot_longer(starts_with("Uefeito")) %>% 
  filter(name=="Uefeito") %>% 
  ggplot(aes(x=k,y=value,color=ext_class)) +
  geom_point() +
  geom_line() +
  geom_smooth(method="gam",se=FALSE) +
  scale_color_manual(values=c("darkred","darkgreen")) +
  facet_wrap(~contraste,ncol=3)
#############################################################
############## ajuste dos gam para os extremos ##############
#############################################################
f_gam <- \(dfi){
  gam(Uefeito~s(k,bs="cr"),data=dfi)
}
l_md <- dlply(df_buf,c("contraste","ext_class"),f_gam)
for(i in names(l_md)){
  p <- gratia::draw(l_md[[i]],residuals = TRUE) + labs(title=i)
  print(p)
  if(!isTRUE(askYesNo("Do you want to see the residuals?"))) break
  p <- gratia::appraise(l_md[[i]])
  print(p)
  if(!isTRUE(askYesNo("next model?"))) break
}
saveRDS(l_md,file=paste0(v_path,"rds/l_md_refU.rds"))
#############################################################
############## remoção da área do buffer ##############
#############################################################
# i) carregar os dados de l_df_pred
l_paths <- paste0(v_path,"rds/l_dfpred_",c("fragtotal","fragperse","areaperse"),".rds")
l_df_pred <- lapply(l_paths,readRDS) %>% 
  lapply(.,"[[","apenas fixo") %>% 
  lapply(.,rename,k=k_cont)
names(l_df_pred) <- c("Frag. total","Frag. per se","Área per se")
# ii) filtrar os valores únicos de k em um novo data frame
l_df_ref <- lapply(l_df_pred,select,k,SiteCode) %>% 
  lapply(.,distinct)
# iii) fazer a predição média do modelo para cada ponto e guardar em um dataframe
l_df_ref <- lapply(names(l_df_ref),\(li){
  lmd <- l_md[grep(li,names(l_md))]
  names(lmd) <- gsub(paste0(li,"."),"",names(lmd)) 
  df_ref <- lapply(names(lmd),\(i){
    md <- lmd[[i]]
    dfr <- l_df_ref[[li]]  
    dfr[[i]] <- predict.gam(md,dfr)
    return(dfr)
  }) %>% Reduce("inner_join",.)
})
names(l_df_ref) <- names(l_df_pred)
# iv) filtrar para cada k os valores entre as predições
library(data.table)
l_df_new <- lapply(names(l_df_ref),\(li){
  df_pred <- l_df_pred[[li]]
  df_ref <- l_df_ref[[li]]
  setDT(df_pred)
  alply(df_ref,1,\(dfi){
    maxv <- dfi[["max"]]
    minv <- dfi[["min"]]
    kv <- dfi[["k"]]
    df_pred[k==kv & between(Uefeito,minv,maxv),]
  }) %>% rbindlist
})
names(l_df_new) <- names(l_df_ref)
# salvamento
l_paths <- gsub("dfpred","dfnew",l_paths)
names(l_paths) <- names(l_df_ref)
lapply(names(l_df_ref),\(i){
  dfrds <- l_df_new[[i]]
  vpath <- l_paths[i]
  saveRDS(dfrds,file=vpath)
})