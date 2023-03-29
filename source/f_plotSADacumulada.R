# Função para plotar as SADs acumuladas observada e preditas por SiteCode e cenário de limitação de dispersão
# INPUT: um data frame com as colunas:
## SiteCode: coluna guia
## p: proporção de cobertura vegetal
## k: cenário de limitação de dispersão 
## nCongKS, nCongDTS, Smed, Ssd, Smin, Smax: variáveis respostas de interesse - métricas que avaliam as SADs replicas
## SADrep.path: caminho para o .csv com o output da simulação de MNEE
## SADobs.path: caminho ara o .csv com a SAD observada
## effort_ha: área da parcela no qual a SAD observada foi amostrada
## Ntotal: número total de indivíduos vivos na SAD observada
## S_obs: número de espécies na SAD observada
#
# um exemplo de output pode ser obtido em apendices/SADs_neutras.Rmd, df_plot, resultado de um inner_join
#
# OUTPUT: .png com as SADs acumuladas
#
library(ggpmisc)
f_plot.SADacumulada <- function(df,parallel=TRUE,repo="figuras/SADs_acumuladas/"){
  # dados
  m_SADrep <- read_csv(df$SADrep.path[1])
  m_SADrep$k <- factor(m_SADrep$k,levels=unique(m_SADrep$k))
  df_SADobs <- read_csv(df$SADobs.path[1])
  # subfunção
  f_SAD <- function(X){
    v_SADrep <- unname(sort(table(as.matrix(X[,-(1:2)])),decreasing = T))
    cbind(X[,1:2],v_SADrep)
  }
  # rotina
  registerDoMC(2)
  df_return <- rename(adply(m_SADrep,1,f_SAD,.expand = F,.parallel = parallel),n=X1,sp.name=Var1)
  # labels
  v_labels <- paste0("k=",round(df$k,2))
  names(v_labels) <- df$k
  # tabelas inseridas
  df_text.table <- df |> 
    select(k:Ssd) |> 
    mutate(Smed = round(Smed,2),
           Ssd = round(Ssd,2)) |> 
    pivot_longer(-k,names_to="variáveis") |> 
    group_by(k) |> 
    nest() |> 
    mutate(x=0.95,y=0.50)
  # title
  v_title = paste0(df$SiteCode[1],
                   ", p=",round(df$p[1],2),
                   ", plot=",df$effort_ha[1],"ha",
                   ", Ntotal=",df$Ntotal[1],
                   ", Sobs=",df$S_obs[1])
  # gráfico
  p <- ggplot(df_return,aes(x=Freq)) +
    stat_ecdf(alpha=0.4,aes(group=n)) +
    geom_table_npc(data=df_text.table,aes(npcx = x, npcy = y, label = data)) +
    facet_wrap(~k,ncol=5,labeller = as_labeller(v_labels)) +
    stat_ecdf(data=df_SADobs,aes(x=N),color="red") +
    coord_cartesian(xlim = range(df_SADobs$N),expand = FALSE) +
    scale_x_continuous(breaks=seq(min(df_SADobs$N),max(df_SADobs$N),by=10)) +
    labs(y="",x="indivíduos",
         title=v_title) +
    theme_minimal() 
  png(paste0(repo,df$SiteCode[1],".png"),width = 13.8, height = 7.33,units="in",res=500)
  print(p)
  dev.off()
}