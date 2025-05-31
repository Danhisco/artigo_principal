library(gt);library(webshot2); library(patchwork)
library(magick)
library(grid)
library(png)
library(grid)
library(gtable)
library(ggpubr)
library(gt)
library(flextable)
library(dagitty)
library(ggdag)
library(gratia)
# library(sads)
library(doMC)
# library(metR)
library(kableExtra)
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
# library(mgcv)
library(lme4)
library(data.table)
library(plyr)
library(dplyr)
library(sf)
probs = c(0.05,0.25, 0.5, 0.75,0.95)
l_dfpred <- readRDS(file = "dados/csv_SoE/rds/l_dfpred_md_cong_absoluta.rds")
vsites105 <- l_dfpred$fixo_e_aleat$SiteCode %>% unique
################################################################################################
################################# criacao GE dados disponiveis #################################
################################################################################################
f_z <- \(x) (x-mean(x))/sd(x)
theme_set(theme_light())
# df_dados_disponiveis <- read_csv(file = "../../dados/df_dados_disponiveis.csv")
df_dados_disponiveis <- read_csv(file = "dados/df_dados_disponiveis.csv")
df_p_extensoes <- read_csv("dados/csv/df_p_extensoes.csv")
df_p <- ddply(df_p_extensoes,"SiteCode",\(dfi){
  vfilter <- dfi$lado_km[which.min(abs(4-dfi$lado_km))]
  filter(dfi,lado_km==vfilter)
}) %>% select(-c(Ntotal:S_obs))
df_plot <- df_dados_disponiveis %>% 
  inner_join(df_p) %>% 
  select(SiteCode, p, effort_ha, Ntotal, S_obs, year_bestProxy, forest_succession, lat:long_correct)
# df_p_extensoes <- read_csv("../../dados/csv/df_p_extensoes.csv")
probs = c(0.05,0.25, 0.5, 0.75,0.95)
df_contrastes <- read_csv(file="dados/csv/taxaU/df_contrastes.csv") 
df_md_Uefeito <- df_contrastes %>% select(SiteCode:p, contains("_logratio")) %>% 
  pivot_longer(-c(SiteCode:p)) %>% 
  mutate(name = gsub("_logratio","",name),
         across(c(p,k),f_z,.names = "{.col}_z"),
         SiteCode = factor(SiteCode)) %>% 
  rename("tp_efeito"="name","v_efeito"=value)
shp_FA <- read_sf("dados/shapefile/biome_border.shp")
##############
l_p <- list()
sf_plot <- st_as_sf(df_plot %>% 
                      mutate(lat = ifelse(is.na(lat),lat_correct,lat),
                             long = ifelse(is.na(long),long_correct,long)), 
                    coords = c("long", "lat"), crs = 4326)
sf_plot <- st_transform(sf_plot,crs=st_crs(shp_FA))
l_p[[1]] <- ggplot() +
  geom_sf(data=shp_FA,fill="lightgreen",color="green") +
  geom_sf(data=sf_plot,aes(color=p),size=3) +
  coord_sf(ylim = c(-32,-5),xlim = c(-56,-34)) +
  scale_color_gradientn("% CF\n4x4 km",colours = c("red","black")) +
  guides(color=guide_colourbar(position="inside")) +
  theme(legend.position.inside = c(0.9,0.35)) +
  labs(subtitle = "a)")
v_labels <- c("effort_ha" = "plot area (ha)",
              "Ntotal" = "# ind",
              "S_obs" = "riqueza",
              "year_bestProxy" = "approx. data year")
df_plot[grep("\\_",df_plot$year_bestProxy),"year_bestProxy"] <- "2012.5"
l_p[[2]] <- df_plot |>
  mutate(year_bestProxy = as.numeric(year_bestProxy)) |> 
  pivot_longer(effort_ha:year_bestProxy) |> 
  mutate(name=case_when(
    name == "effort_ha" ~ "plot area (ha)",
    name == "Ntotal" ~ "# ind",
    name == "S_obs" ~ "S obs",
    name == "year_bestProxy" ~ "year"
  )) %>% 
  ggplot(aes(x=name,y=value)) +
  theme_bw() +
  geom_jitter(alpha=0.4) +
  geom_boxplot() +
  labs(x="",y="",subtitle="b)") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
    # axis.text.x = element_blank(),
    #     axis.ticks.x = element_blank()
    ) +
  facet_wrap(~name,nrow=1,scales="free")
# coord_flip()
l_p[[3]] <- df_p_extensoes |> 
  filter(lado_km<=16.02) |> 
  ggplot(aes(x=lado_km,y=p)) + 
  geom_point(size=0.2,alpha=0.2) +
  geom_line(aes(group=SiteCode),linewidth=0.2,alpha=0.2) +
  geom_vline(xintercept = 2^(-1:4),color="darkred",alpha=0.3) +
  labs(x="lado da paisagem quadrada (km)",y="proporção de cobertura florestal (%CF)",
       subtitle="c)") +
  scale_y_continuous(breaks = seq(0,1,by=0.1)) +
  scale_x_continuous(breaks = 2^(-1:4)) +
  coord_cartesian(expand = F) +
  theme_light()
l_p[[4]] <- df_plot |>
  mutate(pert_class = factor(forest_succession,
                             levels=c("primary",
                                      "primary/secondary",
                                      "secondary",
                                      "capoeira"),
                             labels=c("baixa",
                                      "mediana",
                                      "alta",
                                      "altíssima"))) %>% 
  ggplot(aes(x=pert_class)) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.2, colour = "black") +
  labs(x="classificação de perturbação (tempo e grau)",
       y="contagem",subtitle="d)")
# grid.arrange(grobs=l_p,nrow=2)
p <- arrangeGrob(grobs=l_p,nrow=2)
ggsave(
  filename="figuras/GE_dados_disponiveis.png",
  plot = p,
  width=12,height = 10)
################################
#
################################################################################################
##################################### criação de fig5_SoE ######################################
################################################################################################
library(ggnewscale)
library(scales)
df_plot <- read_csv("dados/csv/df_plotSoE.csv")
v_threshold <- paste0("SoE >= ",
                      str_extract(paste(names(df_plot),collapse = " "),
                                  "(?<=SoE\\_)[0-9]{1}[.][0-9]{1,2}")) 
cols <- c("blue","red") # "#f04546","3 quant"=
names(cols) <- c("avg+-sd",v_threshold )
a <- 0.5
df_SoE_k <- df_plot %>% select(k,SoE_0.95) %>% distinct()
p <- df_plot |> 
  ggplot(aes(x=lado_numeric,y=Pre_Umed)) + #,group=k_factor,color=k
  geom_point(alpha=0.40,aes(color=S.N)) +
  geom_line(aes(group=SiteCode,color=S.N),alpha=0.25) +
  scale_color_distiller(palette = "Spectral",name="S / N") +
  new_scale_color() +
  stat_summary(fun.data = mean_se,
               geom ="line",
               aes(group=k_factor,color="avg+-sd"),
               alpha=0.6,size=a) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               aes(x=lado_numeric,
                   y=Pre_Umed,
                   color="avg+-sd"),
               alpha=0.6,size=a) +
  geom_vline(aes(xintercept=SoE_0.95,
                 color=names(cols)[2]),
             alpha=0.5) +
  scale_color_manual(name="",values=cols) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x)) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  facet_wrap(~k_factor,ncol=5) +
  theme_classic() +
  labs(title = "Média da estimativa média da taxa U predito pelo LMM:",
       subtitle="logit(Umed) ~ log(S) * log(N) * k * scale + (1|Site)",
       x="scale (landscape side, km)",
       y="avg U rate (Umed)",
       color="") +
  theme(legend.position = "top")
ggsave(
  file=paste0(v_path,
              "figuras/fig5_SoE.png"),
  plot=p,
  width = 10.2,
  height = 6.7)
#############################################################################
################################ congruência absoluta #######################
#############################################################################
library(hexbin)
p <- l_dfpred$fixo_e_aleat %>% 
  ggplot(aes(x=k,y=nSAD)) +
  # geom_point(alpha=0.3) +
  geom_line(aes(y=fit,group = SiteCode),alpha=0.3) +
  # geom_boxplot(aes(group=k),alpha=0.6) +
  geom_ribbon(data=l_dfpred$fixo,
              aes(x=k,y=fit,ymin=lower,ymax=upper),
              color="blue",
              fill="lightblue") +
  geom_line(data=l_dfpred$fixo,
            aes(x=k,y=fit),
            color="darkred") +
  geom_hex(bins = 50,alpha=0.5) + 
  scale_fill_gradient("contagem",low = "gray", high = "black", na.value = NA) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(x="Proporção de propágulos na vizinhança imediata (k)",
       y="número de SAD simuladas com boa congruência com o observado (nSAD)",
       title="Descrição da congruência absoluta da SAD simulada nas paisagens hipotéticas",
       subtitle="HGAM: nSAD ~ s(k,by=interaction(paisagem hipotética, classe de perturbação) + s(lat,long) + te(k,Sítio,by=paisagem hipotética)") +
  facet_grid(pert_class ~ land) +
  theme(strip.text = element_text(size=12, #margin=margin(),
                                  face="bold"))
ggsave(filename="figuras/descricao_congruencia_absoluta.png",
       p,height = 7.33,width = 13.8)
library(magick)
img <- image_read("figuras/descricao_congruencia_absoluta.png") %>% 
  image_resize("50%") %>% 
  image_trim()
image_write(img,path = "figuras/descricao_congruencia_absoluta.png")
#
#
###################################################################################
##################################### LogU/U ######################################
###################################################################################
# objetos
df_p <- read_csv("dados/df_p.csv")
df_sim <- read_csv("dados/df_simulacao.csv") |> 
  inner_join(x=df_p,by="SiteCode")
df_md <- read_csv("dados/csv_SoE/df_logOR.csv")
df_logUU <- readRDS(file="dados/csv_SoE/df_contrastes_z_e_padraor.rds") %>% 
  select(-contains("_z")) %>% 
  pivot_longer(cols=contains("logratio"),
               values_to="Uefeito",
               names_to = "contraste") %>% 
  inner_join(.,
             select(df_md,SiteCode,forest_succession) %>% distinct) %>% 
  mutate(contraste = gsub("_logratio","",contraste) %>% 
           gsub("\\.","",.) %>% 
           gsub("area","Área per se",.) %>% 
           gsub("fragtotal","Frag. total",.) %>% 
           gsub("fragperse","Frag. per se",.),
         contraste = factor(contraste,
                            levels=c("Frag. total",
                                     "Frag. per se",
                                     "Área per se")),
         pert_class = factor(forest_succession,
                             levels=c("primary",
                                      "primary/secondary",
                                      "secondary"),
                             labels=c("baixa",
                                      "mediana",
                                      "alta"))) %>% 
  select(-forest_succession)
df_SoE_k <- mutate(df_SoE_k,SoE_0.95 = ifelse(is.na(SoE_0.95),1,SoE_0.95))
df_p_extensoes <- read_csv("dados/csv/df_p_extensoes.csv")
df_pk <- ddply(df_SoE_k,"SoE_0.95",\(dfi){
  vref <- unique(dfi$SoE_0.95)
  ddply(df_p_extensoes,"SiteCode",\(dfi){
    vfilter <- dfi$lado_km[which.min(abs(vref-dfi$lado_km))]
    filter(dfi,lado_km==vfilter)
  }) %>% 
    mutate(SoE_0.95 = round(lado_km)) %>% 
    select(-c(Ntotal:S_obs,lado_km)) %>% 
    inner_join(dfi,by="SoE_0.95",relationship="many-to-many")
}) %>% arrange(k)
p <- df_pk %>% 
  inner_join(distinct(select(df_logUU,SiteCode,pert_class))) %>% 
  ggplot(aes(x=k,y=p,color=factor(SoE_0.95),group=SiteCode)) +
  geom_point(alpha=0.2) +
  geom_line() +
  scale_color_manual("Escala Espacial",values=c("red","blue","green")) +
  theme(legend.position = "top") +
  facet_wrap(~pert_class,nrow=1)
ggsave(filename="./1_to_compile_dissertacao_EM_USO/09_SI/figuras_SI/p_k_pert_class.png",
       plot=p,
       width=12,
       height=5,
       dpi=200)
df_p <- ddply(df_p_extensoes,"SiteCode",\(dfi){
  vfilter <- dfi$lado_km[which.min(abs(4-dfi$lado_km))]
  filter(dfi,lado_km==vfilter)
}) %>% select(-c(Ntotal:S_obs))
df_logUU_pk <- inner_join(select(df_logUU,-p),
                          df_pk,by=c("SiteCode","k")) #%>% 
  #inner_join(distinct(select(df_md,SiteCode,k,k_z)))
saveRDS(df_logUU_pk,file="dados/csv_SoE/rds/df_logUUpk.rds")

#########################################################################################
########################### classificação das paisagens hipotéticas quanto ao número de boas congruências
df_ad <- read_csv(file="dados/csv_SoE/df_congruencia_simulacao.csv") %>% 
  filter(SiteCode %in% vsites105) %>% 
  rename(nSAD = nCongKS,
         land = land_type) %>% 
  select(-c(Smed:Smax)) %>% 
  mutate(k=round(k,2),
         SiteCode = factor(SiteCode),
         land = factor(land,
                       levels=c(
                         "cont",
                         "non_frag",
                         "ideal"),
                       labels=c(
                         "fragmentada",
                         "aglomerada",
                         "prístina")),
         congruencia = case_when(
           nSAD >= 95 ~ "muito alta",
           nSAD %in% 76:94 ~ "alta",
           nSAD %in% 25:75 ~ "mediana",
           nSAD %in% 6:24 ~ "baixa",
           nSAD <= 5 ~ "muito baixa"),
         congruencia = factor(congruencia,
                              levels=c("muito alta",
                                       "alta",
                                       "mediana",
                                       "baixa",
                                       "muito baixa")[5:1],
                              labels=c("muito alta\n>=95",
                                       "alta\n75 ~ 95",
                                       "mediana\n25 ~ 75",
                                       "baixa\n5 ~ 25",
                                       "muito baixa\n<=5")[5:1]))
df_plot <- df_ad %>% 
  group_by(land,congruencia) %>% 
  summarise(n=n(),
            nSites=length(unique(SiteCode))) %>% 
  group_by(land) %>% 
  mutate(perc = round(n * 100/sum(n),2)) %>% 
  ungroup() %>% 
  mutate(percSites = round(nSites * 100/105,2)) %>% 
  pivot_longer(starts_with("perc")) %>% 
  mutate(name=ifelse(name=="perc","% simulações","% sítios (n=105)"))
f_ggplot <- \(dfp,lsize=4,tsize=8,legsize=0.25){
  dfp %>% 
    ggplot(aes(x=land,y=congruencia,fill=value)) +
    geom_tile(color="white") +
    geom_label(aes(label=paste0(value,"%")),fill="white",size=lsize) +
    labs(x="",y="",title="Classes de congruência",subtitle="pelo menos 1 simulação, k≥0.50") +
    scale_x_discrete(expand=c(0,0),position = "bottom") +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_viridis_c("%") +
    facet_wrap(~name,scales="free") +
    theme_classic() +
    theme(text=element_text(size=tsize,face="bold"),
          legend.position = "none",
          legend.key.size = unit(legsize, 'cm'),
          plot.margin = margin(0,0,0,0),
          axis.text = element_text(color="black"),
          legend.text = element_text(angle=90))
}
lp <- dlply(df_plot,"name",f_ggplot)
p <- arrangeGrob(grobs=lp,ncol=2)
ggsave(filename="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/heatmap_nSAD_land.png",
       plot=lp[[2]],
       width = 3.5,
       height = 2.15,
       dpi=200)
saveRDS(p,file="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/heatmap_nSAD_land.rds")
img <- image_read("1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/heatmap_nSAD_land.png") %>% 
  image_trim()
image_write(img,"1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/heatmap_nSAD_land.png")
df_ij <- df_ad %>% 
  filter(k>=0.49999) %>%
  group_by(SiteCode) %>% 
  summarise(nSAD_maiorigual75 = all(nSAD>=75))
v_ij <- c("alta congruência sempre"=sum(df_ij$nSAD_maiorigual75),
          "baixa congruência em algum k"=nrow(df_ij)-sum(df_ij$nSAD_maiorigual75))
library(cowplot)
df_logUU_pk <- readRDS(file="dados/csv_SoE/rds/df_logUUpk.rds")
df_logUU_plot <- inner_join(df_logUU_pk,df_ij) %>% 
  mutate(label = ifelse(nSAD_maiorigual75,
                        paste0("alta cong. sempre (≥75%)",", n sítios=",v_ij[1]),
                        paste0("cong. baixa em alguma sim.",", n sítios=",v_ij[2])))
f_ggplot <- \(dfp,tsize1=15,tsize2=10,fsize=3){
  vcolor = c("red","blue")
  names(vcolor) = unique(dfp$label)
  p2 <- dfp %>% 
    filter(k==0.99,contraste=="Frag. total") %>% 
    ggplot(aes(fill=label,x=p)) +
    geom_histogram(position = "dodge",bins=60) +
    labs(x="",y="",title="%CF: sítios com boa cong. sempre") +
    scale_fill_manual("",values=vcolor) +
    guides(fill=guide_legend(override.aes = list(size = fsize),
                             nrow = 2,byrow=TRUE)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_minimal_grid() +
    theme(legend.position = c(0.01,0.9),
          legend.box = "horizontal",
          legend.key.size = unit(0.1, 'cm'),
          legend.key.height = unit(0.1, 'cm'),
          legend.key.width = unit(0.1, 'cm'),
          legend.title = element_text(size=7),
          legend.text = element_text(size=7),
          text = element_text(size=tsize2),
          axis.text = element_text(size=9),
          legend.margin = margin(),
          plot.margin = margin())
  return(p2)
}
p_boasempre_ounao <- f_ggplot(df_logUU_plot,tsize1 = 10,tsize2 = 7.5,fsize = 0.5)
p <- free(lp[[2]],side = "tb") + p_boasempre_ounao
saveRDS(p,file="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/heatmap_nSAD_land.rds")
ggsave(filename="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/sitios_filtrados.png",
       plot=p,
       width = 6,
       height = 2.15,
       dpi=200)
img <- image_read("1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/sitios_filtrados.png") %>% 
  image_trim()
image_write(img,"1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/sitios_filtrados.png")



f_ggplot <- \(dfp,tsize1=15,tsize2=10){
  p1 <- dfp %>% 
    mutate(label=gsub(", n sítios=","",label) %>% 
             gsub("[1-9]{2}","",.) %>% 
             gsub("k","sim.",.) %>% 
             gsub(" \\(\\≥\\%\\)","",.)) %>% 
    filter(k>=0.49999) %>% 
    ggplot(aes(x=k,y=Uefeito,color=p)) +
    geom_boxplot(aes(group=k)) +
    geom_point(alpha=0.6) +
    labs(y="logU/U",x="proporção de propágulos até a vizinhança imediata (k)") +
    geom_line(aes(group=SiteCode),alpha=0.6) +
    scale_colour_gradient2("% CF",midpoint=0.5,
                           low="red",
                           mid = "yellow",
                           high = "darkgreen") +
    facet_grid(label~contraste) +
    theme_classic() +
    theme(text=element_text(size=tsize1),
          legend.position = "top",
          legend.title = element_text(size=tsize2),
          legend.text = element_text(size=tsize2))
  return(p1)
}
p <- f_ggplot(df_logUU_plot,tsize1 = 15)
ggsave(filename="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/pedacofigfinal_1alinha.png",
       plot=p,
       width = 12,
       height = 7,
       dpi=200)
img <- image_read("1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/pedacofigfinal_1alinha.png") %>% 
  image_trim()
image_write(img,"1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/pedacofigfinal_1alinha.png")
saveRDS(p,file="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/pedacofigfinal_1alinha.rds")
df_md <- select(df_logUU_plot[df_logUU_plot$nSAD_maiorigual75,],
                SiteCode:Uefeito,p)
saveRDS(df_md,"dados/csv_SoE/rds/df_logUUpk.rds")
#########################################################################################
########################### logU/U : área per se e fragmentçaão per se ##################
#########################################################################################
f_gsub <- \(xlab){
  gsub("area","Área per se",xlab) %>% 
    gsub("fragperse","Frag. per se",.) %>% 
    gsub("fragtotal","Frag. total",.)
}

df_plot <- df_logUU_plot %>% filter(contraste!="Frag. total",nSAD_maiorigual75)
label_efeitos="fragmentação per se ~ área per se"
xlab="Área per se"
ylab="Frag. per se"
l_fplot <- list()
l_fplot$efeito_x_efeito <- \(df_plot,
                             label_efeitos=NULL,
                             xlab = NULL,
                             ylab = NULL){
  #@ xylab = "area","fragperse","fragtotal"
  if(is.null(label_efeitos)) stop("precisa fornecer label_efeitos = e.g. fragmentação per se ~ área per se")
  if(is.null(xlab)|is.null(ylab)) stop("precisa fornecer label_efeitos = e.g. fragmentação per se ~ área per se")
  geom_list1 <- list(
    geom_hline(yintercept = 0,color="black"),
    geom_vline(xintercept = 0,color="black"),
    geom_abline(slope=1,intercept=0,color="darkblue",linetype=1),
    geom_abline(slope=-1,intercept=0,color="darkblue",linetype=1),
    geom_hex(bins = 15,alpha=0.5),
    scale_fill_gradient("contagem",low = "yellow", high = "red", na.value = NA),
    # facet_wrap(~pert_class,ncol=3),
    theme_classic(),
    theme(aspect.ratio = 1)
  )
  p <- df_plot %>% 
    pivot_wider(names_from=contraste,values_from=Uefeito) %>%   
    ggplot(aes(x=.data[[xlab]],y=.data[[ylab]])) +
    geom_list1 +
    labs(title=paste0("logU/U :",label_efeitos),
         x=f_gsub(xlab),
         y=f_gsub(ylab))
  # saveRDS(p,file="figuras/logUU_construcao/p_fragperse_areaperse.rds")
  return(p)
}
l_fplot$tabela_contagem_maior_magnitude_absoluta <- \(dfi){
  dfi %>% 
    group_by(class_diff) %>% 
    tally() %>% 
    mutate(n_rel = round(n * 100 / sum(n),2)) %>% 
    ggplot(aes(x=class_diff,y="",fill=n_rel)) +
    geom_raster() +
    labs(x="",y="",
         title="") +
    scale_fill_gradient("",low = "yellow", high = "red", na.value = NA) +
    geom_text(aes(label=n_rel),size=10) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand =  c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=15),
          axis.ticks.y = element_blank(),
          title = element_text(size=15),
          strip.text = element_text(size=15),
          legend.position = "none")
}

dfp <- df_plot %>% 
  mutate(Uefeito_abs = abs(Uefeito)) %>% 
  select(-Uefeito) %>% 
  ddply(.,c("k","SiteCode"),\(dfi){
    dfi$diff_efeito <- with(dfi,{
      Uefeito_abs[contraste=="Área per se"] -
        Uefeito_abs[contraste=="Frag. per se"]
    })
    return(dfi)
  }) %>% pivot_wider(names_from=contraste,values_from=Uefeito_abs) %>%
  mutate(class_diff = ifelse(diff_efeito>0,
                             paste("Área per se", ">", "Frag. per se"),
                             paste("Frag. per se", ">", "Área per se")))
dfi <- dfp


l_fplot$histograma_efeito_menos_efeito <- \(dfi,xlabel=NULL){
  if(is.null(xlabel)) stop("precisa fornecer xlabel = e.g. logU/U: |área per se| - |frag. per se|")
  geom_list2 <- list(
    geom_histogram() ,
    geom_vline(xintercept = 0,color="red",linetype=2) ,
    theme_classic() ,
    scale_x_continuous(expand=c(0,0)) ,
    scale_y_continuous(expand=c(0,0))
  )
p <- dfi %>% 
    ggplot(aes(x=diff_efeito)) +
    geom_list2 +
    labs(x=xlabel,
         y="") 
saveRDS(p,file="figuras/logUU_construcao/p_hist_diff_fragperse_areaperse.rds")
}
#
f_plots <- \(comp_efeitos=c("fragperse ~ area",
                            "fragtotal ~ fragperse",
                            "fragtotal ~ area")){
  #@ y ~ x
  l_s <- list()
  l_s$filtro <- gsub("fragperse ~ area","Frag. total",comp_efeitos) %>% 
    gsub("fragtotal ~ fragperse","Área per se",.) %>% 
    gsub("fragtotal ~ area","Frag. per se",.)
  l_s$labels <- str_split_1(comp_efeitos," \\~ ")
  names(l_s$labels) <- c("y","x")
  l_s$levels <- gsub("area","Área per se",l_s$labels) %>% 
    gsub("fragperse","Frag. per se",.) %>% 
    gsub("fragtotal","Frag. total",.)
  names(l_s$levels) <- c("y","x")
  l_s$axis_values <- tolower(l_s$levels) %>% 
    gsub("área per se","área",.) %>% 
    paste0("|",.,"|")
  names(l_s$axis_values) <- c("y","x")
  # 
  l_df <- list()
  l_df[[1]] <- df_logUU_pk %>% 
    filter(contraste!=l_s$filtro) %>% 
    mutate(contraste = factor(contraste,
                              levels=l_s$levels,
                              labels=l_s$labels))
  l_df[[2]] <- l_df[[1]] %>% 
    mutate(Uefeito_abs = abs(Uefeito)) %>% 
    select(-Uefeito) %>% 
    ddply(.,c("k","SiteCode"),\(dfi){
      dfi$diff_efeito <- with(dfi,{
        Uefeito_abs[contraste==l_s$labels["x"]] -
          Uefeito_abs[contraste==l_s$labels["y"]]
        })
      return(dfi)
    }) %>% pivot_wider(names_from=contraste,values_from=Uefeito_abs) %>%
    mutate(class_diff = ifelse(diff_efeito>0,
                               paste(l_s$axis_values["x"], ">", l_s$axis_values["y"]),
                               paste(l_s$axis_values["y"], ">", l_s$axis_values["x"]))) 
  #
  l_p <- list()
  l_p[[1]] <- l_fplot[[1]](l_df[[1]],
                           label_efeitos = paste0(l_s$levels["y"]," ~ ",l_s$levels["x"]),
                           xlab = l_s$labels["x"],
                           ylab = l_s$labels["y"])
  l_p[[2]] <- l_fplot[[2]](l_df[[2]])
  l_p[[3]] <- l_fplot[[3]](l_df[[2]],
                           xlabel = paste0("logU/U: ",l_s$levels["x"]," - ",l_s$levels["y"]))
  #
  names(l_p) <- names(l_fplot)
  return(l_p)
}
v_comparacao_efeitos <- c("fragperse ~ area",
            "fragtotal ~ fragperse",
            "fragtotal ~ area")
l_p <- lapply(v_comparacao_efeitos,f_plots)
names(l_p) <- v_comparacao_efeitos
saveRDS(l_p,file="./1_to_compile_dissertacao_EM_USO/09_SI/figuras_SI/l_p_comparacao_empirica_efeitos_diferenca_entre_pares.rds")
#####################################################################
################# objetos criados em 'ajuste_GAMM.R"
###################################################################
l_p_apenas_fixo <- readRDS(file="./figuras/logUU_construcao/l_p_apenas_fixo.rds")
l_p_fixo_e_aleat <- readRDS(file="./figuras/logUU_construcao/l_p_fixo_e_aleat.rds")
###
# apenas fixo
###
p <- arrangeGrob(grobs=l_p_apenas_fixo,nrow=1)
ggsave(filename="figuras/todos_efeitos_medios_desconsiderandoSitios.png",plot=p,
       width = 15,
       height = 7)
#####################################################################################
################# objetos criados em 'ajuste_GAMM.R"
# dados dos observados e preditos para os sítios
df_real <- readRDS("dados/csv_SoE/rds/df_obs_pred_plot_logUU_pk.rds") %>% 
  do.call("rbind",.) %>% 
  mutate(
    p_class = case_when(
      p==1 ~ "%CF = 100",
      p<1 & p>=0.80 ~ "80 ≤ %CF < 100",
      p<0.80 & p>=0.60 ~ "60 ≤ %CF < 80",
      p<0.60 & p>=0.30 ~ "30 ≤ %CF < 60",
      p<0.30 ~ "%CF < 30"),
    p_class = factor(p_class,levels=c("%CF < 30",
                                      "30 ≤ %CF < 60",
                                      "60 ≤ %CF < 80",
                                      "80 ≤ %CF < 100",
                                      "%CF = 100")),
    lower = fit - 1.96 * se.fit,
    upper = fit + 1.96 * se.fit
  )
# paisagens sem perda de cobertura florestal
df_refCF1 <- df_real %>% 
  filter(p==1) %>% 
  reframe(
    names = c("min","max","Q5%","Q95%"),
    value = c(min(Uefeito), max(Uefeito),quantile(Uefeito,probs=0.05),quantile(Uefeito,probs=0.95))
  )
# comparação das paisagens sem perda com outras paisagens: elas são maiores do que o quantil de 90% dos %CF=1? 
df_sites <- ddply(df_refCF1,"names",\(dfi){
  #
  df_ref <- df_real %>% filter(p!=1)
  vsig <- ifelse(dfi$names%in%c("min","Q5%"),"<",">")
  #
  dfret <- ddply(df_ref,c("efeito","p_class"),\(dfd){
    vlog <- eval(
      parse(text=paste("dfd$Uefeito",vsig,dfi$value[1]))
    )
    dfd[vlog,c("efeito","p_class","SiteCode")] %>% 
      distinct() %>% 
      mutate(cond_class = paste(vsig,dfi$names))
  })
}) %>% select(-names) %>% relocate(cond_class)
# qual a % de sítios que ultrapassa essa região de ausência de diferença?
df_sites_summary <- df_sites %>% 
  group_by(efeito,p_class,cond_class) %>% 
  tally() %>% 
  inner_join(df_sites %>% 
               group_by(p_class) %>% 
               summarise(nSite = unique(SiteCode) %>% length())) %>% 
  ungroup() %>% 
  mutate(perc = round(n*100/nSite,1),
         cond_class = factor(cond_class,
                             levels=c(
                               "< min","< Q5%",
                               "> Q95%","> max"
                             )),
         efeito = factor(efeito,
                         levels=c(
                           "Frag. total",
                           "Frag. per se",
                           "Área per se"),
                         labels=c(
                           "FT",
                           "FPS",
                           "APS"
                         )))
###################
df_ij <- df_real %>% 
  select(SiteCode,p_class) %>% 
  distinct() %>% 
  group_by(p_class) %>% 
  tally() %>% 
  mutate(label = paste0(p_class,", n Sítios=",n)) %>% 
  select(-n)
df_ij$label <- factor(df_ij$label,
                      levels=unique(df_ij$label))
dfp <- inner_join(df_sites_summary,df_ij)



# Gráficos
p <- dfp %>% 
  ggplot(aes(x=cond_class,y=factor(1))) +
  geom_col(aes(,fill=perc)) +
  scale_fill_gradient("% sítios",low="gray",high="red") +
  # ggnewscale::new_scale_fill() +
  geom_label(aes(label=paste0(perc,"%"),y=0.5)) +
  labs(x="em comparação com logU/U quando %CF=100",y="",
       title="Sítios que ultrapassam os efeitos da paisagem sem perda de cobertura florestal em pelo menos uma sim.",
       caption="FT = frag. total; FPS = frag. per se; APS = área per se") +
  scale_y_discrete(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  facet_grid(efeito~label,scales="free") +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(face="bold"),
        plot.caption = element_text(lineheight = 0.5,size=7),
        plot.title = element_text(hjust=0.5,face="bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(filename = "1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/propSitios_ultrapassam_efeitoCF100.png",
       plot=p,
       width = 12,
       height = 2.23,
       dpi=200)
saveRDS(p,file = "1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/propSitios_ultrapassam_efeitoCF100.rds")
######### obs X predito (fixo + aleatório)
f_obs_predito_bysite <- \(dff,logic_ribbon=FALSE){
  dfpred <- dff %>% 
    mutate(
      p_class = case_when(
        p==1 ~ "%CF = 100",
        p<1 & p>=0.80 ~ "80 ≤ %CF < 100",
        p<0.80 & p>=0.60 ~ "60 ≤ %CF < 80",
        p<0.60 & p>=0.30 ~ "30 ≤ %CF < 60",
        p<0.30 ~ "%CF < 30"),
      p_class = factor(p_class,levels=c("%CF < 30",
                                        "30 ≤ %CF < 60",
                                        "60 ≤ %CF < 80",
                                        "80 ≤ %CF < 100",
                                        "%CF = 100")),
      lower = fit - 1.96 * se.fit,
      upper = fit + 1.96 * se.fit
    )
  df_ij <- dfpred %>% 
    select(SiteCode,p_class) %>% 
    distinct() %>% 
    group_by(p_class) %>% 
    tally() %>% 
    mutate(label = paste0(p_class,", n Sítios=",n)) %>% 
    select(-n)
  df_ij$label <- factor(df_ij$label,
                        levels=unique(df_ij$label))
  dfp <- inner_join(dfpred,df_ij)
  nsites <- filter(dfp,p_class=="%CF = 100") %>% pull(SiteCode) %>% unique %>% length()
  dfrange <- filter(dfp,p_class=="%CF = 100") %>% 
    reframe(
      names = c("amplitude","amplitude","Q 5%-95%","Q 5%-95%"),
      value = c(min(Uefeito), max(Uefeito),quantile(Uefeito,probs=0.05),quantile(Uefeito,probs=0.95))
    )
  # across(Uefeito,
  #            .fns = list(min=min,max=max),
  #            .names = "{.fn}")
  # 
  dfp <- filter(dfp,p_class!="%CF = 100")
  dfp %>% 
    mutate(efeito = factor(efeito,
                           levels=c("Frag. total",
                                    "Frag. per se",
                                    "Área per se"))) %>% 
    ggplot(aes(x=k,y=Uefeito)) +
    geom_boxplot(aes(group=k),alpha=0.6) +
    geom_hline(yintercept = 0,color="black") +
    geom_hline(data = dfrange, aes(yintercept = value,linetype=names),
               color="black",linewidth=0.5) +
    scale_linetype_manual(name = paste0("logU/U %CF=100",", n sítios=",nsites),
                          values=c("dotted","dashed")) +
    {if(logic_ribbon)geom_ribbon(aes(x=k,y=fit,ymin=lower,ymax=upper,group = SiteCode),
                                 fill="lightblue",
                                 alpha=0.2)}+
    geom_line(aes(y=fit,group = SiteCode,color=p),alpha=0.4) + #,color="darkred"
    geom_point(aes(color=p),alpha=0.85) +
    scale_colour_gradient2("% CF",midpoint=0.5,
                           low="red",
                           mid = "yellow",
                           high = "darkgreen") +
    labs(x="prop. de propágulos na vizinhança imediata (k)",y="logU/U",title="") +
    facet_grid(efeito~label,scales="free") +
    theme_classic() +
    theme(plot.margin=unit(c(0,0.2,0,0), "cm"),
          legend.position = "inside",
          legend.position.inside = c(0.72,0.25),
          legend.direction="horizontal") 
}
p_efeitos <- f_obs_predito_bysite(df_real,logic_ribbon = TRUE)
ggsave(filename = "1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/p_efeitos_porclasseCF.png",
       dpi=200,
       plot=p_efeitos,
       width=12,height=8)
saveRDS(p_efeitos,file="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/p_efeitos_porclasseCF.rds")
###############
#
#### Sorteio de sítios para ilustrar os maiores desvios da paisagem
df_sites <- df_real %>% 
  filter(p_class!="%CF = 100") %>% 
  group_by(efeito,p_class,SiteCode) %>% 
  summarise(max_logUU = max(Uefeito))
df_sites_q5_25_75_95 <- ddply(df_sites,c("p_class"),\(dfi){
  vqlogUU <- quantile(dfi$max_logUU,probs=c(0.05,0.25,0.75,0.95))
  vsites <- sapply(vqlogUU,\(x){
    dfi$max_logUU[which.min(abs(x-dfi$max_logUU))]
  })
  filter(dfi,max_logUU%in%vsites) %>% 
    select(-max_logUU)
}) %>% select(-efeito) %>% distinct()
vsites <- df_sites_q5_25_75_95 %>% filter(p_class=="%CF < 30") %>% pull(SiteCode)
vsites <- df_sites %>% 
  filter(p_class=="%CF < 30" & !(SiteCode%in%vsites)) %>% 
  pull(SiteCode) %>% unique %>% sample()
# olhar f_SoE em source/dinamica_coalescente.R
df_sites_q5_25_75_95 <- rbind(df_sites_q5_25_75_95,
                              data.frame(p_class="%CF < 30",
                                         SiteCode=vsites[1])) %>% 
  arrange(p_class)
a_ply(df_sites_q5_25_75_95,1,\(dfi){
  porigem <- paste0("dados/simulacao/",dfi$SiteCode,".txt")
  pdestino <- gsub("simulacao","mapas_sitios_gdrive",porigem)
  file.copy(from = porigem, to = pdestino)
})
saveRDS(df_sites_q5_25_75_95,"dados/csv_SoE/rds/df_sites_q5_25_75_95.rds")
######################################################################################
# 1a versão figura completa:
f_circle_plot <- \(dfp,boxlabelsize=10,tsize=15,lsize=15){
  dfp2 <- mutate(
    dfp,
    region_id = case_when(
      A_maior0 == "Área per se>=0" & F_maior0 == "Frag per se>=0" & abs_AmaiorF == "|Área per se|>|Frag per se|" ~ 1,
      A_maior0 == "Área per se>=0" & F_maior0 == "Frag per se>=0" & abs_AmaiorF == "|Frag per se|>|Área per se|" ~ 2,
      A_maior0 == "Área per se<0" & F_maior0 == "Frag per se>=0" & abs_AmaiorF == "|Frag per se|>|Área per se|" ~ 3,
      A_maior0 == "Área per se<0" & F_maior0 == "Frag per se>=0" & abs_AmaiorF == "|Área per se|>|Frag per se|" ~ 4,
      A_maior0 == "Área per se<0" & F_maior0 == "Frag per se<0" & abs_AmaiorF == "|Área per se|>|Frag per se|" ~ 5,
      A_maior0 == "Área per se<0" & F_maior0 == "Frag per se<0" & abs_AmaiorF == "|Frag per se|>|Área per se|" ~ 6,
      A_maior0 == "Área per se>=0" & F_maior0 == "Frag per se<0" & abs_AmaiorF == "|Frag per se|>|Área per se|" ~ 7,
      A_maior0 == "Área per se>=0" & F_maior0 == "Frag per se<0" & abs_AmaiorF == "|Área per se|>|Frag per se|" ~ 8)) %>% 
    ungroup() %>% 
    # select(-starts_with("perc_")) %>% 
    mutate(perc=round(n*100/sum(n),2))
  #
  angles <- seq(0, 2*pi, length.out = 9)
  #
  circle_data <- data.frame()
  for(i in 1:9){
    segment <- data.frame(
      region_id = i,
      x = c(0, cos(angles[i]), cos(angles[i+1])),
      y = c(0, sin(angles[i]), sin(angles[i+1]))
    )
    circle_data <- rbind(circle_data, segment)
  }
  # 
  plot_data <- circle_data %>%
    left_join(dfp2, by = "region_id")
  #
  segment_data <- data.frame(
    x = 0,
    y = 0,
    xend = cos(angles),
    yend = sin(angles)
  )
  #
  df_label <- plot_data %>% 
    group_by(region_id) %>% 
    summarise(x_lab=mean(x),y_lab=mean(y)) %>% 
    inner_join(select(plot_data,region_id,perc,perc_sites) %>% distinct()) %>% 
    mutate(across(starts_with("perc"),~paste0(.x,"%")))
  f_ggplot <- \(vfill){
    vtitle=ifelse(grepl("site",vfill),"% de sítios","% de simulações")
    plot_data %>% 
      mutate(facet_label = vtitle) %>% 
      ggplot() +
      geom_polygon(
        aes(x = x, y = y, group = region_id, fill = .data[[vfill]]),
        color = "white", linewidth = 0.5
      ) +
      scale_fill_gradient("%",low = "gray", high = "brown") +
      geom_label(data=df_label,
                 aes(label=.data[[vfill]],x=x_lab,y=y_lab),
                 fill="white",
                 size=boxlabelsize) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() +
      facet_wrap(~facet_label) +
      theme(text=element_text(size=tsize),
            plot.margin = margin(0,0,0,0),
            strip.text = element_text(size=lsize,face="bold"),
            legend.position = "none",
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank())  
  }
  lp <- lapply(c("% sim"="perc","% sitios"="perc_sites"),f_ggplot)
  library(patchwork)
  lp[[1]]+lp[[2]]
}

f_circle_plot2 <- \(dfp,df_ref,boxlabelsize=10,tsize=15,lsize=15){
  dfp2 <- mutate(
    dfp,
    region_id = case_when(
      A_maior0 == "Área per se>=0" & F_maior0 == "Frag per se>=0" & abs_AmaiorF == "|Área per se|>|Frag per se|" ~ 1,
      A_maior0 == "Área per se>=0" & F_maior0 == "Frag per se>=0" & abs_AmaiorF == "|Frag per se|>|Área per se|" ~ 2,
      A_maior0 == "Área per se<0" & F_maior0 == "Frag per se>=0" & abs_AmaiorF == "|Frag per se|>|Área per se|" ~ 3,
      A_maior0 == "Área per se<0" & F_maior0 == "Frag per se>=0" & abs_AmaiorF == "|Área per se|>|Frag per se|" ~ 4,
      A_maior0 == "Área per se<0" & F_maior0 == "Frag per se<0" & abs_AmaiorF == "|Área per se|>|Frag per se|" ~ 5,
      A_maior0 == "Área per se<0" & F_maior0 == "Frag per se<0" & abs_AmaiorF == "|Frag per se|>|Área per se|" ~ 6,
      A_maior0 == "Área per se>=0" & F_maior0 == "Frag per se<0" & abs_AmaiorF == "|Frag per se|>|Área per se|" ~ 7,
      A_maior0 == "Área per se>=0" & F_maior0 == "Frag per se<0" & abs_AmaiorF == "|Área per se|>|Frag per se|" ~ 8)) %>% 
    ungroup() %>% 
    # select(-starts_with("perc_")) %>% 
    mutate(perc=round(n*100/sum(n),2))
  #
  angles <- seq(0, 2*pi, length.out = 9)
  #
  circle_data <- data.frame()
  for(i in 1:9){
    segment <- data.frame(
      region_id = i,
      x = c(0, cos(angles[i]), cos(angles[i+1])),
      y = c(0, sin(angles[i]), sin(angles[i+1]))
    )
    circle_data <- rbind(circle_data, segment)
  }
  # 
  plot_data <- circle_data %>%
    left_join(dfp2, by = "region_id")
  #
  segment_data <- data.frame(
    x = 0,
    y = 0,
    xend = cos(angles),
    yend = sin(angles)
  )
  #
  df_label <- plot_data %>% 
    group_by(region_id) %>% 
    summarise(x_lab=mean(x),y_lab=mean(y)) %>% 
    inner_join(select(plot_data,region_id,perc,perc_sites) %>% distinct()) %>% 
    mutate(across(starts_with("perc"),~paste0(.x,"%")))
  f_ggplot <- \(vfill){
    vtitle=ifelse(grepl("site",vfill),"% de sítios","% de simulações")
    plot_data %>% 
      mutate(facet_label = vtitle) %>% 
      ggplot() +
      geom_polygon(
        aes(x = x, y = y, group = region_id, fill = .data[[vfill]]),
        color = "white", linewidth = 0.5
      ) +
      scale_fill_gradient("%",low = "gray", high = "brown") +
      geom_label(data=df_label,
                 aes(label=.data[[vfill]],x=x_lab,y=y_lab),
                 fill="white",
                 size=boxlabelsize) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw() +
      facet_wrap(~facet_label) +
      theme(text=element_text(size=tsize),
            plot.margin = margin(0,0,0,0),
            strip.text = element_text(size=lsize,face="bold"),
            legend.position = "none",
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank())  
  }
  lp <- lapply(c("% sim"="perc","% sitios"="perc_sites"),f_ggplot)
  library(patchwork)
  lp[[1]]+lp[[2]]
}

f_hist_ggplot <- \(dff=df_real,tsize=10,diffrange){
  dff <- dff %>% 
    #select(-fit,-se.fit) %>% 
    filter(efeito!="Frag. total") %>% 
    mutate(Uefeito_abs = abs(Uefeito)) %>% 
    select(-Uefeito) %>% 
    ddply(.,c("SiteCode"),\(dfi){
      dfi$diff_efeito <- with(dfi,{
        Uefeito_abs[efeito=="Área per se"] -
          Uefeito_abs[efeito=="Frag. per se"]
      })
      return(dfi)
    }) %>% 
    pivot_wider(names_from=efeito,values_from=Uefeito_abs) %>%
    mutate(class_diff = ifelse(diff_efeito>0,
                               paste("Área per se", ">", "Frag. per se"),paste("Frag. per se", ">", "Área per se")))
  geom_list2 <- list(
    geom_histogram() ,
    geom_vline(xintercept = 0,color="red",linetype=2) ,
    theme_classic() ,
    scale_x_continuous(expand=c(0,0),limits = diffrange) ,
    scale_y_continuous(expand=c(0,0))
  )
  p <- dff %>% 
    mutate(flab = "|APS| - |FPS|") %>% 
    ggplot(aes(x=diff_efeito)) +
    geom_list2 +
    labs(x="",y="") +
    facet_wrap(~flab) +
    theme(text=element_text(size=tsize,face="bold"),
          plot.margin = margin(0,0,0,0))
  return(p)
}
f_polig_ref <- \(df_ref,df_analise){
  ## 1. Calcular o polígono convexo do grupo de referência
  hull <- df_ref[chull(df_ref[["Área per se"]], df_ref[["Frag. per se"]]), ]
  hull_analise <- df_analise[chull(df_analise[["Área per se"]], df_analise[["Frag. per se"]]), ]
  ## 2. Converter para objeto sf (simple features) para análise espacial
  polygon <- st_as_sf(hull, coords = c("Área per se", "Frag. per se")) %>% 
    summarise(geometry = st_combine(geometry)) %>% 
    st_cast("POLYGON")
  points <- st_as_sf(df_analise, coords = c("Área per se", "Frag. per se"))
  ## 3. Verificar quais pontos estão dentro do polígono
  inside <- st_contains(polygon, points, sparse = FALSE)[1,]
  ## 4. Calcular a porcentagem
  percent_inside <- mean(inside) * 100
  label_text <- sprintf("%.1f%% dentro", percent_inside)
  ## 5. Encontrar o centroide do polígono para posicionar o label
  centroid <- st_centroid(polygon)
  centroid_coords <- st_coordinates(centroid)
  list(
    "ref_hull"=hull,
    "pontos_hull"=hull_analise,
    "centro_coords"=centroid_coords,
    "label"=label_text
  )
}
f_scatterplot <- \(dfp,xrange,yrange,axislabs=TRUE,df_ref,
                   tsize=15){
  label_efeitos="fragmentação per se ~ área per se"
  xlab="Área per se"
  ylab="Frag. per se"
  f_gsub <- \(xlab){
    gsub("area","APS",xlab) %>% 
      gsub("fragperse","FPS",.) %>% 
      gsub("fragtotal","Frag. total",.)
  }
  geom_list1 <- list(
    geom_hline(yintercept = 0,color="black"),
    geom_vline(xintercept = 0,color="black"),
    geom_abline(slope=1,intercept=0,color="darkblue",linetype=1),
    geom_abline(slope=-1,intercept=0,color="darkblue",linetype=1),
    geom_point(alpha=1),
    scale_x_continuous(limits = xrange),
    scale_y_continuous(limits = yrange),
    geom_hex(bins = 60,alpha=0.5),
    scale_fill_gradient("número de\nsimulaçoes",low = "yellow", high = "red", na.value = NA),
    theme_classic(),
    guides(fill=guide_colourbar(position="inside")),
    theme(legend.position = "none",
          legend.position.inside = c(0.7, 0.9),
          legend.direction="horizontal",
          # axis.title = element_blank(),
          aspect_ratio=1,
          plot.margin = margin(0,0,0,0),
          text=element_text(size=tsize,face="bold"))
  )
  # valores de referência
  l_poligons <- f_polig_ref(df_ref=df_ref,df_analise=dfp)
  #
  p <- dfp %>% 
    ggplot(aes(x=.data[[xlab]],y=.data[[ylab]])) +
    geom_polygon(data=l_poligons[["pontos_hull"]],
                 aes(x=.data[[xlab]],y=.data[[ylab]]),
                 fill="#90EE90",color="black")+
    geom_polygon(data=l_poligons[["ref_hull"]],
                 aes(x=.data[[xlab]],y=.data[[ylab]]),
                 fill="darkgray",color=NA,alpha=0.75) +
    labs(x=ifelse(axislabs,f_gsub(xlab),""),
         y=ifelse(axislabs,f_gsub(ylab),"")) +
    geom_list1
  p_reg <- ggExtra::ggMarginal(p, type = "boxplot",
                               fill = "steelblue", col = "darkblue",size=17)
  return(p_reg)
}

library(cowplot)
library(patchwork)



####### versão mais simples 

f_figgeral <- \(df_plot){
  dfp0 <- df_plot %>% 
    select(-`Frag. total`) %>% 
    mutate(k=factor(round(k,2)),
           A_maior0 = ifelse(`Área per se`>=0,"Área per se>=0","Área per se<0"),
           A_maior0 = factor(A_maior0,levels=c("Área per se<0","Área per se>=0")),
           F_maior0 = ifelse(`Frag. per se`>=0,"Frag per se>=0","Frag per se<0"),
           F_maior0 = factor(F_maior0,levels=c("Frag per se>=0","Frag per se<0")),
           abs_AmaiorF = case_when(
             abs(`Área per se`) > abs(`Frag. per se`) ~ "|Área per se|>|Frag per se|",
             abs(`Área per se`) < abs(`Frag. per se`) ~ "|Frag per se|>|Área per se|",
             abs(`Área per se`) == abs(`Frag. per se`) ~ "|frag|=|area|"))
  df_ref <- dfp0 %>% filter(p_class=="%CF = 100")
  df_pclass <-  dfp0 %>% 
    filter(p_class!="%CF = 100",
           k%in%c("0.99","0.75","0.5"))
  xrange <- range(df_pclass$`Área per se`)
  yrange <- range(df_pclass$`Frag. per se`)
  diffrange <- df_pclass %>% 
    select(SiteCode,`Área per se`,`Frag. per se`) %>% 
    pivot_longer(-SiteCode) %>% 
    pull(value) %>% range()
  diffrange[2] <- 0.4
  lp <- dlply(df_pclass,c("p_class","k"),\(dfi){
    # 
    dfp <- dfi %>% 
      select(SiteCode,`Área per se`,`Frag. per se`) %>% 
      pivot_longer(-SiteCode,names_to="efeito",values_to="Uefeito")
      # group_by(A_maior0,F_maior0,abs_AmaiorF) %>% 
      # summarise(n=n(),
      #           nSite=length(unique(SiteCode))) %>% 
      # group_by(A_maior0,F_maior0) %>% 
      # mutate(perc = round(n*100/sum(n),2),
      #        perc_sites = round(nSite*100/67,2))
    # 
    p_reg <- f_scatterplot(dfp = dfi,
                           xrange = xrange, yrange = yrange,
                           df_ref=df_ref,
                           axislabs = TRUE,
                           tsize = 8)
    #p_circ <- f_circle_plot2(dfp,boxlabelsize=3.5,tsize=6,lsize = 8)
    p_hist <- f_hist_ggplot(dff = dfp,
                            tsize = 8,
                            diffrange=diffrange)
    p_final <- wrap_plots(A=p_reg,B=p_hist,widths = c(1.5,1))
    return(p_reg)  
  })
  matnames <- matrix(names(lp),ncol=4,byrow=FALSE)
  pfinal <- wrap_plots(lp,ncol=4,byrow = FALSE,axes = "collect")
  # inclusão dos títulos
  # Definir um layout com:
  # - 4 linhas: 3 para os gráficos + 1 para os títulos das colunas
  # - 4 colunas: 3 para os gráficos + 1 para os títulos das linhas
  layout_matrix <- rbind(
    1:5,   # Títulos das colunas (acima dos gráficos)
    6:10,    # Linha 1 de gráficos
    11:15, # Linha 2 de gráficos
    16:20 # Linha 3 de gráficos
  )
  # Alturas e larguras relativas:
  hstrip <- 0.1
  heights <- c(hstrip, 1, 1, 1)    # 10% para títulos das colunas, resto para gráficos
  widths <- c(1, 1, 1, 1, hstrip)     # 10% para títulos das linhas (direita)
  #
  # Lista de todos os elementos (títulos + gráficos)
  elementos <- list(
    # Títulos das colunas (linha 1 do layout)
    textGrob(vcolunas[1], gp = gpar(fontsize = 12, fontface = "bold")),
    textGrob(vcolunas[2], gp = gpar(fontsize = 12, fontface = "bold")),
    textGrob(vcolunas[3], gp = gpar(fontsize = 12, fontface = "bold")),
    textGrob(vcolunas[4], gp = gpar(fontsize = 12, fontface = "bold")),
    textGrob("", gp = gpar(fontsize = 12)),  # Espaço vazio (canto superior direito)
    # Gráficos + títulos das linhas (à direita)
    ## 1a linha
    lp[[ matnames[1,1] ]], lp[[ matnames[1,2] ]], lp[[ matnames[1,3] ]], lp[[ matnames[1,4] ]],
    textGrob(vlinhas[1], rot = 270, gp = gpar(fontsize = 12, fontface = "bold")),
    ## 2a linha
    lp[[ matnames[2,1] ]], lp[[ matnames[2,2] ]], lp[[ matnames[2,3] ]], lp[[ matnames[2,4] ]],
    textGrob(vlinhas[2], rot = 270, gp = gpar(fontsize = 12, fontface = "bold")),
    ## 3a linha
    lp[[ matnames[3,1] ]], lp[[ matnames[3,2] ]], lp[[ matnames[3,3] ]], lp[[ matnames[3,4] ]],
    textGrob(vlinhas[3], rot = 270, gp = gpar(fontsize = 12, fontface = "bold"))
  )
  # Plotar tudo
  p <- arrangeGrob(
    grobs = elementos,
    layout_matrix = layout_matrix,
    heights = heights,
    widths = widths
  )
}
df_plot <- df_real %>% 
  select(efeito,Uefeito,k,SiteCode,p_class) %>%
  pivot_wider(names_from="efeito",values_from="Uefeito")
p <- f_figgeral(df_plot)
ggsave(filename="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/fragperse_e_areaperse_exploracaofinal.png",
       plot=p,
       width = 12,
       dpi=200)
############################### 1a versão
f_finalNAOUSADA <- \(dfdados){
  df_plot <- dfdados %>% 
    filter(efeito!="Frag. total") %>% 
    select(-fit,-se.fit,-lower,-upper) %>% 
    pivot_wider(names_from = "efeito",values_from="Uefeito") %>% 
    mutate(A_maior0 = ifelse(`Área per se`>=0,"Área per se>=0","Área per se<0"),
           A_maior0 = factor(A_maior0,levels=c("Área per se<0","Área per se>=0")),
           F_maior0 = ifelse(`Frag. per se`>=0,"Frag per se>=0","Frag per se<0"),
           F_maior0 = factor(F_maior0,levels=c("Frag per se>=0","Frag per se<0")),
           abs_AmaiorF = case_when(
             abs(`Área per se`) > abs(`Frag. per se`) ~ "|Área per se|>|Frag per se|",
             abs(`Área per se`) < abs(`Frag. per se`) ~ "|Frag per se|>|Área per se|",
             abs(`Área per se`) == abs(`Frag. per se`) ~ "|frag|=|area|"))
  dfp <- df_plot %>% 
    group_by(A_maior0,F_maior0,abs_AmaiorF) %>% 
    summarise(n=n(),
              nSite=length(unique(SiteCode))) %>% 
    group_by(A_maior0,F_maior0) %>% 
    mutate(perc = round(n*100/sum(n),2),
           perc_sites = round(nSite*100/67,2))  
  #####
  p_reg <- f_scatterplot(dfp = df_plot,xrange = xrange,yrange = yrange)
  p_circ <- f_circle_plot(dfp,boxlabelsize=3.5,tsize=6,lsize = 8)
  p_hist <- f_hist_ggplot(dff = dfdados,tsize = 15)
  design <- "AB
             AB
             AC
             AC
             AC"
  p_final <- wrap_plots(A=p_reg,B=p_circ,C=p_hist,design = design,widths = c(2,1.5))
  return(p_final)
}
