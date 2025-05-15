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
l_dfpred <- readRDS(file = "dados/csv_SoE/rds/l_dfpred_md_cong_absoluta.rds")
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
p <- df_logUU_pk %>%
  mutate(across(c(SiteCode,contraste),factor),
         label = factor(contraste,
                        levels=c("Frag. total",
                                 "Frag. per se",
                                 "Área per se"))) %>% 
  ggplot(aes(x=k,y=Uefeito,group=SiteCode,color=p)) +
  geom_hline(yintercept = 0,color="black",linetype=3) +
  # geom_hline(yintercept = v_hline,color="darkred") + 
  geom_line(alpha=0.75) +
  geom_point(alpha=0.75) +
  geom_boxplot(inherit.aes = FALSE,
               aes(x=k,y=Uefeito,group=k),alpha=0.25) +
  scale_colour_gradient2("% CF",midpoint=0.5,
                         low="red",
                         mid = "yellow",
                         high = "darkgreen") +
  labs(x="Proporção de propágulos na vizinhança imediata (k)",
       y="log(U/U)") +
  scale_y_continuous(expand=c(0.01,0.01)) +
  theme_classic() +
  theme(plot.margin=unit(c(0,0.2,0,0), "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.49,0.92),
        legend.direction="horizontal") +
  facet_grid(pert_class~label)
ggsave("figuras/pedacofigfinal_1alinha.png",
       p,
       dpi=300,
       width = 13.8,
       height = 7.33)
img_obj <- image_read("figuras/pedacofigfinal_1alinha.png") %>% 
  image_trim() %>% 
  image_resize("75%")
image_write(img_obj,"figuras/pedacofigfinal_1alinha.png")
#
#
#########################################################################################
########################### logU/U : área per se e fragmentçaão per se ##################
#########################################################################################
f_gsub <- \(xlab){
  gsub("area","Área per se",xlab) %>% 
    gsub("fragperse","Frag. per se",.) %>% 
    gsub("fragtotal","Frag. total",.)
}

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
    geom_hex(bins = 50,alpha=0.5),
    scale_fill_gradient("contagem",low = "yellow", high = "red", na.value = NA),
    facet_wrap(~pert_class,ncol=3),
    theme_classic()
  )
  df_plot %>% 
    pivot_wider(names_from=contraste,values_from=Uefeito) %>%   
    ggplot(aes(x=.data[[xlab]],y=.data[[ylab]])) +
    geom_list1 +
    labs(title=paste0("logU/U :",label_efeitos),
         x=f_gsub(xlab),
         y=f_gsub(ylab))
}
l_fplot$tabela_contagem_maior_magnitude_absoluta <- \(dfi){
  dfi %>% 
    group_by(pert_class,class_diff) %>% 
    tally() %>% 
    group_by(pert_class) %>% 
    mutate(n_rel = round(n * 100 / sum(n),2)) %>% 
    ggplot(aes(y=class_diff,x=pert_class,fill=n_rel)) +
    geom_raster() +
    labs(x="",y="",
         title="") +
    scale_fill_gradient("",low = "yellow", high = "red", na.value = NA) +
    geom_text(aes(label=n_rel),size=10) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand =  c(0,0)) +
    theme_classic() +
    facet_wrap(~pert_class,ncol=3,scales = "free_x") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=15),
          axis.ticks.y = element_blank(),
          title = element_text(size=15),
          strip.text = element_text(size=15),
          legend.position = "none")
}
l_fplot$histograma_efeito_menos_efeito <- \(dfi,xlabel=NULL){
  if(is.null(xlabel)) stop("precisa fornecer xlabel = e.g. logU/U: |área per se| - |frag. per se|")
  geom_list2 <- list(
    geom_histogram() ,
    geom_vline(xintercept = 0,color="red",linetype=2) ,
    theme_classic() ,
    scale_x_continuous(expand=c(0,0)) ,
    scale_y_continuous(expand=c(0,0)) ,
    facet_wrap(~pert_class,ncol=3)
  )
  dfi %>% 
    ggplot(aes(x=diff_efeito)) +
    geom_list2 +
    labs(x=xlabel,
         y="") 
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
#########################################################################################
########################### classificação das paisagens hipotéticas quanto ao número de boas congruências
df_ad <- read_csv(file="dados/csv_SoE/df_congruencia_simulacao.csv") %>% 
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
                                       "muito baixa\n<=5")[5:1])) %>% 
  filter(k>=0.49999)
df_plot <- df_ad %>% 
  group_by(land,congruencia) %>% 
  tally()
df_plot %>% 
  ggplot(aes(x=land,y=congruencia,fill=n)) +
  geom_tile(color="white") +
  geom_label(aes(label=n),fill="white",size=10) +
  labs(x="",y="",title="Núm. SADs cong.: k>=0.50") +
  scale_x_discrete(expand=c(0,0),position = "top") +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_viridis_c("n") +
  theme_classic() +
  theme(text=element_text(size=15,face="bold"))

df_ij <- df_ad %>% 
  group_by(SiteCode) %>% 
  summarise(nSAD_maiorigual75 = all(nSAD>=75))
v_ij <- c("alta congruência sempre"=sum(df_ij$nSAD_maiorigual75),
          "baixa congruência em algum k"=nrow(df_ij)-sum(df_ij$nSAD_maiorigual75))

df_logUU_plot <- inner_join(df_logUU_pk,df_ij) %>% 
  mutate(label = ifelse(nSAD_maiorigual75,
                        paste0("alta cong. sempre",", n sítios=",v_ij[1]),
                        paste0("cong. baixa em algum k",", n sítios=",v_ij[2])))
p1 <- df_logUU_plot %>% 
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
  theme_classic()
p2 <- df_logUU_plot %>% 
  mutate(label = gsub(", n sítios\\=","",label) %>% 
           gsub("\\d+","",.) %>% 
           gsub(" sempre","",.) %>% 
           gsub(" em algum k","",.)) %>% 
  filter(k==0.99,contraste=="Frag. total") %>% 
  ggplot(aes(fill=label,x=p)) +
  geom_histogram(position = "dodge",bins=60) +
  labs(x="%CF",y="contagem") +
  scale_fill_manual("",values=c("red","blue")) +
  guides(fill=guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_minimal_grid() +
  theme(legend.position = "top",
        text = element_text(size=8,face="bold"),
        axis.text = element_text(size=8),
        legend.margin = margin(),
        plot.margin = margin())
library(cowplot)
p <- ggdraw() +
  draw_plot(p1) +
  draw_plot(p2,
            height = 0.27, width = 0.25,
            x=0.66, y=0.35)
ggsave(filename="1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/pedacofigfinal_1alinha.png",
       plot=p,
       width = ,
       height = ,
       dpi=200)
df_md <- select(df_logUU_plot[df_logUU_plot$nSAD_maiorigual75,],
                SiteCode:Uefeito,p)
saveRDS(df_md,"dados/csv_SoE/rds/df_logUUpk.rds")


############################### anterior ao 3o comitê
# acordamos de remover a parte do logOR! E dessa forma as coisas ficaram bem mais focadas e claras!!!





############################################################################################
############################## diagnóstico dos modelos mais plausíveis 2#####################
############################################################################################
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
# modelos ajustados
l_path <- list()
l_path$te <-  paste0(v_path,"rds/l_md_",c("areaperse","fragperse","fragtotal"),".rds")
l_rds <- gsub("l_md_","l_dfpred_",l_path$te)
names(l_rds) <- str_extract(l_rds,"(?<=dfpred\\_)(.*?)(?=\\.rds)") %>% 
  gsub("areaperse","area",.)
l_dfplot <- lapply(l_rds,readRDS) %>% 
  lapply(.,"[[","fixo e aleat") %>% 
  lapply(.,select, logOR:k_z, SiteCode, Q_0.5, forest_succession)
#####
df_logUU_plot <- df_logUU_pk %>% select(SiteCode:k,k_z,contraste,Uefeito) %>% 
  mutate(contraste = gsub("Área per se","area",contraste) %>% 
           gsub("Frag. per se","fragperse",.) %>% 
           gsub("Frag. total","fragtotal",.)) %>% 
  pivot_wider(names_from="contraste",values_from="Uefeito")
#####
df_plot <- lapply(names(l_dfplot),\(li){
  dfi <- select(l_dfplot[[li]],-Uefeito)
  inner_join(dfi,
             select(df_logUU_plot,SiteCode:k_z,all_of(li)),
             by=c("SiteCode","k_z")) %>% 
    rename("Uefeito" = li) %>% 
    mutate(efeito = li,
           pert_class = factor(forest_succession,
                               levels=c("primary",
                                        "primary/secondary",
                                        "secondary"),
                               labels=c("baixa",
                                        "mediana",
                                        "alta"))) %>% 
    select(-forest_succession)
}) %>% do.call("rbind",.)
f_lmbySite <- \(dfs){
  md <- lm(Q_0.5 ~ logOR,data=dfs)
  coef <- unname(md$coefficients)
  data.frame(slope=coef[2],
             intercept=coef[1])
}
df_lmpars <- ddply(df_plot,c("efeito","SiteCode"),f_lmbySite)
df_plot2 <- inner_join(df_plot,df_lmpars,by=c("SiteCode","efeito"))
f_plot <- \(df_plot){
  #@ xylab = "logOR","logU/U"
  geom_list1 <- list(
    geom_hline(yintercept = 0,color="black"),
    geom_vline(xintercept = 0,color="black"),
    geom_abline(slope=1,intercept=0,color="darkblue",linetype=1),
    geom_abline(slope=-1,intercept=0,color="darkblue",linetype=1),
    geom_hex(bins = 50,alpha=0.5),
    scale_fill_gradient("contagem",low = "yellow", high = "red", na.value = NA),
    facet_grid(pert_class~label),
    theme_classic(),
    scale_x_continuous(expand=c(0,0)),
    scale_y_continuous(expand=c(0,0))
  )
  df_plot <- df_plot %>% 
    mutate(label = factor(efeito,
                          levels=c("fragtotal",
                                   "fragperse",
                                   "area"),
                          labels=c("Frag. total",
                                   "Frag. per se",
                                   "Área per se")))
  
  l_p <- list()
  l_p[[1]] <- df_plot %>% 
    ggplot(aes(x=Uefeito,y=logOR)) + # ,color=intercept,group=SiteCode
    #geom_point(alpha=0.1) +
    #geom_line(alpha=0.1) +
    geom_list1 +
    labs(x="logU/U",
         y="logOR")
  df_plot <- mutate(df_plot,diff_logs=abs(logOR)-abs(Uefeito))
  geom_list2 <- list(
    geom_histogram() ,
    geom_vline(xintercept = 0,color="red",linetype=2) ,
    theme_classic() ,
    scale_x_continuous(expand=c(0,0)) ,
    scale_y_continuous(expand=c(0,0)) ,
    facet_grid(pert_class~label)
  )
  l_p[[2]] <- df_plot %>% 
    ggplot(aes(x=diff_logs)) +
    geom_list2 +
    labs(x="|logOR| - |logU/U|",y="")
  l_p[[3]] <- df_plot %>% 
    # pivot_longer(c(logOR,Uefeito)) %>% 
    # mutate(name=gsub("Uefeito","logU/U",name)) %>% 
    ggplot(aes(x=slope,y=logOR)) +
    geom_vline(xintercept = 1,color="black",alpha=0.25) +
    # geom_point(alpha=0.1) +
    geom_hex(bins = 50,alpha=0.75) +
    scale_fill_gradient("contagem",low = "darkblue", high = "red", na.value = NA) +
    labs(x="inclinação: lm(obs ~ pred)",y="logOR") +
    # scale_color_manual("",values=c("black","red")) +
    facet_grid(pert_class~label)
  l_p[[4]] <- ggplot(df_plot,aes(x=logOR,y=Q_0.5,color=slope)) +
    geom_point() +
    scale_colour_gradient2("inclinação",midpoint=1,
                           low="red",
                           mid = "blue",
                           high = "green") +
    facet_grid(pert_class~label)
  names(l_p) <- c("logOR ~ loU/U","histograma: logOR - logU/U","logOR ~ inclinação","logOR_pred ~ logOR_obs")
  return(l_p)
}
l_p <- f_plot(df_plot2)
saveRDS(l_p,file="./1_to_compile_dissertacao_EM_USO/09_SI/figuras_SI/l_p_comparacao_empirica_logOR_logUU.rds")
###############################################################################################
####################################### diagnóstico final artigo

## funções de ajuste e de plot
source("source/2samples_testes.R")
source("source/general_tools.R")
source("source/GAMMtools.R")
source("source/fig_tools.R")
f_imagefunc <- \(ggobj){
  tempfile_path <- tempfile(fileext = ".png")
  ggsave(tempfile_path,ggobj,width = 7.17,height = 7.20)
  img <- image_trim( image_read(tempfile_path))
  image_write(img,tempfile_path)
}
f_resize_2rectangle <- \(ref_img,toresize_img,ref_side,tore_side){
  tore_info <- image_info(toresize_img)
  ref_info <- image_info(ref_img)
  refside_info <- ref_info[[ref_side]]
  v_command <- ifelse(tore_side=="height",
                      paste0("x",refside_info),
                      paste0(refside_info,"x"))
  image_resize(toresize_img,v_command)
}
image_title <- \(imgobj,
                 vtitle,
                 vheight=150,
                 bgcolor="white",
                 vsize=80,
                 vgrav="north",
                 vcolor="black",
                 vloc="+0+20"){
  library(magick)
  image_info <- image_info(imgobj)
  image_width <- image_info$width
  white_canvas <- image_blank(
    width = image_width, 
    height = vheight, 
    color = bgcolor)
  imgobj <- image_append(c(white_canvas, 
                           imgobj), 
                         stack = TRUE)
  image_annotate(
    imgobj,
    text = vtitle,
    size = vsize,
    gravity = vgrav,
    color = vcolor,
    location = vloc
  )# %>% image_trim
}
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/dados/csv_SoE/"
l_path <- paste0(v_path,"rds/l_dfpred_",
                 c("areaperse",
                   "fragperse",
                   "fragtotal"),
                 ".rds")
l_df <- lapply(l_path,readRDS)
names(l_df) <- gsub("areaperse","Área per se",
                    str_extract(l_path,"(?<=dfpred_)(.*?)(?=\\.rds)")) %>% 
  gsub("fragperse","Frag. per se",.) %>% 
  gsub("fragtotal","Frag. total",.)
l_df <- lapply(l_df,"[[","fixo e aleat")
######### gráficos diagnósticos para o artigo principal
library(hexbin)
f_plot <- \(dfp,bypert=FALSE){
  dfp <- dfp %>% 
    pivot_longer(starts_with("Q_"),
                 names_to="quantile",
                 values_to = "preditor") %>% 
    mutate(quantile=gsub("Q_","",quantile) %>% 
             as.numeric(),
           quantile=as.character(quantile*100))
  fggplot <- \(dfi){
    df_slopes <- ddply(dfi,"SiteCode",\(dfii){
      mdl <- as.data.frame(lm(preditor~logOR,data=dfii)$coefficients)
      mdl$coef = rownames(mdl) %>% gsub("\\(","",.) %>% gsub("\\)","",.)
      rownames(mdl) <- NULL
      names(mdl)[1] <- "valor"
      pivot_wider(mdl,names_from="coef",values_from="valor")
    }) %>% pivot_longer(-SiteCode) %>% 
      mutate(name=ifelse(name=="logOR","preditor",name),
             vline = ifelse(name=="preditor",1,0),
             name = gsub("Intercept","intercepto",name) %>% 
               gsub("preditor","inclinação",.),
             name=factor(name,levels=c("intercepto","inclinação")))
    l_p <- list()
    l_p$coef <- ggplot(df_slopes,aes(x=value)) +
      geom_histogram(aes(y=after_stat(density)),bins=50,fill="skyblue",color="black") +
      geom_density(color="darkgreen",linewidth=1) +
      geom_vline(aes(xintercept=vline),color="red",linetype=2,linewidth=1) +
      scale_y_continuous("densidade",expand=c(0,0)) +
      xlab("estimativa") +
      facet_wrap(~name,ncol=1,scales="free") +
      theme(strip.text = element_text(size=13,margin = margin()),
            plot.margin = margin(0,0,0,0))
    l_p$predobs <- dfi %>% 
      mutate(quantile=paste0("Quantil: ",quantile,"%")) %>% 
      ggplot(aes(x=logOR,y=preditor)) +
      geom_smooth(aes(group=SiteCode),
                  method="lm",
                  se=FALSE,
                  color="#4682B4",
                  alpha=0.5) +
      geom_hex(bins = 50,alpha=0.5) + 
      geom_abline(slope=1,intercept=0,color="black",linewidth=1.2,linetype=2) +
      scale_fill_gradient("contagem",low = "yellow", high = "red", na.value = NA) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      labs(x="logOR observado",y="logOR predito") +
      facet_wrap(~quantile) +
      theme(legend.position = c(0.90,0.25),
            strip.text = element_text(size=13,margin = margin()),
            plot.title = element_text(size=10,
                                      hjust=0.5,
                                      face="bold"),
            axis.title.x = element_text(margin=margin(t=0,b=0,l=0,r=0)),
            axis.title.y = element_text(margin=margin(t=0,b=0,l=0,r=0))
      )
    if(dfi$quantile[1]=="50"){
      p <- ggdraw() +
        draw_plot(l_p[["predobs"]]) +
        draw_plot(l_p[["coef"]],
                  height = 0.4, width = 0.25,
                  x=0.08, y=0.55)
      vpath <- f_imagefunc(p)
      return(vpath)
    }else{
      p <- l_p$predobs +
        theme(plot.tile=element_text(size=))
      vpath <- f_imagefunc(l_p$predobs)
      return(vpath)
    }
  }
  if(bypert){
    dfpath <- ddply(filter(dfp,quantile=="50"),"forest_succession",fggplot)
  }else{
    dfpath <- fggplot(filter(dfp,quantile=="50"))
  }
  # lp <- dlply(dfp,"quantile",fggplot)
  # limg <- lapply(lp,\(li) image_trim(image_read(li)))
  # base_img <- image_append(do.call("c",limg[c("5","95")]),stack=FALSE)
  # base_img <- f_resize_2rectangle(ref_img = limg[["50"]],
  #                                 toresize_img = base_img,
  #                                 ref_side = "width",
  #                                 tore_side = "width")
  # img_final <- image_append(c(limg[["50"]],base_img),stack = TRUE) %>% 
  #   image_resize("50%")
  l_img <- alply(dfpath,1,\(li){
    vpath <- ifelse(is.data.frame(li),li$V1,li)
    image_read(vpath) %>% 
      image_trim() %>% 
      image_resize("50%")
  })
  if(is.data.frame(dfpath)){
    names(l_img) <- gsub("^primary$","baixa",dfpath$forest_succession) %>% 
      gsub("^secondary","alta",.) %>% 
      gsub("primary/secondary","mediana",.)
    l_img <- lapply(names(l_img),\(li){
      image_title(imgobj = l_img[[li]],
                  vtitle = paste0("classe de pert.: ",li),
                  vsize = 60,
                  vgrav = "north",
                  vheight = 120)
    })
    img_final <- image_append(do.call("c",l_img),stack = TRUE)
    return(img_final)
  }else{
    img <- image_read(dfpath) %>% 
      image_trim() %>% 
      image_resize("50%")
    return(img)  
  }
}
l_img_diag <- lapply(l_df,f_plot)
l_img_diag <- lapply(names(l_img_diag),\(li){
  img <- l_img_diag[[li]]
  image_title(imgobj = img, 
              vtitle = li,
              vsize = 50,
              vheight = 90) %>% 
    image_trim
})
names(l_img_diag) <- paste0(
  v_path,
  "figuras/diagfinal_",
  case_when(
    grepl("Área",names(l_df)) ~ "areaperse",
    grepl("Frag. per se",names(l_df)) ~ "fragperse",
    grepl("total",names(l_df)) ~ "fragtotal",
  ),
  ".jpeg"
)
lapply(names(l_img_diag),\(li) image_write(l_img_diag[[li]],path=li))
##### colar as figuras
l_img <- lapply(names(l_img_diag),image_read)
names(l_img) <- names(l_img_diag)
vi <- sapply(c("fragtotal","fragperse","areaperse"),\(x) grep(x,names(l_img)))
img_final <- image_append(do.call("c",l_img[vi]),stack = FALSE) %>% 
  image_resize("50%")
image_write(img_final,
            "dados/csv_SoE/figuras/diagfinal_3efeitos.jpeg",
            format="jpeg")
image_write(img_final,"figuras/diagfinal_3efeitos.jpeg",
            format="jpeg")
