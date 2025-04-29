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


