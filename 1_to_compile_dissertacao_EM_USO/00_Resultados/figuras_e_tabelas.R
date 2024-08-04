v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
################################
# criacao GE dados disponiveis #
################################
theme_set(theme_light())
df_dados_disponiveis <- read_csv(file = "dados/df_dados_disponiveis.csv")
df_plot <- df_dados_disponiveis %>% 
  inner_join(df_p) %>% 
  select(SiteCode, p, effort_ha, Ntotal, S_obs, year_bestProxy, forest_succession, lat:long_correct)
df_p_extensoes <- read_csv("dados/csv/df_p_extensoes.csv")
l_p <- list()
l_p[[1]] <- df_plot |> 
  ggplot(aes(y=lat,x=long,fill=p)) +
  geom_point(alpha=0.85,shape = 21,size=3,color="black") +
  scale_fill_gradientn("% CF\n4 km",colours = terrain.colors(10)) +
  guides(fill=guide_colourbar(position="inside")) +
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
  ggplot(aes(x=forest_succession)) +
  geom_bar() +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.2, colour = "black") +
  labs(x="classificação de perturbação (tempo e grau)",
       y="contagem",subtitle="d)")
# grid.arrange(grobs=l_p,nrow=2)
p <- arrangeGrob(grobs=l_p,nrow=2)
ggsave(
  filename=paste0(v_path,
                  "figuras/GE_dados_disponiveis.png"),
  plot = p,
  width=12,height = 10)
################################
#
################################
##### criação de fig5_SoE ######
################################
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
################################
#
##################################
# criação de GE_taxaU_contrastes #
##################################
df_contrastes <- read_csv(file="dados/csv/taxaU/df_contrastes.csv") |> 
  inner_join(df_sim |> select(SiteCode,Ntotal:S_obs) |> distinct(),
             by="SiteCode") |> 
  rename(N=Ntotal,S=S_obs) |> 
  mutate(across(N:S,log,.names="log.{.col}"),
         across(c(p,k,log.N:log.S),f_z,.names = "{.col}_z"),
         SiteCode = factor(SiteCode)) |> 
  select(-c(N:log.S))
l_p <- list()
df_plot <- df_contrastes %>% select(SiteCode:p, contains("_logratio")) %>% 
  pivot_longer(-c(SiteCode:p)) %>% 
  mutate(name = gsub("_logratio","",name))
v_sites_RefNulo <- df_plot %>% filter(p>=0.975) %>% pull(SiteCode) %>% unique
v_range_RefNulo <- df_plot %>% filter(p>=0.975) %>% pull(value) %>% range
theme_set(theme_bw())
df_text <- df_plot %>% 
  mutate(maior_q = value > max(v_range_RefNulo),
         menor_q = value < min(v_range_RefNulo)) %>% 
  group_by(name) %>% 
  summarise(prop_maiorq = round(100 * sum(maior_q) / n(),2),
            prop_menorq = round(100 * sum(menor_q) / n(),2)) %>% 
  mutate(v_geomtext = case_when(grepl("area",name) ~ "Contraste Área:\n",
                                grepl("frag.perse",name) ~ "Contraste Frag. per se:\n",
                                grepl("frag.total",name) ~ "Contraste Frag. Total:\n"),
         v_geomtext = paste0(v_geomtext,prop_maiorq,"% a direita ",prop_menorq,"% a esquerda"))
l_p$fig1 <- df_plot %>% 
  mutate(
    label = case_when(
      grepl("area",name) ~ "Contraste Área per se: log(U aglomerado / U prístino)",
      grepl("total",name) ~ "Contraste Frag. Total: log(U contemporâneo / U prístino)",
      grepl("perse",name) ~ "Contraste Frag. per se: log(U contemporâneo / U aglomerado)")) %>%
  left_join(.,df_text) %>%
  ggplot(aes(x=value)) +
  geom_vline(xintercept = 0,color="darkblue",linetype=3) +
  geom_vline(xintercept = v_range_RefNulo,color="darkred") + 
  geom_histogram(bins = 120) +
  geom_boxplot(aes(y=-17.5),width=30,alpha=0.4) +
  labs(y="",x="",
       tag="a)",
       caption="contraste em paisagens com %CF>=0.95") +
  geom_text(aes(x=1.5,y=300,label=v_geomtext),alpha=0.01,size=4) + 
  theme(plot.margin = unit(c(0.1,0.25,0.25,0), "cm"),
        strip.text = element_text(size=10),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0.2, face= "italic",vjust=10)) +
  facet_wrap(~label,ncol=1)
##
df_md <- df_contrastes %>% select(SiteCode:p, contains("_logratio")) %>% 
  pivot_longer(-c(SiteCode:p)) %>% 
  mutate(name = gsub("_logratio","",name),
         across(c(p,k),f_z,.names = "{.col}_z"),
         SiteCode = factor(SiteCode))
l_p$fig2 <- df_md %>% 
  mutate(label=case_when(
    grepl("area",name) ~ "Área per se",
    grepl("total",name) ~ "Frag. Total",
    grepl("perse",name) ~ "Frag. per se")
  ) %>% 
  ggplot(aes(x=k,y=value,group=SiteCode,color=p)) +
  geom_boxplot(inherit.aes = FALSE,
               aes(x=k,y=value,group=k)) +
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  scale_colour_gradient2("% CF",midpoint=0.5,
                         low="darkred",
                         mid = "blue",
                         high = "darkgreen") +
  labs(x="k (prop. de propágulos até a vizinhança imediata)",
       y="log(U/U)",
       tag="b)") +
  facet_wrap(~label,ncol=3)
# grid.arrange(grobs=l_p,ncol=1)
p <- arrangeGrob(grobs=l_p,
                 # ncol=1,
                 layout_matrix = rbind(c(1),
                                       c(1),
                                       c(1),
                                       c(2),
                                       c(2)))
ggsave(
  filename = paste0(v_path,"figuras/GE_taxaU_contrastes.png"),
  plot = p,
  width = 10,height=7)
