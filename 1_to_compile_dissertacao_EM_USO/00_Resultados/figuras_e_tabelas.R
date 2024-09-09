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
df_p <- read_csv("dados/df_p.csv")
v_path <- "/home/danilo/Documentos/mestrado_Ecologia/artigo_principal/1_to_compile_dissertacao_EM_USO/00_Resultados/"
setwd(v_path)
probs = c(0.05,0.25, 0.5, 0.75,0.95)
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
f_z <- function(x) (x-mean(x))/sd(x)
df_sim <- read_csv("dados/df_simulacao.csv") |> 
  inner_join(x=df_p,by="SiteCode")
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
v_hline <- filter(df_plot,p>=0.95) %>% 
  pull(value) %>% quantile(.,probs)
theme_set(theme_bw())
df_text <- df_plot %>% 
  mutate(maior_q = value > max(v_hline),
         menor_q = value < min(v_hline)) %>% 
  group_by(name) %>% 
  summarise(prop_maiorq = round(100 * sum(maior_q) / n(),2),
            prop_menorq = round(100 * sum(menor_q) / n(),2)) %>% 
  mutate(v_geomtext = case_when(grepl("area",name) ~ "Contraste Área:\n",
                                grepl("frag.perse",name) ~ "Contraste Frag. per se:\n",
                                grepl("frag.total",name) ~ "Contraste Frag. Total:\n"),
         v_geomtext = paste0(v_geomtext,prop_maiorq,"% a direita, ",prop_menorq,"% a esquerda"))
# l_p$fig1 <- df_plot %>% 
#   mutate(
#     label = case_when(
#       grepl("area",name) ~ "Contraste Área per se: log(U aglomerado / U prístino)",
#       grepl("total",name) ~ "Contraste Frag. Total: log(U contemporâneo / U prístino)",
#       grepl("perse",name) ~ "Contraste Frag. per se: log(U contemporâneo / U aglomerado)")) %>%
#   left_join(.,df_text) %>%
#   ggplot(aes(x=value)) +
#   geom_vline(xintercept = 0,color="blue",linetype=3) +
#   geom_vline(xintercept = v_hline,color="darkred") + 
#   geom_histogram(bins = 120) +
#   geom_boxplot(aes(y=-17.5),width=30,alpha=0.4) +
#   labs(y="",x="",
#        tag="a)",
#        caption="quantils (0.05, 0.25, 0.75, 0.95) dos contraste em paisagens com %CF>=0.95") +
#   geom_text(aes(x=1.5,y=300,label=v_geomtext),alpha=0.01,size=4) + 
#   theme(strip.text = element_text(size=10),
#         plot.caption.position = "plot",
#         plot.caption = element_text(hjust = 0.2, face= "italic",vjust=10),
#         plot.margin=unit(c(0,0.2,-0.5,0), "cm")) +
#   facet_wrap(~label,ncol=1)
# Distribuição da taxa U em função da perturbação dos sítios
df_plot <- inner_join(
  df_md_Uefeito,
  df_dados_disponiveis %>% select(SiteCode,forest_succession)
) %>% filter(forest_succession!="capoeira") %>% 
  mutate(tp_efeito = gsub("area",
                          "Área per se: aglomerado / prístino",tp_efeito) %>% 
           gsub("frag.perse",
                "Frag. per se: contemporâneo / aglomerado",.) %>% 
           gsub("frag.total",
                "Frag. total: contemporâneo / prístino",.),
         forest_succession = gsub("primary$","1a",forest_succession) %>% 
           gsub("^primary/secondary$","1a-2a",.) %>% 
           gsub("^secondary","2a",.))
df_plot <- rbind(
  df_plot %>% mutate(forest_succession = "todos"),
  df_plot
)
l_p$fig1 <- df_plot %>% 
  ggplot(aes(x=forest_succession,y=v_efeito)) +
  labs(x="classe de preservação (grau de perturbação e tempo de recuperação)",
       y="log(U/U)",
       tag="a)") +
  geom_hline(yintercept = 0,color="blue",linetype=3) +
  geom_hline(yintercept = v_hline,color="darkred") + 
  geom_jitter() +
  geom_boxplot(alpha=0.4) +
  # scale_x_discrete(position = "top") +
  # coord_flip() +
  facet_wrap(~tp_efeito,nrow=1)
##
df_md <- df_contrastes %>% select(SiteCode:p, contains("_logratio")) %>% 
  pivot_longer(-c(SiteCode:p)) %>% 
  mutate(name = gsub("_logratio","",name),
         across(c(p,k),f_z,.names = "{.col}_z"),
         SiteCode = factor(SiteCode))
l_p$fig2 <- df_md %>% 
  mutate(label=case_when(
    grepl("area",name) ~ "Área per se: aglomerado / prístino",
    grepl("total",name) ~ "Frag. Total: contemporâneo / prístino",
    grepl("perse",name) ~ "Frag. per se: contemporâneo / aglomerado")
  ) %>% 
  ggplot(aes(x=k,y=value,group=SiteCode,color=p)) +
  geom_boxplot(inherit.aes = FALSE,
               aes(x=k,y=value,group=k)) +
  geom_hline(yintercept = 0,color="blue",linetype=3) +
  geom_hline(yintercept = v_hline,color="darkred") + 
  geom_line(alpha=0.5) +
  geom_point(alpha=0.5) +
  scale_colour_gradient2("% CF",midpoint=0.5,
                         low="red",
                         mid = "black",
                         high = "#69874c") +
  labs(x="k (prop. de propágulos até a vizinhança imediata)",
       y="log(U/U)",
       tag="b)") +
  theme(plot.margin=unit(c(-0.2,0.2,0,0), "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.49,0.9),
        legend.direction="horizontal") +
  # guides(color=guide_legend(direction='horizontal')) +
  facet_wrap(~label,ncol=3)
# grid.arrange(grobs=l_p,ncol=1)
p <- arrangeGrob(grobs=l_p,
                 ncol=1,
                 # layout_matrix = rbind(c(1),
                 #                       c(1),
                 #                       c(1),
                 #                       c(2),
                 #                       c(2))
                 )
ggsave(
  filename = paste0(v_path,"figuras/GE_taxaU_contrastes.png"),
  plot = p,
  width = 11,height=7.7)
#
# DAG do conjunto de variáveis de controle:
library(dagitty)
library(ggdag)
library(ggplot2)
library(dplyr)
cov_dag <- dagify(
  local ~ paisagem + biogeo + preservacao_local,
  paisagem ~ biogeo + U,
  preservacao_local ~ U,
  labels = c(
    "local" = "biodiv. local",
    "paisagem" = "efeito paisagem",
    "preservacao_local" = "preserv. local",
    "biogeo" = "contexto biogeo.",
    "U" =  "unobs. cov."
  ),
  latent = "U",
  exposure = "paisagem",
  outcome = "local",
  coords = list(
    x = c(paisagem = 2.5, biogeo = 7.5, preservacao_local = 9, local = 5, U = 9),
    y = c(paisagem = 5, biogeo = 3.5, preservacao_local = 2.5, local = 1, U = 5.5)
  )
)
# cov_dag %>% 
#   ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
#         geom_dag_point(color = "orange") +
#         geom_dag_edges_arc(edge_color = "blue", curvature = 0) +
#         geom_dag_label_repel(aes(label = label), colour = "black", show.legend = FALSE) +
#         theme_dag()
# Tidy the DAG and set the node class
tidy_cov_dag <- tidy_dagitty(cov_dag) %>%
  mutate(node_class = case_when(
    name == "paisagem" ~ "exposure",
    name == "local" ~ "outcome",
    name == "U" ~ "latent",
    TRUE ~ "adjusted"
  ))
# Create a new column to classify edge type
tidy_cov_dag <- tidy_cov_dag %>%
  mutate(edge_type = case_when(
    name %in% c("biogeo", "preservacao_local") ~ "closed path",
    TRUE ~ "open path"
  ))
# Split edge data based on parent node for custom coloring
edges_gray <- filter(tidy_cov_dag$data, name %in% c("biogeo", "preservacao_local"))
edges_blue <- filter(tidy_cov_dag$data, !name %in% c("biogeo", "preservacao_local"))
# criação de objeto para inclusão da legenda
legend_data <- data.frame(
  x = c(2.5, 2.5),
  y = c(2.5, 2.5),
  xend = c(2.5, 2.5),
  yend = c(2.5, 2.5),
  edge_type = c("closed path", "open path")
)
# Plot the DAG with custom arrow colors
p <- ggplot(data = tidy_cov_dag) +
  geom_dag_edges_arc(data = edges_gray, 
                     mapping = (aes(x=x,y=y,xend=xend,yend=yend)),
                     curvature = 0, edge_color = "gray") +
  geom_dag_edges_arc(data = edges_blue, 
                     mapping = (aes(x=x,y=y,xend=xend,yend=yend)),
                     curvature = 0, edge_color = "black") +
  geom_dag_point(aes(x = x, y = y, fill = node_class, size = node_class), 
                 shape = 21) +
  geom_dag_label_repel(aes(x = x, 
                           y = y, 
                           label = label), 
                       colour = "black", 
                       show.legend = FALSE,
                       nudge_x = 0.1,
                       nudge_y = 0.5,
                       box.padding = 0.5,
                       point.padding = 0.5) +
  geom_segment(data = legend_data,
               aes(x = x, y = y, xend = xend, yend = yend, color = edge_type),
               arrow = arrow(length = unit(0.000001, "inches"), type = "closed"),
               size = 1) +
  scale_color_manual("Edge Type", values = c(
    "closed path" = "gray", 
    "open path" = "black")) +
  scale_fill_manual("",values = c(
    "exposure" = "darkblue", 
    "outcome" = "darkgreen", 
    "latent" = "black", 
    "adjusted" = "darkred")) +
  scale_size_manual("",values = c(
    "exposure" = 10, 
    "outcome" = 10, 
    "latent" = 7, 
    "adjusted" = 7)) +
  guides(color = guide_legend("")) +
  theme_dag_blank()
ggsave(
  filename = paste0(v_path,"figuras/DAG_var_de_controle.png"),
  plot = p,
  width = 11,height=7.7,
  bg = "white")

img <- image_read(paste0(v_path,
                         "figuras/DAG_var_de_controle.png"))
trimmed_img <- image_trim(img)
print(trimmed_img)
image_write(trimmed_img, path = paste0(v_path,"figuras/DAG_var_de_controle.png"))

#
#
# TABELAS #
#
#
# quantis observados dos contrastes da taxa U
probs = c(0.05,0.25, 0.5, 0.75,0.95)
df_contrastes <- read_csv(file="dados/csv/taxaU/df_contrastes.csv") 
df_md_Uefeito <- df_contrastes %>% select(SiteCode:p, contains("_logratio")) %>% 
  pivot_longer(-c(SiteCode:p)) %>% 
  mutate(name = gsub("_logratio","",name),
         across(c(p,k),f_z,.names = "{.col}_z"),
         SiteCode = factor(SiteCode)) %>% 
  rename("tp_efeito"="name","v_efeito"=value)
f_quant <- function(x) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}
df_quantisobs_Uefeito <- df_md_Uefeito %>% 
  inner_join(
    .,
    select(df_plot,SiteCode,forest_succession) %>% 
      filter(forest_succession!="todos") %>%
      distinct)
df_quantisobs_Uefeito <- rbind(
  df_quantisobs_Uefeito,
  df_quantisobs_Uefeito %>% mutate(forest_succession="todos")
)
df_quantisobs_Uefeito <- df_quantisobs_Uefeito %>%  
  group_by(tp_efeito,forest_succession) %>% 
  reframe(f_quant(v_efeito)) %>% 
  pivot_wider(names_from = "quant",
              values_from = "val",
              names_prefix = "Q: ") %>% 
  mutate(tp_efeito = tp_efeito %>% 
           gsub("area","Área per se",.) %>% 
           gsub("frag.perse","Frag. per se",.) %>% 
           gsub("frag.total","Frag. total",.),
         across(starts_with("Q:"),exp,.names = "exp({.col})")) %>% 
  rename("Contraste log(U/U)"="tp_efeito",
         "Grau de preservação" = "forest_succession")
v_cols <- lapply(probs,\(i) grep(i,names(df_quantisobs_Uefeito),value = TRUE)) %>% 
  do.call("c",.) %>% c(names(df_quantisobs_Uefeito)[1:2],.)
df_quantisobs_Uefeito <- select(df_quantisobs_Uefeito,all_of(v_cols))
write_csv(df_quantisobs_Uefeito,
          file=paste0(v_path,"tabelas/df_quantisobs_Uefeito.csv"))
##
# tabela de comparação de estruturas hierarquicas do modelo
df_tabsel <- read_csv("rds/tabsel_simples.csv")
# df_tabelaSelecao <- read_csv(file=paste0(v_path,"tabelas/df_tabelaSelecaologOR_cgi.csv"))
# df_tabelaSelecao <- df_tabelaSelecao %>% 
#   mutate(`Prática Ontológica` = ifelse(grepl("contemp-ideal",pair),
#                                        "Interdependete","Independente"),
#          pair = pair %>% 
#            gsub("contemp-ideal","Frag. total",.) %>%
#            gsub("contemp-non_frag","Frag. per se",.) %>%
#            gsub("non_frag-ideal","Área per se",.),
#          modelo = modelo %>% 
#            gsub("por paisagem$","por paisagem 1",.) %>% 
#            gsub("1","- gs",.) %>% 
#            gsub("2","- gi",.) %>% 
#            gsub("com ","",.) %>% 
#            gsub("por paisagem","por sítio",.) %>% 
#            gsub("de paisagem","fixo",.)) %>% 
#   select(-`Prática Ontológica`) %>% 
#   rename(Contraste = pair,
#          `Dev. explained` = "dev.expl") %>% 
#   as.data.frame
f_tabelaselecao_com_plot0 <- \(dfi,
                              folder_pattern="figuras/tabsel_c_plot_"){
  v_subtitle <- dfi$contraste[1] %>% 
    gsub("Frag. total","contemp. / príst.",.) %>% 
    gsub("Frag. per se","contemp. / aglomer.",.) %>% 
    gsub("Área per se","aglomer. / príst.",.)
  table <- dfi %>%
    select(modelo:dAICc) %>% 
    gt() %>%
    # tab_header(title = md(dfi$contraste[1]),
    #            subtitle = v_subtitle) %>%
    fmt_number(
      columns = "dAICc",decimals = 2
    ) %>%
    tab_options(
      table.font.size = "small",
      table.align = "left"
    ) %>% 
    cols_label(
      modelo = "GAHM",
      dAICc = "ΔAICc"
    )
  plot <- dfi %>% 
    mutate(rank = 1:n()) %>% 
    ggplot(aes(x = rank)) +
    geom_bar(aes(y = weight), stat = "identity", fill = "gray") +
    geom_point(aes(y = `p-value`), shape = 15, size = 5) + #fill = `Moran I statistic (res)`), 
    geom_segment(
      aes(x = rank - 0.4, xend = rank + 0.4, 
          y = dev.expl, yend = dev.expl,
          color="Deviance\nExplained"), 
      size = 1.5) +
    geom_hline(yintercept = 0.05,alpha=0.2,color="darkgreen") +
    scale_y_continuous("weight, deviance explained, p-value",
                       limits = c(0, 1),
                       expand=c(0,0),
                       breaks = c(0.05,0.25,0.50,0.75,1.00)) +
    # scale_fill_gradientn(colours = c("cyan", "black", "red"),
    #                      values = c(-1,0,1),
    #                      name = "Moran's I\nStatistics") +
    scale_color_manual(values = c("Deviance\nExplained" = "blue"), name = "") +
    labs(x = "", fill = "Moran I statistic",) +
    scale_x_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 90),
      legend.position = "none"
    ) #+ coord_flip()
  # saving the objects
  gtsave(table, paste0(v_path,"tabelas/table_reciclagem.png"),
         vwidth = 800, vheight = 200)
  ggsave(paste0(v_path,"figuras/plot_reciclagem.png"),
         plot = plot,bg = "white",
         width = 3, height = 3, units = "in", dpi = 300)
  img <- image_read(paste0(v_path,"figuras/plot_reciclagem.png")) %>% 
    image_rotate(90)
  # trimmed_img <- image_trim(img) %>% 
    # image_rotate(90)
  image_write(trimmed_img, path = paste0(v_path,"figuras/plot_reciclagem.png"))
  # combinação dos dois:
  # Load the images
  table_img <- image_read(paste0(v_path,"tabelas/table_reciclagem.png"))
  plot_img <- image_read(paste0(v_path,"figuras/plot_reciclagem.png"))
  plot_img <- image_resize(plot_img, geometry = paste0(round(image_info(table_img)$height*0.95,0), "x"))
  height_diff <- image_info(table_img)$height - image_info(plot_img)$height
  padding <- image_blank(width = image_info(plot_img)$width, height = height_diff, color = "none")
  plot_img <- image_append(c(padding, plot_img), stack = TRUE)
  combined_img <- image_append(c(table_img,plot_img), stack = FALSE)
  # Save the combined image
  v_name <- case_when(
    grepl("Total",v_title) ~ "fragtotal.png",
    grepl("Frag. per",v_title) ~ "fragperse.png",
    grepl("Área per",v_title) ~ "areaperse.png",
  )
  image_write(combined_img, path = paste0(v_path,folder_pattern,v_name), format = "png")
}
f_tabelaselecao_com_plot <- \(dff,
                              group_by="Contraste",
                              path_name = paste0(v_path,"figuras/tabsel_aleat.png")){
  folder_pattern <- formals(f_tabelaselecao_com_plot0)[[2]] %>% 
    gsub("figuras/","",.)
  d_ply(dff,group_by,f_tabelaselecao_com_plot0)
  paths <- list.files(path = paste0(v_path,"figuras"),
                      pattern = folder_pattern,
                      full.names = TRUE)
  l_png <- lapply(paths,image_read)
  names(l_png) <- str_extract(paths,"(?<=plot_)(.*?)(?=\\.png)")  
  tabela_final <- image_append(
    do.call("c",l_png),stack = TRUE
  )
  # creating the footnote
  image_write(tabela_final, 
              path = path_name, 
              format = "png")
  file.remove(paths)
  file.remove(paste0(v_path,"figuras/plot_reciclagem.png"))
  file.remove(paste0(v_path,"tabelas/table_reciclagem.png"))
}
f_tabelaselecao_com_plot(df_tabelaSelecao)
# tabela de seleção da comparação de modelos cheios:
library(magick)
f_gt_table <- \(dfi,
                v_name,
                vw = 800,
                vh = 200,
                v_folder="tabelas/table_reciclagem_"){
  v_title <- gsub("Frag. total","Frag. total: contemp. / prist.",v_name) %>%
    gsub("Frag. per se","Frag. per se: contemp. / aglomer.",.) %>%
    gsub("Área per se","Área per se: aglomer. / prist.",.)
  v_names <- names(dfi)
  f_gsub <- \(vchar){
    gsub("rank","Rank",vchar) %>% 
      gsub("modelo","GAHM",.) %>% 
      gsub("df","est. coef.",.) %>% 
      gsub("dAICc","ΔAICc",.) %>% 
      gsub("weight" , "Weight (ΔAICc)",.) %>% 
      gsub("dev.expl" , "Dev. Exp.",.) %>% 
      gsub("Moran I statistic (res)", "Moran's I",.) %>% 
      gsub("p-value" , "p-value",.)
  }
  table <- dfi %>%
    gt() %>%
    tab_header(title = md(v_title)) %>%
    fmt_number(
      columns = v_names[-grep("modelo|rank",v_names)],decimals = 2
    ) %>%
    tab_options(
      table.font.size = "small",
      table.align = "left"
    ) %>% 
    cols_label(
      modelo = "GAHM",
      dAICc = "ΔAICc",
      df = "est. coef.",
      weight = "weight",
      dev.expl = "Dev. Exp.",
      `Moran I statistic (res)` = "Moran's I",
      `p-value` = "p-value"
    )
    # cols_label_with(fn=f_gsub)
  v_name <- gsub("Frag. total","fragtotal",v_name) %>%
    gsub("Frag. per se","fragperse",.) %>%
    gsub("Área per se","areaperse",.)
  gtsave(table,
         paste0(v_folder,v_name,".png"),
         vwidth = vw, vheight = vh)
}
# tabela de seleção dos modelos usados
df_tabsel <- read_csv("rds/tabsel_simples.csv")
f_tabsel_png <- \(dff){
  vpaths <- daply(dff,"contraste",\(dfi){
    vpath <- f_gt_table(dfi=select(dfi,-contraste),
                        v_name = dfi$contraste[1],
                        vw = 800)
    return(vpath)
  })
}
v_log <- f_tabsel_png(dff = df_tabsel)
if(all(v_log)) print("figura criada: figuras/tabsel_simples.png")
#
# predição a posteriori: fixo e fixo + aleatório
l_df <- readRDS(paste0(v_path,"rds/l_dfpred_simples.rds"))
# nefeito <- names(l_df)[[1]]
f_plotPI <- \(nefeito){
  # objeto para o gráfico
  v_range_x <- sapply(l_df,\(li){
    range(li[["apenas fixo"]]$Uefeito)
  }) %>% range()
  v_range_y <- sapply(l_df,\(li){
    range(select(li[["apenas fixo"]],starts_with("Q_")))
  }) %>% range()
  f_geom_legend <- list(
    # guide name
    annotate("text", x = 0.01, y =-5.5, 
             label = "Posterior Prediction Interval", hjust = 0) ,
    # mediana
    annotate("segment", x = 0.01, xend = 0.06, y = -6, yend = -6, 
             color = "black", linewidth = 1),
    annotate("text", x = 0.08, y =-6, 
             label = "median", hjust = 0) ,
    # quantile range
    annotate("rect", xmin = 0.01, xmax = 0.06, ymin = -6.75, ymax = -6.25,
             fill = "gray", alpha = 0.7),
    annotate("text", x = 0.08, y = -6.5,
             label = "90% quant. range", hjust = 0),
    # Site
    annotate("rect", xmin = 0.01, xmax = 0.06, ymin = -7.5, ymax = -7,
             fill = "#986868", alpha = 0.7),
    annotate("text", x = 0.08, y = -7.25,
             label = "s(log(U/U)) + s(log(U/U))|Site", hjust = 0),
    # overall
    annotate("rect", xmin = 0.01, xmax = 0.06, ymin = -8.25, ymax = -7.75,
             fill = "darkgreen", alpha = 0.7),
    annotate("text", x = 0.08, y = -8,
             label = "s(log(U/U)) + 0|Site", hjust = 0)
  )
  # padronização dos dados
  l_df[[nefeito]][["apenas fixo"]]$SiteCode <- "apenas fixo"
  ldfi <- lapply(l_df[[nefeito]],mutate,
                 contraste = nefeito,
                 SiteCode = factor(SiteCode)) %>% 
    lapply(.,select,-any_of("logOR"))
  # gráfico
  p <- ldfi[["fixo e aleat"]] %>% 
    ggplot(aes(x = Uefeito, group = SiteCode)) +
    # fixo e aleat
    geom_ribbon(aes(ymin = Q_0.05, ymax = Q_0.95, 
                    fill = "Quantile range", group = SiteCode), 
                alpha = 0.2) +
    geom_line(aes(y = Q_0.5, color = "Median", group = SiteCode), 
              alpha = 0.5) +
    scale_fill_manual(values = c("Quantile range" = "#986868")) +
    scale_color_manual(values = c("Median" = "#C04000")) +
    # apenas fixo
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale_fill() +
    geom_ribbon(data=ldfi[["apenas fixo"]],
                aes(ymin = Q_0.05, ymax = Q_0.95, 
                    fill = "Quantile range", group = SiteCode), 
                alpha = 0.3) +
    geom_line(data=ldfi[["apenas fixo"]],
              aes(y = Q_0.5, color = "Median", group = SiteCode), 
              alpha = 0.7) +
    scale_fill_manual(values = c("Quantile range" = "darkgreen")) +
    scale_color_manual(values = c("Median" = "black")) +
    # 0 x 0 
    geom_hline(yintercept = 0,color="darkgray",alpha=0.75) +
    geom_vline(xintercept = 0,color="darkgray",alpha=0.75) +
    # ajustes
    scale_x_continuous(expand = expansion(add = c(0,0))) +
    scale_y_continuous(expand = expansion(add = c(0,0))) +
    labs(x="log(U / U)",y="log Odds Ratio (goodness-of-fit: SAD sim. - obs.)",
         title=nefeito) +
    theme_classic() +
    theme(
      axis.title.x = element_text(hjust = 0.5, vjust = 0.5,margin = margin(t = -2.5)),
      axis.title.y = element_text(hjust = 0.5, vjust = 0.5,margin = margin(t = -30)),
      legend.position = "none",
      aspect.ratio = 1
    )
  if(nefeito=="Área per se"){
    p <- p + f_geom_legend
  }
  return(p)
}
l_p <- lapply(names(l_df),f_plotPI)
names(l_p) <- names(l_df)
p <- arrangeGrob(grobs=l_p,nrow=1)
ggsave("figuras/figura_final.png",p,
       width = 16,height = 7,
       units = "in", dpi = 300)
image_trim <- image_read("figuras/figura_final.png") %>% 
  image_trim()
image_write(image_trim,"figuras/figura_final.png","png")
################################
############ antigo ############
################################


# tabela de seleção das perguntas
vpaths <- list.files(path=paste0(v_path,"rds/"),
                     pattern=".csv",full.names = T)
lapply(vpaths,\(li){
  dfsel <- read_csv(li)
  d_ply(dfsel,"contraste",\(dfi){
    f_gt_table(dfi=select(dfi,-contraste),
               v_name = dfi$contraste[1],
               f_cols_label_with = "f_gsub",
               vw = 800)
  })
  vpath <- list.files(path = paste0(v_path,"tabelas"),
                      pattern = "table_reci",
                      full.names = TRUE)
  l_png <- lapply(vpath,image_read)
  tabela_final <- image_append(
    do.call("c",l_png),stack = TRUE
  )
  # v_name <- str_extract(li,"(?<=tabsel\\_)(.*?)(?=\\_)")
  image_write(tabela_final, 
              path = paste0(v_path,"figuras/tabsel_","simples",".png"), 
              format = "png")
  file.remove(vpath)  
})
# 
#
# figura final dos efeitos
df_avgpred




#######################################
############# antigo #############
#######################################
df_avgpred <- read_csv(paste0(v_path,"rds/df_avgpred.csv")) %>% 
  mutate(contraste = gsub("fragtotal","Frag. total",contraste) %>% 
           gsub("fragperse","Frag. per se",.))

df_tabsel_logOR_aud <- read_csv(file="dados/csv/df_tabelaSelecao_logOR_cpert.csv")
f_modcheio_table <- \(dff){
  d_ply(dff,"pair",f_gt_table)
  paths <- list.files(path = paste0(v_path,"tabelas"),
                      pattern = "table_reci",
                      full.names = TRUE)
  l_png <- lapply(paths,image_read)
  tabela_final <- image_append(
    do.call("c",l_png),stack = TRUE
  )
  image_write(tabela_final, 
              path = paste0(v_path,"figuras/tabsel_modeloscheios.png"), 
              format = "png")
  file.remove(paths)  
}
f_modcheio_table(df_tabsel_logOR_aud)
#
# Qual o melhor conjunto de covariáveis?
paths_lmd <- list.files("dados/Rdata","l_md_logOR_2a_cpert_",full.names = T)
names(paths_lmd) <- str_extract(paths_lmd,"(?<=cpert\\_)(.*?)(?=\\.rds)") %>% 
  gsub("contemp-ideal","fragtotal",.) %>% 
  gsub("contemp-non_frag","fragperse",.) %>% 
  gsub("non_frag-ideal","areaperse",.)
l_df <- lapply(paths_lmd,\(li){
  li <- readRDS(li)
  f_TabSelGAMM(li)
})
f_gsub <- \(vchar){
  gsub("rank","Rank",vchar) %>% 
    gsub("modelo","GAHM",.) %>% 
    gsub("df","est. coef.",.) %>% 
    gsub("dAICc","ΔAICc",.) %>% 
    gsub("weight" , "Weight (ΔAICc)",.) %>% 
    gsub("dev.expl" , "Dev. Exp.",.) %>% 
    gsub("Moran I statistic (res)", "Moran's I",.) %>% 
    gsub("p-value" , "p-value",.)
}
f_gt_table_final <- \(lpaths,v_stack=TRUE){
  l_image <- lapply(lpaths,image_read)
  file.remove(do.call("c",lpaths))
  image_append(do.call("c",l_image), stack = v_stack)
}
image_table <- lapply(names(l_df),\(v_){
  f_gt_table(dfi = l_df[[v_]] %>% mutate(rank = 1:n()),
             v_title = gsub("fragtotal","Frag. total",v_) %>% 
               gsub("fragperse","Frag. per se",.) %>%
               gsub("areaperse","Área per se",.),
             v_name = v_,
             f_cols_label_with = "f_gsub")
  
}) %>% f_gt_table_final
image_write(image_table,path=paste0(v_path,"figuras/tab_sel_cov_quant_efeito.png"))


# 
# distribuição do logOR em função da taxa U por classe de perturbação

# df_plot <- df_logOR %>%
#   select(contraste,SiteCode,itself,Uefeito,forest_succession) %>% 
#   rename(logOR = itself)
# df_plot %>% 
#   filter(forest_succession!="capoeira") %>% 
#   ggplot(aes(x=Uefeito,y=logOR,color=forest_succession)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap(~contraste,ncol=1)
#   

 
