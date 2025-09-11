library(lme4)
df_md <- readRDS("0_apendices/EfeitoEscala/figuras/df_md.rds")
l_md <- list()
l_md[[1]] <- lmer(logit_Umed ~ log_S_obs_z * log_Ntotal_z * k_factor * lado_factor + (1|SiteCode),
                  data=df_md)
l_md[[2]] <- lmer(logit_Umed ~ log_S_obs_z * k_factor * lado_factor + (1|SiteCode),
                  data=df_md)
l_md[[3]] <- lmer(logit_Umed ~ S_obs_z * Ntotal_z * k_factor * lado_factor  + (1|SiteCode),
                  data=df_md)
l_md[[4]] <- lmer(logit_Umed ~ S_obs_z * k_factor * lado_factor + (1|SiteCode),
                  data=df_md)
saveRDS(l_md,"0_apendices/EfeitoEscala/figuras/l_md.rds")

  # mutate(across(c(k_factor,lado_factor),~as.numeric(as.character(.x))))

df_plot0 %>% 
  mutate(lado = as.numeric(as.character(lado_factor))) %>% 
  ggplot(aes(x=lado_factor,y=Umed)) +
  geom_point() +
  geom_line(group=1) +
  scale_y_continuous(labels = scales::scientific_format(digits = 2)) +
  facet_wrap(~k_factor,ncol=5,scales="free_y") +
  theme_classic()

