f_fig_contraste1 <- \(df_plot){
  theme_set(theme_bw())
  l_p <- lapply(unique(df_plot$name),\(i){
    label <- case_when(grepl("area",i) ~ "Contraste Área per se: log(U_aglomerado / U_prístino)",
                    grepl("total",i) ~ "Contraste Frag. Total: log(U_Contemporâneo / U_prístino)",
                    grepl("perse",i) ~ "Contraste Frag. per se: log(U_Contemporâneo / U_aglomerado)")
    range_value <- range(df_plot$value)
    p1 <- df_plot %>% 
      filter(name==i) %>% 
      mutate(label = label) %>% 
      ggplot(aes(x=value)) +
      geom_histogram(bins = 120) +
      geom_boxplot(aes(y=-17.5),width=30) +
      labs(y="",x="") +
      scale_x_continuous(limits = range_value) +
      theme(plot.margin = unit(c(0.1,0.25,0.25,0), "cm")) +
      facet_wrap(~label)
    p2 <- df_plot %>% 
      filter(name==i) %>% 
      ggplot(aes(x=value)) +
      geom_boxplot() +
      theme_classic() +
      theme(plot.margin = unit(c(0,0,0,0), "cm")) +
      scale_x_continuous(limits = range_value) 
    #
    lims <- c(min(layer_scales(p1)$x$range$range, layer_scales(p2)$x$range$range),
              max(layer_scales(p1)$x$range$range, layer_scales(p2)$x$range$range))
    cowplot::plot_grid(p1 + ylim(lims),p2 + ylim(lims),ncol=1,rel_heights = c(16,3))
  })
  
}