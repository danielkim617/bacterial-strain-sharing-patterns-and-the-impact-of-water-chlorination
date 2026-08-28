source(file.path(here::here(),"0-config.R"))
#Import data
sig.diff.all = read.csv(file.path(here::here(), "data/stats/sig.diff.all1.csv"))
Kraken.phy.genus.filter = readRDS(file = file.path(here::here(), "data/Kraken.phy.genus.filter.rds"))

tax_rank = tax_table(Kraken.phy.genus.filter)

sig.diff.all.rank = merge(sig.diff.all, tax_rank, by.x = "Taxonomy", by.y = "genus", all.x=T)

diff.genera.fig2 = 
  ggplot(sig.diff.all.rank, aes(x = Abundance_Estimate, y = forcats::fct_rev(Taxonomy))) +
  geom_point(size = 3, color = "blue") +  
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +  
  facet_grid(area+order~., scales = "free", space = "free") +
  theme_bw() +
  labs(
    title = "",
    x = "Treatment coefficient (95% CI)",
    y = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(file.path(here::here(),"data/figures/diff.genera.fig21.pdf"), dpi = 300, scale = 0.3, width = 400, height = 650, units = "mm", plot = diff.genera.fig2)
  