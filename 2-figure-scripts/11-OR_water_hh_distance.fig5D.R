# Figure 5D
source(file.path(here::here(),"0-config.R"))
#Load data
all.OR = readRDS(file = file.path(here::here(), "data/all.OR.rds"))

all.OR$area = factor(all.OR$area, levels = c("rural", "urban"), labels = c("WASH B study","Nairobi study"))

all.OR.fig = 
  ggplot(all.OR, aes(x = Odds_Ratio, y = var)) +
  geom_point(size = 3, color = "blue") +  # Add points for odds ratios
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, color = "black") +  # Add error bars
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +  # Add vertical reference line
  facet_grid(target~area) +
  theme_bw() +
  labs(
    title = "Odds Ratios of strain-sharing between households",
    x = "Odds Ratio",
    y = ""
  ) +
  #theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

#Save into files
ggsave(file.path(here::here(),"data/figures/all.OR.fig.pdf"), dpi = 300, scale = 0.3, width = 500, height = 400, units = "mm", plot = all.OR.fig)

