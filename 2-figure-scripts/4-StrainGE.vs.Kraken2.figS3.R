# Figure S3
source(file.path(here::here(),"0-config.R"))

Kraken.phy.genus.RA.tab = readRDS(file.path(here::here(), "data/Kraken.phy.genus.RA.tab.rds"))

#Plot figures
e.coli.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Escherichia"), aes(x=Abundance, y=strainge_RA)) +
  geom_point(alpha = 0.5, col = "red") +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Escherichia") +
  theme_minimal() +
  guides(col = "none")

camp.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Campylobacter"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Campylobacter") +
  theme_minimal() +
  guides(col = "none")

camp.mag.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Campylobacter_MAG"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Campylobacter w/ MAG") +
  theme_minimal() +
  guides(col = "none")

enterob.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Enterobacter"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Enterobacter") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  guides(col = "none")

enteroc.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Enterococcus"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Enterococcus") +
  #coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  guides(col = "none")

kleb.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Klebsiella"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Klebsiella") +
  #coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal() +
  guides(col = "none")

staph.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Staphylococcus"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Staphylococcus") +
  theme_minimal() +
  guides(col = "none")

pseudo.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Pseudomonas"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 0.5)) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Pseudomonas") +
  theme_minimal() +
  guides(col = "none")

acineto.comp.fig =  ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Acinetobacter"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Acinetobacter") +
  theme_minimal() +
  guides(col = "none")

sal.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Salmonella"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Salmonella") +
  theme_minimal() +
  guides(col = "none")

bac.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Bacteroides"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  #coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Bacteroides") +
  theme_minimal() +
  guides(col = "none")

bifi.comp.fig = ggplot(Kraken.phy.genus.RA.tab %>% filter(Genus == "Bifidobacterium"), aes(x=Abundance, y=strainge_RA, col=Genus)) +
  geom_point(alpha = 0.5) +
  #coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.3)) +
  labs(x="Relative abundance (%)\n by Kraken2", y="Relative abundance (%)\n by StrainGST", title = "Bifidobacterium") +
  theme_minimal() +
  guides(col = "none")

total.comp.fig = ggarrange(bac.comp.fig, bifi.comp.fig, e.coli.comp.fig, enteroc.comp.fig, kleb.comp.fig, enterob.comp.fig, camp.mag.comp.fig, pseudo.comp.fig, sal.comp.fig, staph.comp.fig, acineto.comp.fig)

#Save into a file
ggsave(file.path(here::here(),"data/figures/StrainGE.vs.Kraken2.pdf"), dpi = 300, scale = 0.3, width = 900, height = 650, units = "mm", plot = total.comp.fig)
