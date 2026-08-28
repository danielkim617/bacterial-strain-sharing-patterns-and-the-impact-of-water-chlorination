source(file.path(here::here(),"0-config.R"))
#Load R object data
load(file.path(here::here(), "data/strain.freq.RData"))

# Western Kenya and Nairobi # of detected strains
# ecoli: 166 and 179
# camp: 9 and 13
# enterob: 21 and 24
# entero: 20 and 43
# kleb: 30 and 40
# bact: 46 and 53
# bifi: 41 and 75

# Figure S1, S2, S4
ecoli.freq.fig = make_split_heatmap(ecoli.freq, split_n = 4)
camp.mag.freq.fig = make_split_heatmap(camp.mag.freq, split_n = 1)
enterob.freq.fig = make_split_heatmap(enterob.freq %>% rowwise() %>% mutate(strain = gsub("Ente_", "Enterobacter_", strain)) %>% mutate(), split_n = 1)
enterco.freq.fig = make_split_heatmap(enterco.freq %>% rowwise() %>% mutate(strain = gsub("Ente_", "Enterococcus_", strain)), split_n = 1)
kleb.freq.fig = make_split_heatmap(kleb.freq, split_n = 2)
bact.freq.fig = make_split_heatmap(bact.freq, split_n = 2)
bifi.freq.fig = make_split_heatmap(bifi.freq, split_n = 2)

# Save into files
ggsave(file.path(here::here(),"data/figures/ecoli.freq.fig.pdf"), dpi = 300, scale = 0.3, width = 2000, height = 1000, units = "mm", plot = ecoli.freq.fig)
ggsave(file.path(here::here(),"data/figures/camp.mag.freq.fig.pdf"), dpi = 300, scale = 0.3, width = 500, height = 1000, units = "mm", plot = camp.mag.freq.fig)
ggsave(file.path(here::here(),"data/figures/enterob.freq.fig.pdf"), dpi = 300, scale = 0.3, width = 500, height = 1000, units = "mm", plot = enterob.freq.fig)
ggsave(file.path(here::here(),"data/figures/enterco.freq.fig.pdf"), dpi = 300, scale = 0.3, width = 500, height = 1000, units = "mm", plot = enterco.freq.fig)
ggsave(file.path(here::here(),"data/figures/kleb.freq.fig.pdf"), dpi = 300, scale = 0.3, width = 1000, height = 1000, units = "mm", plot = kleb.freq.fig)
ggsave(file.path(here::here(),"data/figures/bact.freq.fig.pdf"), dpi = 300, scale = 0.3, width = 1000, height = 1000, units = "mm", plot = bact.freq.fig)
ggsave(file.path(here::here(),"data/figures/bifi.freq.fig.pdf"), dpi = 300, scale = 0.3, width = 1000, height = 1000, units = "mm", plot = bifi.freq.fig)


