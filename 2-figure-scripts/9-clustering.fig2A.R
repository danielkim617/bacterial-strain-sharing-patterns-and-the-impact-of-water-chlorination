source(file.path(here::here(),"0-config.R"))
#Import data
#R21_sample_list = readRDS(file.path(here::here(), "data/R21_sample_list_no_gps.rds"))
#Kraken.phy.genus.filter = readRDS(file = file.path(box.path, "Kenya_bacterial_strain_sharing/analysis/scripts/scripts_DK/Kraken.phy.genus.filter.rds"))
# 
# #DirichletMultinomial modeling (DMM) - Unsupervised clustering
# ######### Run on Savio #########################################################################################################################
# Kraken.phy.genus.filter <- readRDS(file = "/global/home/users/danielkim617/Kraken.phy.genus.filter.rds")
# 
# # Function to filter taxa by its prevalence
# retain_taxa_over = function(x, t = 0.20, rel_abundance_threshold = 0.0001) {
#   # Convert the OTU table to relative abundances (normalize within each sample)
#   relative_abundance <- otu_table(x) / colSums(otu_table(x))
#   
#   # Mark OTUs as present in a sample if relative abundance > rel_abundance_threshold
#   presence_absence <- relative_abundance > rel_abundance_threshold
#   
#   # Calculate prevalence as the proportion of samples (columns) where each taxon is present
#   taxa_prevalence <- rowSums(presence_absence) / ncol(otu_table(x))
#   
#   # Retain taxa with prevalence >= threshold t
#   taxa_to_keep <- taxa_prevalence >= t
#   
#   # Prune taxa to keep only the relevant ones
#   y <- prune_taxa(taxa_to_keep, x)
#   
#   return(y)
# }
# 
# Kraken.phy.genus.filter.sim <- retain_taxa_over(Kraken.phy.genus.filter, t = 0.30, rel_abundance_threshold = 0.0001)
# saveRDS(Kraken.phy.genus.filter.sim, file = file.path(box.path, "Kenya_bacterial_strain_sharing/analysis/scripts/scripts_DK/Kraken.phy.genus.filter.sim.rds"))
# 
# # To speed up, only consider the core taxa
# # that are prevalent at 0.1% relative abundance in 50% of the samples
# # (note that this is not strictly correct as information is being discarded; one alternative would be to aggregate rare taxa)
# 
# # Pick the OTU count matrix
# # and convert it into samples x taxa format
# dat <- abundances(Kraken.phy.genus.filter)
# count <- as.matrix(t(dat))
# #
# fit_genus_list = vector("list",5)
# #Use all taxa
# seed   =  617
# set.seed(seed); seeds=sample(1:1000, 5)
# for(i in 1:5) {
#   
#   set.seed(seeds[i])
#   
#   fit_genus <- mclapply(c(1:30), dmn, count=count, verbose=FALSE, mc.cores=24)
#   
#   fit_genus_list[[i]] = fit_genus
#   
#   print(i) 
#   
# }
# # Save list into a file
# saveRDS(fit_genus_list, file = "/global/home/users/danielkim617/fit_genus_list.rds")
# ##################################################################################################################################
# fit_genus_list_new = readRDS(file = "/Users/danielkim617/Library/CloudStorage/Box-Box/Kenya_bacterial_strain_sharing/analysis/scripts/scripts_DK/fit_genus_list_new.rds")
# 
# lplc1 <- base::sapply(fit_genus_list_new[[1]], DirichletMultinomial::laplace) # AIC / BIC / Laplace
# lplc2 <- base::sapply(fit_genus_list_new[[2]], DirichletMultinomial::laplace) # AIC / BIC / Laplace
# lplc3 <- base::sapply(fit_genus_list_new[[3]], DirichletMultinomial::laplace) # AIC / BIC / Laplace
# lplc4 <- base::sapply(fit_genus_list_new[[4]], DirichletMultinomial::laplace) # AIC / BIC / Laplace
# lplc5 <- base::sapply(fit_genus_list_new[[5]], DirichletMultinomial::laplace) # AIC / BIC / Laplace
# 
# lplc = data.frame(lplc1, lplc2, lplc3, lplc4, lplc5)
# lplc$k = 1:30
# lplc.m = melt(lplc, id = "k")
# 
# #Plot model fit
# lplodel.fit = ggplot(lplc.m, aes(x=k, y=value, col=variable)) +
#   geom_point() +
#   geom_line(alpha = 0.5) +
#   scale_x_continuous(breaks = seq(0,30,1)) +
#   labs(x="Number of Dirichlet Components", y="Model Fit (Laplace)", title="Dirichlet Multinomial Mixtures Model fit", col="Iteration") +
#   theme_bw() +
#   theme(text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16))
# 
# #Find the number of k for the best fit
# itr1 = fit_genus_list_new[[1]]
# itr2 = fit_genus_list_new[[2]]
# itr3 = fit_genus_list_new[[3]]
# itr4 = fit_genus_list_new[[4]]
# itr5 = fit_genus_list_new[[5]]
# #
# best1 <- itr1[[which.min(unlist(lplc1))]]
# best2 <- itr2[[which.min(unlist(lplc2))]]
# best3 <- itr3[[which.min(unlist(lplc3))]]
# best4 <- itr4[[which.min(unlist(lplc4))]]
# best5 <- itr5[[which.min(unlist(lplc5))]]
# 
# best = c(best1, best2, best3, best4, best5)
# 
# # k is determined as 6
# mixturewt(best1)
# 
# part <- apply(mixture(best1), 1, which.max)
# 
# for (k in seq(ncol(fitted(best1)))) {
#   d <- melt(fitted(best1))
#   colnames(d) <- c("OTU", "cluster", "value")
#   d <- subset(d, cluster == k) %>%
#     # Arrange OTUs by assignment strength
#     arrange(value) %>%
#     mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
#     # Only show the most important drivers
#     filter(abs(value) > quantile(abs(value), 0.8))     
#   
#   p <- ggplot(d, aes(x = OTU, y = value)) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     labs(title = paste("Top drivers: community type", k))
#   print(p)
# }
# 
# library(phateR)
# #data(tree.data)
# #plot(prcomp(tree.data$data)$x, col=tree.data$branches)
# #tree.phate <- phate(tree.data$data, gamma = 0)
# part1 = part %>% as.data.frame()
# colnames(part1) = "partition"
# part1$sample_id = rownames(part1)
# 
# Kraken.genus.phate <- phate(t(otu_table(Kraken.phy.genus.filter)), gamma = 0.5)
# 
# Kraken.genus.phate.tab = Kraken.genus.phate$embedding %>% as.data.frame()
# Kraken.genus.phate.tab$sample_id = rownames(Kraken.genus.phate.tab)
# Kraken.genus.phate.tab = merge(Kraken.genus.phate.tab, part1, by = "sample_id")
# Kraken.genus.phate.tab = merge(Kraken.genus.phate.tab, R21_sample_list, by = "sample_id")
# # tsne.bray = tax.bray.tsne.all$tsne
# # tsne.bray = merge(tsne.bray, part1, by = "sample_id")
# # tsne.bray = merge(tsne.bray, Kraken.phy.genus.alpha[,1:10], by = "sample_id")
# #tsne.bray = merge(tsne.bray, shortbred.alpha, by = "sample_id")
# 
# phate.partition = ggplot(Kraken.genus.phate.tab, aes(x=PHATE1, y=PHATE2)) +
#   geom_point(aes(color = as.factor(partition), shape = tr), size = 2, alpha = 0.6) +
#   labs(x="PHATE 1", y="PHATE 2", title="DMM", col="Partition") +
#   scale_shape_manual(values = c("Control" = 21, "Water" = 17)) +
#   scale_color_manual(values=c("#ffb703", "#023047", "#FB8500", "#219EBC", "#606C38", "#F94144")) + 
#   theme_bw() +
#   theme(text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16)) +
#   coord_fixed()
# 
# phate.age = ggplot(Kraken.genus.phate.tab, aes(x=PHATE1, y=PHATE2, color = as.factor(Age_group), fill = as.factor(Age_group), shape = tr)) +
#   geom_point(size = 2, alpha = 0.6) +
#   labs(x="PHATE 1", y="PHATE 2", title="DMM", col="Partition") +
#   scale_shape_manual(values = c("Control" = 21, "Water" = 1)) +
#   theme_bw() +
#   theme(text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16)) +
#   coord_fixed()
# 
# phate.shannon = ggplot(Kraken.genus.phate, aes(x=PHATE1, y=PHATE2, col=tsne.bray$Shannon)) +
#   geom_point() +
#   scale_color_viridis_c(option = "D", direction = -1) +
#   labs(x="PHATE 1", y="PHATE 2", title="DMM (Shannon)", col="Shannon") +
#   theme_bw() +
#   theme(text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16)) +
#   coord_fixed()
# 
# #ggplot(tree.phate, aes(x=PHATE1, y=PHATE2, col=tsne.bray$Shannon.y)) +
# #  geom_point() +
# #  scale_color_viridis_c(option = "D", direction = -1)
# 
# phate.area = ggplot(Kraken.genus.phate, aes(x=PHATE1, y=PHATE2, col=tsne.bray$area)) +
#   geom_point() +
#   labs(x="PHATE 1", y="PHATE 2", title="DMM (Area)", col="Area") +
#   theme_bw() +
#   theme(text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16)) +
#   coord_fixed()
# 
# 
# phate.treatment =  ggplot(Kraken.genus.phate, aes(x=PHATE1, y=PHATE2, col=tsne.bray$tr)) +
#   geom_point() +
#   labs(x="PHATE 1", y="PHATE 2", title="DMM (Treatment)", col="Treatment") +
#   theme_bw() +
#   theme(text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=16)) +
#   coord_fixed()
# 
# ggarrange(phate.partition, phate.age, phate.shannon, phate.area, phate.treatment, align = c("hv"))
# 

######## Clustering based on MASH distance
#Import data
R21_sample_list = readRDS(file.path(here::here(), "data/R21_sample_list_no_gps.rds"))
urban <- R21_sample_list %>% filter(area=="urban") %>% pull(sample_id) %>% setdiff(., "stools_12093101_SM-NA1P5")
rural <- R21_sample_list %>% filter(area=="rural") %>% pull(sample_id)
# MASH
mash.dist = parseDistanceDF(file.path(here::here(), "data/mash_report.dist"))
mash.dist[upper.tri(mash.dist)] = t(mash.dist)[upper.tri(mash.dist)]
mash.dist = as.matrix(mash.dist)
diag(mash.dist)[is.na(diag(mash.dist))] = 0
mash.dist = mash.dist[rownames(mash.dist) != "stools_12093101_SM-NA1P5", colnames(mash.dist) != "stools_12093101_SM-NA1P5"]
mash.dist.urban <- mash.dist[urban, urban]
mash.dist.rural <- mash.dist[rural, rural]

set.seed(333)
mash.umap = umap(mash.dist)
mash.umap.cluster = Mclust(mash.umap)
#urban
set.seed(123)
mash.umap.urban = umap(mash.dist.urban)
mash.umap.urban.cluster = Mclust(mash.umap.urban)
#rural
set.seed(617)
mash.umap.rural = umap(mash.dist.rural)
mash.umap.rural.cluster = Mclust(mash.umap.rural)

colors <- c(
  "#d9c87f",  # Cluster 1 (light yellow)
  "#a2312e",  # Cluster 8 (red)
  "#e38837",  # Cluster 2 (orange)
  "#629ec7",  # Cluster 3 (light blue)
  "#f4cb6f",  # Cluster 4 (golden yellow)
  "#6bb5d5",  # Cluster 5 (sky blue)
  "#507a5a",  # Cluster 6 (dark green)
  "#7b8434",  # Cluster 7 (olive green)
  "#5e5a4a"   # Cluster 9 (brownish-gray)
)

#Graph clusters based on dim reductions by UMAP,
mash.umap.plot.table = as.matrix(mash.umap.cluster$classification) %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  dplyr::rename(Cluster =V1) %>% mutate(Cluster = as.factor(Cluster)) %>%
  inner_join(as.matrix(mash.umap) %>% 
               as.data.frame() %>% 
               rownames_to_column("ID") %>% 
               dplyr::rename(X = V1, Y = V2)) %>% 
  inner_join(R21_sample_list %>% 
               dplyr::rename(ID = sample_id) %>% mutate(ID = as.character(ID)))
# urban
mash.umap.urban.plot.table = as.matrix(mash.umap.urban.cluster$classification) %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  dplyr::rename(Cluster =V1) %>% mutate(Cluster = as.factor(Cluster)) %>%
  inner_join(as.matrix(mash.umap.urban) %>% 
               as.data.frame() %>% 
               rownames_to_column("ID") %>% 
               dplyr::rename(X = V1, Y = V2)) %>% 
  inner_join(R21_sample_list %>% 
               dplyr::rename(ID = sample_id) %>% mutate(ID = as.character(ID)))
# rural
mash.umap.rural.plot.table = as.matrix(mash.umap.rural.cluster$classification) %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  dplyr::rename(Cluster =V1) %>% mutate(Cluster = as.factor(Cluster)) %>%
  inner_join(as.matrix(mash.umap.rural) %>% 
               as.data.frame() %>% 
               rownames_to_column("ID") %>% 
               dplyr::rename(X = V1, Y = V2)) %>% 
  inner_join(R21_sample_list %>% 
               dplyr::rename(ID = sample_id) %>% mutate(ID = as.character(ID)))

# Plot data
mash.plot = mash.umap.plot.table %>% 
  ggplot(aes(x = X, y = Y, group = Cluster)) + 
  geom_point(aes(color = tr, shape = Age_group), size = 3,alpha = 0.5) +
  #scale_color_manual(values = c("#F9C74F", "#90BE6D", "#577590","#F94144")) + 
  scale_color_manual(values = c("#F94144","#577590")) + 
  scale_shape_manual(values = c(16, 17, 15, 3)) +
  geom_mark_hull(concavity = 5, expand = 0.001, radius = 0) + 
  theme(
    text = element_text(size = 20),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) + 
  labs(x = "UMAP1", y = "UMAP2", col = "Age groups", shape = "treatment")
#urban
mash.urban.plot = mash.umap.urban.plot.table %>%
  ggplot(aes(x = X, y = Y, group = Cluster)) + 
  geom_point(aes(color = tr), size = 3, alpha = 0.5) +
  #scale_shape_manual(values = c(21, 17)) + 
  scale_color_manual(values = c("#F94144","#577590")) + 
  #scale_color_manual(values = colors) + 
  scale_shape_manual(values = c(16, 17, 15, 3)) +
  geom_mark_hull(concavity = 5, expand = 0.001, radius = 0) + 
  theme(
    text = element_text(size = 20),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) + 
  labs(x = "UMAP1", y = "UMAP2", col = "Cluster", shape = "treatment")
#rural
mash.rural.plot = mash.umap.rural.plot.table %>% 
  ggplot(aes(x = X, y = Y, group = Cluster)) + 
  geom_point(aes(color = tr), size = 3, alpha = 0.5) +
  scale_shape_manual(values = c(21, 17)) + 
  #scale_color_manual(values = colors) + 
  #scale_color_manual(values = c("#F9C74F", "#90BE6D", "#577590","#F94144")) + 
  scale_color_manual(values = c("#F94144","#577590")) + 
  #scale_shape_manual(values = c(17, 15)) +
  geom_mark_hull(concavity = 5, expand = 0.001, radius = 0) + 
  theme(
    text = element_text(size = 20),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) + 
  labs(x = "UMAP1", y = "UMAP2", col = "Cluster", shape = "treatment") 

#Save the figure into a file
# ggsave(file.path(here::here(),"data/figures/mash.plotp.pdf"), dpi = 300, scale = 0.3, width = 900, height = 600, units = "mm", plot = mash.plot)
ggsave(file.path(here::here(),"data/figures/mash.plotp.bystudy.pdf"), dpi = 300, scale = 0.3, width = 900, height = 300, units = "mm", plot = ggarrange(mash.rural.plot, mash.urban.plot))

# color by Age
#urban
mash.urban.plot.age = mash.umap.urban.plot.table %>%
  ggplot(aes(x = X, y = Y, group = Cluster)) + 
  geom_point(aes(color = Age_group), size = 3, alpha = 0.5) +
  #scale_shape_manual(values = c(21, 17)) + 
  #scale_color_manual(values = c("#F94144","#577590")) + 
  scale_color_manual(values = c("#F9C74F", "#90BE6D", "#577590","#F94144")) + 
  #scale_color_manual(values = colors) + 
  #scale_shape_manual(values = c(16, 17, 15, 3)) +
  geom_mark_hull(concavity = 5, expand = 0.001, radius = 0) + 
  theme(
    text = element_text(size = 20),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) + 
  labs(x = "UMAP1", y = "UMAP2", col = "Cluster", shape = "treatment")
#rural
mash.rural.plot.age = mash.umap.rural.plot.table %>% 
  ggplot(aes(x = X, y = Y, group = Cluster)) + 
  geom_point(aes(color = Age_group), size = 3, alpha = 0.5) +
  scale_shape_manual(values = c(21, 17)) + 
  #scale_color_manual(values = colors) + 
  scale_color_manual(values = c("#90BE6D", "#577590")) + 
  #scale_color_manual(values = c("#F94144","#577590")) + 
  #scale_shape_manual(values = c(17, 15)) +
  geom_mark_hull(concavity = 5, expand = 0.001, radius = 0) + 
  theme(
    text = element_text(size = 20),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) + 
  labs(x = "UMAP1", y = "UMAP2", col = "Cluster", shape = "treatment") 

ggarrange(mash.rural.plot.age, mash.urban.plot.age) # Figure 2A
# Save
ggsave(file.path(here::here(),"data/figures/mash.plot.bystudy.age.pdf"), dpi = 300, scale = 0.3, width = 1000, height = 300, units = "mm", plot = ggarrange(mash.rural.plot.age, mash.urban.plot.age))

# Urban
# Statistical testing
tr.table.urban = mash.umap.urban.plot.table %>% select(Cluster, tr) %>% table %>% as.data.frame()
tr.table.urban = xtabs(Freq ~ Cluster + tr, data = tr.table.urban)
# Chi-squared test
chisq.test(tr.table.urban)
# Age
age.table.urban = mash.umap.urban.plot.table %>% select(Cluster, Age_group) %>% table %>% as.data.frame()
age.table.urban = xtabs(Freq ~ Cluster + Age_group, data = age.table.urban)
# Chi-squared test
chisq.test(age.table.urban)
# Rural
# Statistical testing
tr.table.rural = mash.umap.rural.plot.table %>% select(Cluster, tr) %>% table %>% as.data.frame()
tr.table.rural = xtabs(Freq ~ Cluster + tr, data = tr.table.rural)
# Chi-squared test
chisq.test(tr.table.rural)
# Age
age.table.rural = mash.umap.rural.plot.table %>% select(Cluster, Age_group) %>% filter(Age_group %in% c("19 - 30 months", "> 30 months")) %>%
  droplevels() %>%   # Drops unused factor levels
  table() %>%
  as.data.frame()
age.table.rural = xtabs(Freq ~ Cluster + Age_group, data = age.table.rural)
# Chi-squared test
chisq.test(age.table.rural)





