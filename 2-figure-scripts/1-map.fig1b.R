#Figure 1B
source(file.path(here::here(),"0-config.R"))

R21_sample_list = readRDS(file.path(here::here(), "data/R21_sample_list_no_gps.rds"))
#Get map with GPS coordinates marked on it
library(ggmap)
# register_stadiamaps([key]) <- enter API key

# Define bounding box for Kenya
bbox <- c(left = 33.8, bottom = -4.7, right = 41.9, top = 5.3)
# Get map from Stadia (or use source = "stamen")
# kenya_map1 <- get_stadiamap(bbox = bbox, zoom = 7, maptype = "alidade_smooth")
kenya_map1 <- get_stadiamap(bbox = bbox, zoom = 6, maptype = "outdoors")

# Plot the map with GPS points
entire_map1 <- ggmap(kenya_map1) +
  geom_point(data = R21_sample_list, aes(x = gps_hhlongitude, y = gps_hhlatitude, color = study), size = 2, alpha = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Longitude", y = "Latitude", title = "Study locations in Kenya", color = "Study") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("rect", xmin = 34.0, xmax = 34.898, ymin = -4.5, ymax = -4.4, fill = "black") +  # 100 km scale bar
  annotate("text", x = 34.45, y = -4.55, label = "100 km", size = 3)

coord1 = c(34.0, -4.5)
coord2 = c(34.898 + 0.0898, -4.4)
distance_m <- distHaversine(coord1, coord2)

#Save into files
ggsave(file.path(here::here(),"data/figures/entire_map1.pdf"), dpi = 300, scale = 0.3, width = 650, height = 450, units = "mm", plot = entire_map1)

######################################################################################################################################################
# map for rural area
kenya_map.rural <- get_stadiamap(bbox = c(left = 34.36421-0.05, bottom = 0.2663711-0.05, right = 34.85575+0.05, top = 0.5652816+0.05), zoom = 11, maptype = "outdoors")

ggmap(kenya_map.rural) + geom_point(data = R21_sample_list[R21_sample_list$area == "rural",], aes(x=gps_hhlongitude, y=gps_hhlatitude, col=location, shape = tr )) + labs(x="Longitude", y="Latitude") + theme_bw()
#Plot map for WASHB study
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

rural_map = ggmap(kenya_map.rural) +
  geom_point(data = R21_sample_list %>% filter(study == "WASHB"), aes(x = gps_hhlongitude, y = gps_hhlatitude, color = tr), size = 2, alpha = 0.5) +
  scale_color_manual(values = my_colors) + # Use manual color palette
  scale_shape_manual(values = c(16, 17, 18, 19)) + # Customize shapes if needed
  labs(x = "Longitude", y = "Latitude", title = "WASHB study area", subtitle = "",shape = "Treatment", color = "Block") +
  theme_bw() + # Clean black and white theme
  theme(legend.position = "right", # Adjust legend position
        plot.title = element_text(hjust = 0.5), # Center plot title
        plot.subtitle = element_text(hjust = 0.5), # Center plot subtitle
        text = element_text(size = 12), # Adjust text size
        axis.title = element_text(face = "bold"), # Bold axis titles
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  annotate("rect",
           xmin = 34.33, xmax = 34.33 + 0.0898,  # 10 km wide
           ymin = 0.23, ymax = 0.235,            # short height bar
           fill = "black") +
  annotate("text",
           x = 34.33 + 0.0898/2, y = 0.227,       # centered below the bar
           label = "10 km",
           size = 3)

coord1 = c(34.33, 0.23)
coord2 = c(34.33 + 0.0898, 0.235)
distance_m <- distHaversine(coord1, coord2)

#Save into files
ggsave(file.path(here::here(),"data/figures/rural_map.pdf"), dpi = 300, scale = 0.3, width = 650, height = 450, units = "mm", plot = rural_map)

# map for urban area
kenya_map.urban <- get_stadiamap(bbox = c(left = 36.73734-0.01, bottom = -1.316038-0.02, right = 36.79044+0.01, top = -1.298977+0.02), zoom = 13, maptype = "outdoors")

urban_map = ggmap(kenya_map.urban) +
  geom_point(data = R21_sample_list %>% filter(study == "Nairobi"), aes(x = gps_hhlongitude, y = gps_hhlatitude, color = tr), size = 2, alpha = 0.5) +
  scale_color_manual(values = my_colors) + # Use manual color palette
  scale_shape_manual(values = c(16, 17, 18, 19)) + # Customize shapes if needed
  labs(x = "Longitude", y = "Latitude", title = "Nairobi study area", subtitle = "", color = "Treatment") +
  theme_bw() + # Clean black and white theme
  theme(legend.position = "right", # Adjust legend position
        plot.title = element_text(hjust = 0.5), # Center plot title
        plot.subtitle = element_text(hjust = 0.5), # Center plot subtitle
        text = element_text(size = 12), # Adjust text size
        axis.title = element_text(face = "bold"), # Bold axis titles
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()) + # Remove minor grid lines
  annotate("rect",
           xmin = 36.73, xmax = 36.73 + 0.0090,  # ~1 km width
           ymin = -1.33, ymax = -1.3285,         # thin bar height
           fill = "black") +
  annotate("text",
           x = 36.7345, y = -1.331,              # centered below the bar
           label = "1 km",
           size = 3)
coord1 = c(36.73, -1.33)
coord2 = c(36.73 + 0.0090, -1.3285)
distance_m <- distHaversine(coord1, coord2)

#Save into files
ggsave(file.path(here::here(),"data/figures/urban_map.pdf"), dpi = 300, scale = 0.3, width = 650, height = 450, units = "mm", plot = urban_map)

