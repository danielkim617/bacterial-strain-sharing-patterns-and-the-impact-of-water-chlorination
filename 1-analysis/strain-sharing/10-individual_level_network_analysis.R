# This scripts produce individual strain sharing network files for visualization in Gephi. 
source(file.path(here::here(),"0-config.R"))
#Load data
# Set of edge tables
load(file.path(here::here(), "data/all.ML.networks.RData"))
# Load sample data and network tables
R21_sample_list.rural = readRDS(file.path(here::here(), "data/R21_sample_list.rural.rds"))
R21_sample_list.urban = readRDS(file.path(here::here(), "data/R21_sample_list.urban.rds"))
# All sample pair combination and sharing info
urban.all.pairs = readRDS(file.path(here::here(), "data/urban.all.pairs.rds"))
rural.all.pairs = readRDS(file.path(here::here(), "data/rural.all.pairs.rds"))

# Creating graphs (Only connections between households)
comm_genera = c("Bact", "Bifi")
noncomm_genera = c("Esch", "Camp_MAG", "Kleb", "Entc", "Entb")

# Rural area
# All genera
g_rural_all <- graph_from_data_frame(rural_bw_All %>% filter(HH_type == "Between"), vertices = R21_sample_list.rural %>% select(sample_id, everything()), directed = FALSE)
# Commensals
g_rural_comm <- graph_from_data_frame(rural_bw_All %>% filter(HH_type == "Between" & genus %in% comm_genera), vertices = R21_sample_list.rural %>% select(sample_id, everything()), directed = FALSE)
# Non-Commensals
g_rural_noncomm <- graph_from_data_frame(rural_bw_All %>% filter(HH_type == "Between" & genus %in% noncomm_genera), vertices = R21_sample_list.rural %>% select(sample_id, everything()), directed = FALSE)

# Urban area
# All genera
g_urban_all <- graph_from_data_frame(urban_bw_All %>% filter(HH_type == "Between"), vertices = R21_sample_list.urban %>% select(sample_id, everything()), directed = FALSE)
# Commensals
g_urban_comm <- graph_from_data_frame(urban_bw_All %>% filter(HH_type == "Between" & genus %in% comm_genera), vertices = R21_sample_list.urban %>% select(sample_id, everything()), directed = FALSE)
# Non-Commensals
g_urban_noncomm <- graph_from_data_frame(urban_bw_All %>% filter(HH_type == "Between" & genus %in% noncomm_genera), vertices = R21_sample_list.urban %>% select(sample_id, everything()), directed = FALSE)

# Save graphs into files
# Rural
# All
write_graph(g_rural_all, file = file.path(here::here(), "data/networks/g_rural_all.graphml"), format = "graphml")
# Comm
write_graph(g_rural_comm, file = file.path(here::here(), "data/networks/g_rural_comm.graphml"), format = "graphml")
# Non-comm
write_graph(g_rural_noncomm, file = file.path(here::here(), "data/networks/g_rural_noncomm.graphml"), format = "graphml")

# Urban
# All
write_graph(g_urban_all, file = file.path(here::here(), "data/networks/g_urban_all.graphml"), format = "graphml")
# Comm
write_graph(g_urban_comm, file = file.path(here::here(), "data/networks/g_urban_comm.graphml"), format = "graphml")
# Non-comm
write_graph(g_urban_noncomm, file = file.path(here::here(), "data/networks/g_urban_noncomm.graphml"), format = "graphml")



