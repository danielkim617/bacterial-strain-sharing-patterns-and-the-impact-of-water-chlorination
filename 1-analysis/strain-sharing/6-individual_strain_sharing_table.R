# This script produces individual level strain sharing tables
source(file.path(here::here(),"0-config.R"))
#Load data
staingr.org.ML = read.table(file.path(here::here(), "data/StrainGE_ML/staingr.org.ML.tsv"), sep = "\t", header = T)

#Get network tables
#Urban
urban_bw_All = staingr.org.ML %>% filter(area == "urban" & genus != "Camp")
#by genera
urban_bw_Bact = staingr.org.ML %>% filter(area == "urban" & genus == "Bact") 
urban_bw_Bifi = staingr.org.ML %>% filter(area == "urban" & genus == "Bifi") 
urban_bw_Esch = staingr.org.ML %>% filter(area == "urban" & genus == "Esch") 
urban_bw_Camp_MAG = staingr.org.ML %>% filter(area == "urban" & genus == "Camp_MAG") 
urban_bw_Entc = staingr.org.ML %>% filter(area == "urban" & genus == "Entc") 
urban_bw_Entb = staingr.org.ML %>% filter(area == "urban" & genus == "Entb") 
urban_bw_Kleb = staingr.org.ML %>% filter(area == "urban" & genus == "Kleb") 
#No duplicated connections between same sample pairs
urban_bw_All_no_dup = staingr.org.ML %>% filter(area == "urban"& genus != "Camp") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_comm_no_dup = staingr.org.ML %>% filter(area == "urban" & genus %in% c("Bact", "Bifi")) %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_noncomm_no_dup = staingr.org.ML %>% filter(area == "urban" & genus %in% c("Esch", "Camp_MAG", "Kleb", "Entc", "Entb")) %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_Bact_no_dup = staingr.org.ML %>% filter(area == "urban"& genus == "Bact") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_Bifi_no_dup = staingr.org.ML %>% filter(area == "urban"& genus == "Bifi") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_Esch_no_dup = staingr.org.ML %>% filter(area == "urban"& genus == "Esch") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_Camp_MAG_no_dup = staingr.org.ML %>% filter(area == "urban"& genus == "Camp_MAG") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_Entc_no_dup = staingr.org.ML %>% filter(area == "urban"& genus == "Entc") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_Entb_no_dup = staingr.org.ML %>% filter(area == "urban"& genus == "Entb") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
urban_bw_Kleb_no_dup = staingr.org.ML %>% filter(area == "urban"& genus == "Kleb") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique

#Rural
rural_bw_All = staingr.org.ML %>% filter(area == "rural" & genus != "Camp")
#by genera
rural_bw_Bact = staingr.org.ML %>% filter(area == "rural" & genus == "Bact") 
rural_bw_Bifi = staingr.org.ML %>% filter(area == "rural" & genus == "Bifi") 
rural_bw_Esch = staingr.org.ML %>% filter(area == "rural" & genus == "Esch") 
rural_bw_Camp_MAG = staingr.org.ML %>% filter(area == "rural" & genus == "Camp_MAG") 
rural_bw_Entc = staingr.org.ML %>% filter(area == "rural" & genus == "Entc") 
rural_bw_Entb = staingr.org.ML %>% filter(area == "rural" & genus == "Entb") 
rural_bw_Kleb = staingr.org.ML %>% filter(area == "rural" & genus == "Kleb") 
#No duplicated connections between same sample pairs
rural_bw_All_no_dup = staingr.org.ML %>% filter(area == "rural"& genus != "Camp") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_comm_no_dup = staingr.org.ML %>% filter(area == "rural" & genus %in% c("Bact", "Bifi")) %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_noncomm_no_dup = staingr.org.ML %>% filter(area == "rural" & genus %in% c("Esch", "Camp_MAG", "Kleb", "Entc", "Entb")) %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_Bact_no_dup = staingr.org.ML %>% filter(area == "rural"& genus == "Bact") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_Bifi_no_dup = staingr.org.ML %>% filter(area == "rural"& genus == "Bifi") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_Esch_no_dup = staingr.org.ML %>% filter(area == "rural"& genus == "Esch") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_Camp_MAG_no_dup = staingr.org.ML %>% filter(area == "rural"& genus == "Camp_MAG") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_Entc_no_dup = staingr.org.ML %>% filter(area == "rural"& genus == "Entc") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_Entb_no_dup = staingr.org.ML %>% filter(area == "rural"& genus == "Entb") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique
rural_bw_Kleb_no_dup = staingr.org.ML %>% filter(area == "rural"& genus == "Kleb") %>% select(sample1, sample2,sample1_study,sample1_area,sample1_hhid,sample1_tr,sample1_gps_lat,sample1_gps_lon,sample1_location,sample2_study,sample2_area,sample2_hhid,sample2_tr,sample2_gps_lat,sample2_gps_lon,sample2_location,distance,HH_type,tr,area,location,combined_rf_call,sample1_Age_group,sample2_Age_group) %>% unique

#Write network table into a file
write.csv(urban_bw_All, file.path(here::here(), "data/StrainGE_ML/networks/urban_bw_All.csv"), quote = F, row.names = F)
write.csv(rural_bw_All, file.path(here::here(), "data/StrainGE_ML/networks/rural_bw_All.csv"), quote = F, row.names = F)

#Write all data into a RData
save(urban_bw_All, urban_bw_All_no_dup, urban_bw_comm_no_dup, urban_bw_noncomm_no_dup, urban_bw_Bact, urban_bw_Bifi, 
     urban_bw_Esch, urban_bw_Camp_MAG, urban_bw_Entc, urban_bw_Entb, urban_bw_Kleb, rural_bw_All, rural_bw_All_no_dup, 
     rural_bw_comm_no_dup, rural_bw_noncomm_no_dup, rural_bw_Bact, rural_bw_Bifi, rural_bw_Esch, rural_bw_Camp_MAG, 
     rural_bw_Entc, rural_bw_Entb, rural_bw_Kleb, urban_bw_Bact_no_dup, urban_bw_Bifi_no_dup, urban_bw_Esch_no_dup, 
     urban_bw_Camp_MAG_no_dup, urban_bw_Entc_no_dup, urban_bw_Entb_no_dup, urban_bw_Kleb_no_dup, rural_bw_Bact_no_dup, rural_bw_Bifi_no_dup, rural_bw_Esch_no_dup, 
     rural_bw_Camp_MAG_no_dup, rural_bw_Entc_no_dup, rural_bw_Entb_no_dup, rural_bw_Kleb_no_dup, 
     file = file.path(here::here(), "data/all.ML.networks.RData"))

