# This script is to calculate strain sharing rates at the household level using the output from the random forest model.
source(file.path(here::here(),"0-config.R"))
#Load the output from ML strain sharing calls
ML_test=read.table(file.path(here::here(), "data/StrainGE_ML/ML_Esch_Bact_Camp_Bifi_Camp_MAG_Kleb_Entb_Entc.tsv"), header = T, sep = "\t")

R21_sample_list = readRDS(file.path(here::here(), "data/Kraken.phy.genus.alpha.rds"))
#Add Age information
ML_test$sample1_Age_group = R21_sample_list[match(ML_test$sample1, R21_sample_list$sample_id),"Age_group"]
ML_test$sample2_Age_group = R21_sample_list[match(ML_test$sample2, R21_sample_list$sample_id),"Age_group"]

#Remove one of the duplicate samples (removed stools_12093101_SM-NA1P5)
ML_test = ML_test %>% filter(sample1 != "stools_12093101_SM-NA1P5" & sample2 != "stools_12093101_SM-NA1P5")

# Retain pairs with strain-sharing call by a model
staingr.org.ML = ML_test %>% filter(combined_rf_call == 1)
write.table(staingr.org.ML, file.path(here::here(), "data/StrainGE_ML/staingr.org.ML.tsv"), sep = "\t", row.names = F, quote = F)

# Calculate HH level strain-sharing rates
hh.comb.ML.All = get_hh_rate("All", y=R21_sample_list, z=staingr.org.ML)
hh.comb.ML.comm = get_hh_rate("All", y=R21_sample_list, z=staingr.org.ML %>% filter(genus %in% c("Bifi", "Bact"))) %>% mutate(., genus = "comm")
hh.comb.ML.non_comm = get_hh_rate("All", y=R21_sample_list, z=staingr.org.ML %>% filter(genus %in% c("Esch", "Camp_MAG", "Entc", "Entb","Kleb"))) %>% mutate(., genus = "non_comm")

hh.comb.ML.Esch = get_hh_rate("Esch", y=R21_sample_list, z=staingr.org.ML)
hh.comb.ML.Bact = get_hh_rate("Bact", y=R21_sample_list, z=staingr.org.ML)
hh.comb.ML.Camp = get_hh_rate("Camp", y=R21_sample_list, z=staingr.org.ML)
hh.comb.ML.Camp.mag = get_hh_rate("Camp_MAG", y=R21_sample_list, z=staingr.org.ML)
hh.comb.ML.Bifi = get_hh_rate("Bifi", y=R21_sample_list, z=staingr.org.ML)
hh.comb.ML.Entc = get_hh_rate("Entc", y=R21_sample_list, z=staingr.org.ML)
hh.comb.ML.Entb = get_hh_rate("Entb", y=R21_sample_list, z=staingr.org.ML)
hh.comb.ML.Kleb = get_hh_rate("Kleb", y=R21_sample_list, z=staingr.org.ML)

save(hh.comb.ML.All, hh.comb.ML.comm, hh.comb.ML.non_comm, hh.comb.ML.Esch, hh.comb.ML.Bact,hh.comb.ML.Camp,hh.comb.ML.Camp.mag,hh.comb.ML.Bifi, hh.comb.ML.Entc, hh.comb.ML.Entb, hh.comb.ML.Kleb, file = file.path(here::here(), "data/HH.rates.RData"))
load(file.path(here::here(), "data/HH.rates.RData"))

############################################################################################################################
#Compare rates (Within vs. Between)
wit.bet.test.All = get.wit.bet.test(hh.comb.ML.All)
wit.bet.test.comm = get.wit.bet.test(hh.comb.ML.comm)
wit.bet.test.non_comm = get.wit.bet.test(hh.comb.ML.non_comm)
wit.bet.test.Esch = get.wit.bet.test(hh.comb.ML.Esch)
wit.bet.test.Bact = get.wit.bet.test(hh.comb.ML.Bact)
wit.bet.test.Camp = get.wit.bet.test(hh.comb.ML.Camp)
wit.bet.test.Camp.mag = get.wit.bet.test(hh.comb.ML.Camp.mag)
wit.bet.test.Bifi = get.wit.bet.test(hh.comb.ML.Bifi)
wit.bet.test.Entc = get.wit.bet.test(hh.comb.ML.Entc)
wit.bet.test.Entb = get.wit.bet.test(hh.comb.ML.Entb)
wit.bet.test.Kleb = get.wit.bet.test(hh.comb.ML.Kleb)

wit.bet.test.All$genus = "All"
wit.bet.test.comm$genus = "comm"
wit.bet.test.non_comm$genus = "non_comm"
wit.bet.test.Esch$genus = "Esch"
wit.bet.test.Bact$genus = "Bact"
wit.bet.test.Camp$genus = "Camp"
wit.bet.test.Camp.mag$genus = "Camp_MAG"
wit.bet.test.Bifi$genus = "Bifi"
wit.bet.test.Entc$genus = "Entc"
wit.bet.test.Entb$genus = "Entb"
wit.bet.test.Kleb$genus = "Kleb"

#stat test (Within vs. Between)
wit.bet.test = data.frame(rbind(wit.bet.test.All,wit.bet.test.comm, wit.bet.test.non_comm,wit.bet.test.Esch, wit.bet.test.Bact, wit.bet.test.Camp, wit.bet.test.Camp.mag, wit.bet.test.Bifi,wit.bet.test.Entc,wit.bet.test.Entb, wit.bet.test.Kleb))
#Save into a file
write.csv(wit.bet.test, file.path(here::here(), "data/stats/wit.bet.test.csv"), row.names = F, quote = F)

#Compare rates (Within vs. Between) with CI
wit.bet.test.All.CI = get.wit.bet.test.CI(hh.comb.ML.All)
wit.bet.test.comm.CI = get.wit.bet.test.CI(hh.comb.ML.comm)
wit.bet.test.non_comm.CI = get.wit.bet.test.CI(hh.comb.ML.non_comm)
wit.bet.test.Esch.CI = get.wit.bet.test.CI(hh.comb.ML.Esch)
wit.bet.test.Bact.CI = get.wit.bet.test.CI(hh.comb.ML.Bact)
wit.bet.test.Camp.CI = get.wit.bet.test.CI(hh.comb.ML.Camp)
wit.bet.test.Camp.mag.CI = get.wit.bet.test.CI(hh.comb.ML.Camp.mag)
wit.bet.test.Bifi.CI = get.wit.bet.test.CI(hh.comb.ML.Bifi)
wit.bet.test.Entc.CI = get.wit.bet.test.CI(hh.comb.ML.Entc)
wit.bet.test.Entb.CI = get.wit.bet.test.CI(hh.comb.ML.Entb)
wit.bet.test.Kleb.CI = get.wit.bet.test.CI(hh.comb.ML.Kleb)

wit.bet.test.All.CI$genus = "All"
wit.bet.test.comm.CI$genus = "comm"
wit.bet.test.non_comm.CI$genus = "non_comm"
wit.bet.test.Esch.CI$genus = "Esch"
wit.bet.test.Bact.CI$genus = "Bact"
wit.bet.test.Camp.CI$genus = "Camp"
wit.bet.test.Camp.mag.CI$genus = "Camp_MAG"
wit.bet.test.Bifi.CI$genus = "Bifi"
wit.bet.test.Entc.CI$genus = "Entc"
wit.bet.test.Entb.CI$genus = "Entb"
wit.bet.test.Kleb.CI$genus = "Kleb"

#stat test (Within vs. Between)
wit.bet.test.CI = data.frame(rbind(wit.bet.test.All.CI,wit.bet.test.comm.CI, wit.bet.test.non_comm.CI,wit.bet.test.Esch.CI, wit.bet.test.Bact.CI, wit.bet.test.Camp.CI, wit.bet.test.Camp.mag.CI, wit.bet.test.Bifi.CI,wit.bet.test.Entc.CI,wit.bet.test.Entb.CI, wit.bet.test.Kleb.CI))
#Save into a file
write.csv(wit.bet.test.CI, file.path(here::here(), "data/stats/wit.bet.test.CI.csv"), row.names = F, quote = F)


############################################################################################################################################
#Compare rates (Rural vs. Urban)
rural.urban.test.All = get.rural.urban.test(hh.comb.ML.All)
rural.urban.test.comm = get.rural.urban.test(hh.comb.ML.comm)
rural.urban.test.non_comm = get.rural.urban.test(hh.comb.ML.non_comm)
rural.urban.test.Esch = get.rural.urban.test(hh.comb.ML.Esch)
rural.urban.test.Bact = get.rural.urban.test(hh.comb.ML.Bact)
rural.urban.test.Camp = get.rural.urban.test(hh.comb.ML.Camp)
rural.urban.test.Camp.mag = get.rural.urban.test(hh.comb.ML.Camp.mag)
rural.urban.test.Bifi = get.rural.urban.test(hh.comb.ML.Bifi)
rural.urban.test.Entc = get.rural.urban.test(hh.comb.ML.Entc)
rural.urban.test.Entb = get.rural.urban.test(hh.comb.ML.Entb)
rural.urban.test.Kleb = get.rural.urban.test(hh.comb.ML.Kleb)

rural.urban.test.All$genus = "All"
rural.urban.test.comm$genus = "comm"
rural.urban.test.non_comm$genus = "non_comm"
rural.urban.test.Esch$genus = "Esch"
rural.urban.test.Bact$genus = "Bact"
rural.urban.test.Camp$genus = "Camp"
rural.urban.test.Camp.mag$genus = "Camp_MAG"
rural.urban.test.Bifi$genus = "Bifi"
rural.urban.test.Entc$genus = "Entc"
rural.urban.test.Entb$genus = "Entb"
rural.urban.test.Kleb$genus = "Kleb"

rural.urban.test = data.frame(rbind(rural.urban.test.All,rural.urban.test.comm, rural.urban.test.non_comm ,rural.urban.test.Esch, rural.urban.test.Bact, rural.urban.test.Camp, rural.urban.test.Camp.mag, rural.urban.test.Bifi, rural.urban.test.Entc, rural.urban.test.Entb, rural.urban.test.Kleb))
#Save into a file
write.csv(rural.urban.test, file.path(here::here(), "data/stats/rural.urban.test.csv"), row.names = F, quote = F)

rural.urban.test.CI.All = get.rural.urban.test.CI(hh.comb.ML.All)
rural.urban.test.CI.comm = get.rural.urban.test.CI(hh.comb.ML.comm)
rural.urban.test.CI.non_comm = get.rural.urban.test.CI(hh.comb.ML.non_comm)
rural.urban.test.CI.Esch = get.rural.urban.test.CI(hh.comb.ML.Esch)
rural.urban.test.CI.Bact = get.rural.urban.test.CI(hh.comb.ML.Bact)
rural.urban.test.CI.Camp = get.rural.urban.test.CI(hh.comb.ML.Camp)
rural.urban.test.CI.Camp.mag = get.rural.urban.test.CI(hh.comb.ML.Camp.mag)
rural.urban.test.CI.Bifi = get.rural.urban.test.CI(hh.comb.ML.Bifi)
rural.urban.test.CI.Entc = get.rural.urban.test.CI(hh.comb.ML.Entc)
rural.urban.test.CI.Entb = get.rural.urban.test.CI(hh.comb.ML.Entb)
rural.urban.test.CI.Kleb = get.rural.urban.test.CI(hh.comb.ML.Kleb)

rural.urban.test.CI.All$genus = "All"
rural.urban.test.CI.comm$genus = "comm"
rural.urban.test.CI.non_comm$genus = "non_comm"
rural.urban.test.CI.Esch$genus = "Esch"
rural.urban.test.CI.Bact$genus = "Bact"
rural.urban.test.CI.Camp$genus = "Camp"
rural.urban.test.CI.Camp.mag$genus = "Camp_MAG"
rural.urban.test.CI.Bifi$genus = "Bifi"
rural.urban.test.CI.Entc$genus = "Entc"
rural.urban.test.CI.Entb$genus = "Entb"
rural.urban.test.CI.Kleb$genus = "Kleb"

rural.urban.test.CI = data.frame(rbind(rural.urban.test.CI.All,rural.urban.test.CI.comm, rural.urban.test.CI.non_comm ,rural.urban.test.CI.Esch, rural.urban.test.CI.Bact, rural.urban.test.CI.Camp, rural.urban.test.CI.Camp.mag, rural.urban.test.CI.Bifi, rural.urban.test.CI.Entc, rural.urban.test.CI.Entb, rural.urban.test.CI.Kleb))
#Save into a file
write.csv(rural.urban.test.CI, file.path(here::here(), "data/stats/rural.urban.test.CI.csv"), row.names = F, quote = F)

#Compare rates (Water vs. Control)
wat.cont.test.All = get.wat.cont.test(hh.comb.ML.All)
wat.cont.test.comm = get.wat.cont.test(hh.comb.ML.comm)
wat.cont.test.non_comm = get.wat.cont.test(hh.comb.ML.non_comm)
wat.cont.test.Esch = get.wat.cont.test(hh.comb.ML.Esch)
wat.cont.test.Bact = get.wat.cont.test(hh.comb.ML.Bact)
wat.cont.test.Camp = get.wat.cont.test(hh.comb.ML.Camp)
wat.cont.test.Camp.mag = get.wat.cont.test(hh.comb.ML.Camp.mag)
wat.cont.test.Bifi = get.wat.cont.test(hh.comb.ML.Bifi)
wat.cont.test.Entc= get.wat.cont.test(hh.comb.ML.Entc)
wat.cont.test.Entb = get.wat.cont.test(hh.comb.ML.Entb)
wat.cont.test.Kleb = get.wat.cont.test(hh.comb.ML.Kleb)

wat.cont.test.All$genus = "All"
wat.cont.test.comm$genus = "comm"
wat.cont.test.non_comm$genus = "non_comm"
wat.cont.test.All$genus = "All"
wat.cont.test.Esch$genus = "Esch"
wat.cont.test.Bact$genus = "Bact"
wat.cont.test.Camp$genus = "Camp"
wat.cont.test.Camp.mag$genus = "Camp_MAG"
wat.cont.test.Bifi$genus = "Bifi"
wat.cont.test.Entc$genus = "Entc"
wat.cont.test.Entb$genus = "Entb"
wat.cont.test.Kleb$genus = "Kleb"

wat.cont.test = data.frame(rbind(wat.cont.test.All, wat.cont.test.comm, wat.cont.test.non_comm, wat.cont.test.Esch, wat.cont.test.Bact, wat.cont.test.Camp, wat.cont.test.Camp.mag, wat.cont.test.Bifi,wat.cont.test.Entc,wat.cont.test.Entb,wat.cont.test.Kleb))
#Save into a file
write.csv(wat.cont.test, file.path(here::here(), "data/stats/wat.cont.test.csv"), row.names = F, quote = F)

################################################################################################################################################
# For comparisons between studies
################################################################################################################################################

#Load the output from ML strain sharing calls
ML_test=read.table(file.path(here::here(), "data/StrainGE_ML/ML_Esch_Bact_Camp_Bifi_Camp_MAG_Kleb_Entb_Entc.tsv"), header = T, sep = "\t")
ML_test_1 =  ML_test %>% filter(area != "") # To remove redundancy with the new append data
ML_bw=read.table(file.path(here::here(), "data/StrainGE_ML/bw_ML_results.tsv"), header = T, sep = "\t") # comparisons between studies
ML_test_2 = rbind(ML_test_1, ML_bw) # Combine dataset

R21_sample_list = readRDS(file.path(here::here(), "data/Kraken.phy.genus.alpha.rds")) # load sample metadata
#Add Age information
ML_test_2$sample1_Age_group = R21_sample_list[match(ML_test_2$sample1, R21_sample_list$sample_id),"Age_group"]
ML_test_2$sample2_Age_group = R21_sample_list[match(ML_test_2$sample2, R21_sample_list$sample_id),"Age_group"]

#Remove one of the duplicate samples (removed stools_12093101_SM-NA1P5)
ML_test_2 = ML_test_2 %>% filter(sample1 != "stools_12093101_SM-NA1P5" & sample2 != "stools_12093101_SM-NA1P5")

# Retain pairs with strain-sharing call by a model
staingr.org.ML_all = ML_test_2 %>% filter(combined_rf_call == 1)
write.table(staingr.org.ML_all, file.path(here::here(), "data/StrainGE_ML/staingr.org.ML_all.tsv"), sep = "\t", row.names = F, quote = F)

# Calculate HH level strain-sharing rates (all - including pairs bw studies)
bw_hh.comb.ML.All = get_hh_rate("All", y=R21_sample_list, z=staingr.org.ML_all)
bw_hh.comb.ML.comm = get_hh_rate("All", y=R21_sample_list, z=staingr.org.ML_all %>% filter(genus %in% c("Bifi", "Bact"))) %>% mutate(., genus = "comm")
bw_hh.comb.ML.non_comm = get_hh_rate("All", y=R21_sample_list, z=staingr.org.ML_all %>% filter(genus %in% c("Esch", "Camp_MAG", "Entc", "Entb","Kleb"))) %>% mutate(., genus = "non_comm")

bw_hh.comb.ML.Esch = get_hh_rate("Esch", y=R21_sample_list, z=staingr.org.ML_all)
bw_hh.comb.ML.Bact = get_hh_rate("Bact", y=R21_sample_list, z=staingr.org.ML_all)
bw_hh.comb.ML.Camp = get_hh_rate("Camp", y=R21_sample_list, z=staingr.org.ML_all)
bw_hh.comb.ML.Camp.mag = get_hh_rate("Camp_MAG", y=R21_sample_list, z=staingr.org.ML_all)
bw_hh.comb.ML.Bifi = get_hh_rate("Bifi", y=R21_sample_list, z=staingr.org.ML_all)
bw_hh.comb.ML.Entc = get_hh_rate("Entc", y=R21_sample_list, z=staingr.org.ML_all)
bw_hh.comb.ML.Entb = get_hh_rate("Entb", y=R21_sample_list, z=staingr.org.ML_all)
bw_hh.comb.ML.Kleb = get_hh_rate("Kleb", y=R21_sample_list, z=staingr.org.ML_all)
# save into files
save(bw_hh.comb.ML.All, bw_hh.comb.ML.comm, bw_hh.comb.ML.non_comm, bw_hh.comb.ML.Esch, bw_hh.comb.ML.Bact, bw_hh.comb.ML.Camp, bw_hh.comb.ML.Camp.mag, bw_hh.comb.ML.Bifi, bw_hh.comb.ML.Entc, bw_hh.comb.ML.Entb, bw_hh.comb.ML.Kleb, file = file.path(here::here(), "data/bw_HH.rates.RData"))
load(file.path(here::here(), "data/bw_HH.rates.RData"))

# Between HH rate comparisons
bw.studies.test.All = get.bw.studies.test(bw_hh.comb.ML.All) %>% mutate(genus = "All")
bw.studies.test.comm = get.bw.studies.test(bw_hh.comb.ML.comm) %>% mutate(genus = "comm")
bw.studies.test.non_comm = get.bw.studies.test(bw_hh.comb.ML.non_comm) %>% mutate(genus = "non_comm")
#
bw.studies.test.Esch = get.bw.studies.test(bw_hh.comb.ML.Esch) %>% mutate(genus = "Esch")
bw.studies.test.Bact = get.bw.studies.test(bw_hh.comb.ML.Bact) %>% mutate(genus = "Bact")
bw.studies.test.Camp.mag = get.bw.studies.test(bw_hh.comb.ML.Camp.mag) %>% mutate(genus = "Camp_MAG")
bw.studies.test.Bifi = get.bw.studies.test(bw_hh.comb.ML.Bifi) %>% mutate(genus = "Bifi")
bw.studies.test.Entc = get.bw.studies.test(bw_hh.comb.ML.Entc) %>% mutate(genus = "Entc")
bw.studies.test.Entb = get.bw.studies.test(bw_hh.comb.ML.Entb) %>% mutate(genus = "Entb")
bw.studies.test.Kleb = get.bw.studies.test(bw_hh.comb.ML.Kleb) %>% mutate(genus = "Kleb")

bw.studies.test = data.frame(rbind(bw.studies.test.All, bw.studies.test.comm, bw.studies.test.non_comm, bw.studies.test.Esch, bw.studies.test.Bact,
                                   bw.studies.test.Camp.mag, bw.studies.test.Bifi, bw.studies.test.Entc, bw.studies.test.Entb, bw.studies.test.Kleb))
#Save into a file
write.csv(bw.studies.test, file.path(here::here(), "data/stats/bw.studies.test.csv"), row.names = F, quote = F)

get.bw.studies.test.CI.All = get.bw.studies.test.CI(bw_hh.comb.ML.All) %>% mutate(genus = "All")
get.bw.studies.test.CI.comm = get.bw.studies.test.CI(bw_hh.comb.ML.comm) %>% mutate(genus = "comm")
get.bw.studies.test.CI.non_comm = get.bw.studies.test.CI(bw_hh.comb.ML.non_comm) %>% mutate(genus = "non_comm")
#
get.bw.studies.test.CI.Esch = get.bw.studies.test.CI(bw_hh.comb.ML.Esch) %>% mutate(genus = "Esch")
get.bw.studies.test.CI.Bact = get.bw.studies.test.CI(bw_hh.comb.ML.Bact) %>% mutate(genus = "Bact")
get.bw.studies.test.CI.Camp.mag = get.bw.studies.test.CI(bw_hh.comb.ML.Camp.mag) %>% mutate(genus = "Camp_MAG")
get.bw.studies.test.CI.Bifi = get.bw.studies.test.CI(bw_hh.comb.ML.Bifi) %>% mutate(genus = "Bifi")
get.bw.studies.test.CI.Entc = get.bw.studies.test.CI(bw_hh.comb.ML.Entc) %>% mutate(genus = "Entc")
get.bw.studies.test.CI.Entb = get.bw.studies.test.CI(bw_hh.comb.ML.Entb) %>% mutate(genus = "Entb")
get.bw.studies.test.CI.Kleb = get.bw.studies.test.CI(bw_hh.comb.ML.Kleb) %>% mutate(genus = "Kleb")

bw.studies.test.CI = data.frame(rbind(get.bw.studies.test.CI.All, get.bw.studies.test.CI.comm, get.bw.studies.test.CI.non_comm, get.bw.studies.test.CI.Esch, get.bw.studies.test.CI.Bact,
                                   get.bw.studies.test.CI.Camp.mag, get.bw.studies.test.CI.Bifi, get.bw.studies.test.CI.Entc, get.bw.studies.test.CI.Entb, get.bw.studies.test.CI.Kleb))
#Save into a file
write.csv(bw.studies.test.CI, file.path(here::here(), "data/stats/bw.studies.test.CI.csv"), row.names = F, quote = F)


##########################################################################################################################################################
## All areas, restricted to age groups common to rural and urban (excludes 0-19 months & mothers, which exist only in Nairobi). #############################################################################
filtered_urban = R21_sample_list %>% filter(!(Age_group %in% c("0 - 19 months", "Adults > 15 years")))
filtered_staingr.org.ML = staingr.org.ML %>% filter(!(sample1_Age_group %in% c("0 - 19 months", "Adults > 15 years")) &
                                                      !(sample2_Age_group %in% c("0 - 19 months", "Adults > 15 years")))

hh.comb.ML.All.filter = get_hh_rate("All", y=filtered_urban, z=filtered_staingr.org.ML)
hh.comb.ML.comm.filter = get_hh_rate("All", y=filtered_urban, z=filtered_staingr.org.ML %>% filter(genus %in% c("Bifi", "Bact"))) %>% mutate(., genus = "comm")
hh.comb.ML.non_comm.filter = get_hh_rate("All", y=filtered_urban, z=filtered_staingr.org.ML %>% filter(genus %in% c("Esch", "Camp_MAG", "Entc", "Entb","Kleb"))) %>% mutate(., genus = "non_comm")

hh.comb.ML.Esch.filter = get_hh_rate("Esch", y=filtered_urban, z=filtered_staingr.org.ML)
hh.comb.ML.Bact.filter = get_hh_rate("Bact", y=filtered_urban, z=filtered_staingr.org.ML)
hh.comb.ML.Camp.filter = get_hh_rate("Camp", y=filtered_urban, z=filtered_staingr.org.ML)
hh.comb.ML.Camp.mag.filter = get_hh_rate("Camp_MAG", y=filtered_urban, z=filtered_staingr.org.ML)
hh.comb.ML.Bifi.filter = get_hh_rate("Bifi", y=filtered_urban, z=filtered_staingr.org.ML)
hh.comb.ML.Entc.filter = get_hh_rate("Entc", y=filtered_urban, z=filtered_staingr.org.ML)
hh.comb.ML.Entb.filter = get_hh_rate("Entb", y=filtered_urban, z=filtered_staingr.org.ML)
hh.comb.ML.Kleb.filter = get_hh_rate("Kleb", y=filtered_urban, z=filtered_staingr.org.ML)


#Compare rates (Within vs. Between)
wit.bet.test.All.filter = get.wit.bet.test(hh.comb.ML.All.filter)
wit.bet.test.comm.filter = get.wit.bet.test(hh.comb.ML.comm.filter)
wit.bet.test.non_comm.filter = get.wit.bet.test(hh.comb.ML.non_comm.filter)
wit.bet.test.Esch.filter = get.wit.bet.test(hh.comb.ML.Esch.filter)
wit.bet.test.Bact.filter = get.wit.bet.test(hh.comb.ML.Bact.filter)
wit.bet.test.Camp.filter = get.wit.bet.test(hh.comb.ML.Camp.filter)
wit.bet.test.Camp.mag.filter = get.wit.bet.test(hh.comb.ML.Camp.mag.filter)
wit.bet.test.Bifi.filter = get.wit.bet.test(hh.comb.ML.Bifi.filter)
wit.bet.test.Entc.filter = get.wit.bet.test(hh.comb.ML.Entc.filter)
wit.bet.test.Entb.filter = get.wit.bet.test(hh.comb.ML.Entb.filter)
wit.bet.test.Kleb.filter = get.wit.bet.test(hh.comb.ML.Kleb.filter)

wit.bet.test.All.filter$genus = "All"
wit.bet.test.comm.filter$genus = "comm"
wit.bet.test.non_comm.filter$genus = "non_comm"
wit.bet.test.Esch.filter$genus = "Esch"
wit.bet.test.Bact.filter$genus = "Bact"
wit.bet.test.Camp.filter$genus = "Camp"
wit.bet.test.Camp.mag.filter$genus = "Camp_MAG"
wit.bet.test.Bifi.filter$genus = "Bifi"
wit.bet.test.Entc.filter$genus = "Entc"
wit.bet.test.Entb.filter$genus = "Entb"
wit.bet.test.Kleb.filter$genus = "Kleb"

#stat test (Within vs. Between)
wit.bet.test.filter = data.frame(rbind(wit.bet.test.All.filter,wit.bet.test.comm.filter, wit.bet.test.non_comm.filter,wit.bet.test.Esch.filter, wit.bet.test.Bact.filter, wit.bet.test.Camp.filter, wit.bet.test.Camp.mag.filter, wit.bet.test.Bifi.filter,wit.bet.test.Entc.filter,wit.bet.test.Entb.filter, wit.bet.test.Kleb.filter))
#Save into a file
write.csv(wit.bet.test.filter, file.path(here::here(), "data/stats/wit.bet.test.filter.csv"), row.names = F, quote = F)

#Compare rates (Rural vs. Urban)
rural.urban.test.All.filter = get.rural.urban.test(hh.comb.ML.All.filter)
rural.urban.test.comm.filter = get.rural.urban.test(hh.comb.ML.comm.filter)
rural.urban.test.non_comm.filter = get.rural.urban.test(hh.comb.ML.non_comm.filter)
rural.urban.test.Esch.filter = get.rural.urban.test(hh.comb.ML.Esch.filter)
rural.urban.test.Bact.filter = get.rural.urban.test(hh.comb.ML.Bact.filter)
rural.urban.test.Camp.filter = get.rural.urban.test(hh.comb.ML.Camp.filter)
rural.urban.test.Camp.mag.filter = get.rural.urban.test(hh.comb.ML.Camp.mag.filter)
rural.urban.test.Bifi.filter = get.rural.urban.test(hh.comb.ML.Bifi.filter)
rural.urban.test.Entc.filter = get.rural.urban.test(hh.comb.ML.Entc.filter)
rural.urban.test.Entb.filter = get.rural.urban.test(hh.comb.ML.Entb.filter)
rural.urban.test.Kleb.filter = get.rural.urban.test(hh.comb.ML.Kleb.filter)

rural.urban.test.All.filter$genus = "All"
rural.urban.test.comm.filter$genus = "comm"
rural.urban.test.non_comm.filter$genus = "non_comm"
rural.urban.test.Esch.filter$genus = "Esch"
rural.urban.test.Bact.filter$genus = "Bact"
rural.urban.test.Camp.filter$genus = "Camp"
rural.urban.test.Camp.mag.filter$genus = "Camp_MAG"
rural.urban.test.Bifi.filter$genus = "Bifi"
rural.urban.test.Entc.filter$genus = "Entc"
rural.urban.test.Entb.filter$genus = "Entb"
rural.urban.test.Kleb.filter$genus = "Kleb"

rural.urban.test.filter = data.frame(rbind(rural.urban.test.All.filter,rural.urban.test.comm.filter, rural.urban.test.non_comm.filter ,rural.urban.test.Esch.filter, rural.urban.test.Bact.filter, rural.urban.test.Camp.filter, rural.urban.test.Camp.mag.filter, rural.urban.test.Bifi.filter, rural.urban.test.Entc.filter, rural.urban.test.Entb.filter, rural.urban.test.Kleb.filter))
#Save into a file
write.csv(rural.urban.test.filter, file.path(here::here(), "data/stats/rural.urban.test.filter.csv"), row.names = F, quote = F)
################################################################################################################################################################################
####################################################################################################################################

#Compare rates (Within vs. Between) with CI
wit.bet.test.All.CI.filter = get.wit.bet.test.CI(hh.comb.ML.All.filter)
wit.bet.test.comm.CI.filter = get.wit.bet.test.CI(hh.comb.ML.comm.filter)
wit.bet.test.non_comm.CI.filter = get.wit.bet.test.CI(hh.comb.ML.non_comm.filter)
wit.bet.test.Esch.CI.filter = get.wit.bet.test.CI(hh.comb.ML.Esch.filter)
wit.bet.test.Bact.CI.filter = get.wit.bet.test.CI(hh.comb.ML.Bact.filter)
wit.bet.test.Camp.CI.filter = get.wit.bet.test.CI(hh.comb.ML.Camp.filter)
wit.bet.test.Camp.mag.CI.filter = get.wit.bet.test.CI(hh.comb.ML.Camp.mag.filter)
wit.bet.test.Bifi.CI.filter = get.wit.bet.test.CI(hh.comb.ML.Bifi.filter)
wit.bet.test.Entc.CI.filter = get.wit.bet.test.CI(hh.comb.ML.Entc.filter)
wit.bet.test.Entb.CI.filter = get.wit.bet.test.CI(hh.comb.ML.Entb.filter)
wit.bet.test.Kleb.CI.filter = get.wit.bet.test.CI(hh.comb.ML.Kleb.filter)

wit.bet.test.All.CI.filter$genus = "All"
wit.bet.test.comm.CI.filter$genus = "comm"
wit.bet.test.non_comm.CI.filter$genus = "non_comm"
wit.bet.test.Esch.CI.filter$genus = "Esch"
wit.bet.test.Bact.CI.filter$genus = "Bact"
wit.bet.test.Camp.CI.filter$genus = "Camp"
wit.bet.test.Camp.mag.CI.filter$genus = "Camp_MAG"
wit.bet.test.Bifi.CI.filter$genus = "Bifi"
wit.bet.test.Entc.CI.filter$genus = "Entc"
wit.bet.test.Entb.CI.filter$genus = "Entb"
wit.bet.test.Kleb.CI.filter$genus = "Kleb"

#stat test (Within vs. Between)
wit.bet.test.CI.filter = data.frame(rbind(wit.bet.test.All.CI.filter,wit.bet.test.comm.CI.filter, wit.bet.test.non_comm.CI.filter,wit.bet.test.Esch.CI.filter, wit.bet.test.Bact.CI.filter, wit.bet.test.Camp.CI.filter, wit.bet.test.Camp.mag.CI.filter, wit.bet.test.Bifi.CI.filter,wit.bet.test.Entc.CI.filter,wit.bet.test.Entb.CI.filter, wit.bet.test.Kleb.CI.filter))
#Save into a file
write.csv(wit.bet.test.CI.filter, file.path(here::here(), "data/stats/wit.bet.test.CI.filter.csv"), row.names = F, quote = F)
