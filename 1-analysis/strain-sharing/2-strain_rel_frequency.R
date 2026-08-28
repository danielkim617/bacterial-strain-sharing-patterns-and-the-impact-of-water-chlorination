# this script is to calculate the relative frequency of strains identified from straingst results.
source(file.path(here::here(),"0-config.R"))
#Load R object
R21_sample_list = readRDS(file.path(here::here(), "data/Kraken.phy.genus.alpha.rds"))
#Loading straingst results
ecoli.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.ecoli.tsv"), sep = "\t")
camp.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter.tsv"), sep = "\t")
camp_all.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter_all.tsv"), sep = "\t")
camp_MAG.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Campylobacter_MAG.tsv"), sep = "\t")
enterob.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Enterobacter.tsv"), sep = "\t")
enterco.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Enterococcus.tsv") , sep = "\t")
kleb.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Klebsiella.tsv") , sep = "\t")
staph.strains = read.table(file.path(here::here(), "data/straingst_comb/comb_straingst.Staphylococcus.tsv") , sep = "\t")
#From Colin's DB
bacteroides.strains = read.table(file.path(here::here(), "data/colin_DB_straingst/bacteroides_straingst.tsv"), sep = "\t")
bifidobacterium.strains = read.table(file.path(here::here(),"data/colin_DB_straingst/bifidobacterium_straingst.tsv"), sep = "\t")

all = rbind(ecoli.strains, camp.strains, camp_all.strains, camp_MAG.strains, 
            enterob.strains, enterco.strains, kleb.strains, staph.strains, bacteroides.strains, bifidobacterium.strains)

saveRDS(all, file.path(here::here(), "data/straingst_comb/all_strainGST.tsv"))

#Run command to calculate strain frequency for each target organism and make figures
ecoli.freq = get_strain_freq(ecoli.strains, R21_sample_list)
camp.freq = get_strain_freq(camp.strains, R21_sample_list)
camp.mag.freq = get_strain_freq(camp_MAG.strains, R21_sample_list)
enterob.freq = get_strain_freq(enterob.strains, R21_sample_list)
enterco.freq = get_strain_freq(enterco.strains, R21_sample_list)
kleb.freq = get_strain_freq(kleb.strains, R21_sample_list)
staph.freq = get_strain_freq(staph.strains, R21_sample_list)
bact.freq = get_strain_freq(bacteroides.strains, R21_sample_list)
bifi.freq = get_strain_freq(bifidobacterium.strains, R21_sample_list)

# Save the objects to a file named "my_data.RData"
save(ecoli.freq, camp.freq, camp.mag.freq, enterob.freq, enterco.freq, kleb.freq, staph.freq, bact.freq, bifi.freq, file = file.path(here::here(), "data/strain.freq.RData"))
