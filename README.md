# README

## Overview
This project explores bacterial strain-sharing patterns between humans in urban and rural communities at low resource settings.

## Basic settings
You need to set filepaths based on your local environment to run the codes. You can do this by setting variables using `usethis::edit_r_environ()`.

## Description
Two folders, run in order. File names are numbered roughly in the order they should be run, and most carry a suffix pointing to the figure they correspond to (`fig2D`, `figS9`, etc). `data/` holds the raw inputs these scripts read from (Kraken2/Bracken taxonomy, StrainGE/StrainGST/StrainGR strain calls, sample metadata, reference lists) - the bioinformatics pipeline that produced them is not part of this repo.

1. `1-analysis` - R scripts that take the raw bioinformatics/metadata inputs in `data/` and turn them into the tables/RDS objects used downstream. `strain-sharing/` is the part of this that's specific to the strain-sharing pipeline.
2. `2-figure-scripts` - reads what `1-analysis` produced and makes the actual figures (PDFs).

### 1-analysis
- `1-alpha_diversity.fig2D.R` - builds the phyloseq object from Kraken/Bracken and gets alpha diversity (richness, Shannon) - Fig 2D
- `2-Prevalence_and_num_strains_StrainGE.fig2E.R` - strain prevalence and # of strains per sample from StrainGE - Fig 2E
- `3-StrainGE.vs.Kraken2.R` - checks StrainGE relative abundance against Kraken2
- `4-microbial_comm_DAA_p_top_abundance.R` - differential abundance + top taxa
- `5-HH_distance_comparison.FigS9.R` - household geographic distance comparisons - Fig S9

### 1-analysis/strain-sharing
- `1-straingst_pair.R` - finds sample pairs that share the same strain (StrainGST) - feeds everything below
- `2-strain_rel_frequency.R` - relative frequency of each strain
- `3-load_StrainGR4ML.R` - loads StrainGR output, tags species, preps the ML input
- `4-ML_hh_level_strain_sharing.R` - household-level sharing rates from the random forest calls (within/between HH, water vs control, rural vs urban)
- `5-age_group_strain-level-sharing.R` - proportion sharing per strain, by age group - feeds the heatmap
- `6-individual_strain_sharing_table.R` - individual-level sharing tables / network edges
- `7-BW_sharing_clusters_or_villages.R` - within vs. between cluster/village sharing rates
- `8-Comparison_w_permutated_strain_sharing_Savio.R` - permutes HH distances to test distance-decay (this one's meant to run on Savio, not locally)
- `9-strain_sharing_age_groups.R` - odds ratios for distance/treatment vs sharing, plus age-group sharing
- `10-individual_level_network_analysis.R` - writes network files for Gephi
- `11-Prevalence_p_sharing_within_HH.figS5.R` - probability of sharing within a household given someone's positive - Fig S5

### 2-figure-scripts
- `1-map.fig1b.R` - study site map - Fig 1B
- `2-alpha_diversity.fig2BC.R` - alpha diversity by age/treatment/area - Fig 2B-C
- `3-strain_frequency.figS1S2S4.R` - strain frequency by genus - Fig S1, S2, S4
- `4-StrainGE.vs.Kraken2.figS3.R` - StrainGE vs Kraken2 abundance plot - Fig S3
- `5-ML_ACNI_Gap.figS16.R` - ACNI gap distribution from the ML calls - Fig S16
- `6-HH_level_strain-sharing_rates.fig3AB5C.R` - within vs between HH sharing rates - Fig 3A, 3B, 5C
- `7-BW_studies_BW_HH_rates.figS12.R` - between-study vs between-HH rates - Fig S12
- `8-only_19_30_and_30_months.figS6.R` - same thing but restricted to the 19/30-month timepoints - Fig S6
- `9-clustering.fig2A.R` - DMM clustering of community composition - Fig 2A
- `10-Age_sharing_percentage_matrix.Fig3C5ESX.R` - age-pair sharing matrices - Fig 3C, 5E, S
- `11-OR_water_hh_distance.fig5D.R` - OR for sharing vs water treatment/HH distance - Fig 5D
- `12-DAA.figS13.R` - differential abundance results - Fig S13
- `13-shared_strain_heatmap.fig4BS7S8S10S11.R` - sharing heatmaps + trees - Fig 4B, S7, S8, S10, S11
