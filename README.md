# Ptarmi-GUTT: Cecal Microbiomes in Ptarmigan
Long-term patterns in cecal microbial communities in Icelandic Rock Ptarmigan 
=====

This repository contains open-source code, data, & text files for the project exploring microbial communities in the cecal microbiome of Icelandic rock ptarmigan 

## Project Goals

* **Aim 1)** Characterize the cecal microbiome of Icelandic rock ptarmigan over a 12-year study period (2007-2018)

* **Aim 2)** Determine if the cecal microbiome predicts metrics of foraging behavior, diet digestibility, and body condition over time

### Repo Contents

* **Analyses:**
- *Ptarmi-GUTT_processing.r*: R script for analyzing microbiome data 
- *Ptarmi-GUTT_figures.R*: R script for generating figures 

* **Data:**
  1. Raw metadata:
    - Only_ceca_subsampling_2022.csv -- Ptarmigan metadata
    - PtarmiganHealthData_Dec_2023.csv -- Ptarmigan health data
    - 2018GizzardFFDM.csv -- Gizzard fat free dry mass from 2018
    - Diet_diversity.csv -- Crop contents (diversity of diet)
    - IcelandPtarmiganID_GoodPoorFoodPercent_BCI.csv -- Crop contents (% food type)

  2. Raw mothur output data:
    - GUTT.trim.contigs.good.unique.good.filter.good.unique.precluster.denovo.vsearch.pick.ASV.abund.asv.shared.gz - ASV table. Note: needs to be uncompressed
    - GUTT.trim.contigs.good.unique.good.filter.good.unique.precluster.denovo.vsearch.pick.ASV.abund.nr_v138_1.wang.taxonomy -- taxonomy file

  3. Processed data:
    - asv_rel.csv.gz* - Processed ASV relative abundance table. Note: needs to be uncompressed
    - diversity.csv - Microbiome alpha diversity metrics
    - metadata_diet.csv - Crop contents paired to metadata
    metadata_health - Health data paired to metadata
    - tax_top30.csv - Top 30 ASVs (Table 1)
    - asv_reps_rel.csv - Processed ASV relative abundance table for replicates (QC)
    - metadata_reps - Processed metadata for replicates (QC)



* **Figures:**
Fig 1A
Fig 1B
Fig 2A
Fig 2B
Fig 3A
Fig 3B
Fig 3C
Fig 4
Fig S1
Fig S2

* **Bin:**
- *DiversityFunctions.R*: R functions for calculating alpha diversity
- *MothurTools.R*: R tools to analyze data output from the mothur software pipeline

### Software Versions and Dependencies


## Contributors

[Amanda Stromecki](https://muscarellalab.github.io/people/): PhD Student, Department of Biology and Wildlife, University of Alaska Fairbanks

[Dr. Mario Muscarella](https://muscarellalab.github.io/people/):Principle Investigator, Assistant Professor, Institute of Arctic Biology, University of Alaska Fairbanks. Head of the [Muscarella Lab](https://muscarellalab.github.io/).
