## Filtering function

#filtering process
# 1) genlight basic filtering and then genind
# 2) null filtering on genid and then genlight
# 3) further indiviudal subsetting on genind
# 4) clustered filtered on adults after PCA


# load libraries ----------------------------------------------------------
#install.packages('dartr')
library(dartR)
library(dartR.popgen)
library(PopGenReport)
library(adegenet)
library(tictoc)
library(HardyWeinberg)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggrepel)
library(hierfstat)
library(ape)
library(poppr)
library(pegas)
library(dbscan)
library(sp)
library(rgdal)
library(clustertend)
library(cluster)
library(plotly)
library(pheatmap)
library(dendextend)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek1") # set theme in code
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek3") # set theme in code


# import data -------------------------------------------------------------
## File description SNP_2: SNP 2 Rows Format: Each allele scored in a binary fashion ("1"=Presence and "0"=Absence). Heterozygotes are therefore 
#scored as 1/1 (presence for both alleles/both rows)
# This means first row is ref, second is variant. A 1 in ref means both alleles are ref, a '1' in variant means both are variants. 1/1 means one ref one variant. 
# data_gl <- gl.read.dart(filename = "./data/Report_DPlatyg23-7805_SNP_2 - Copy corrected.csv", ind.metafile = "./data/meta_platy_ordered.csv", topskip = 6)
# #data_gl <- gl.read.dart(filename = "./data/Report_DPlatyg23-7805_SNP_mapping_2 - Copy.csv", ind.metafile = "./data/meta_platy_ordered.csv", topskip = 6)
# data_gl <- gl.reassign.pop(data_gl, as.pop = "stage")
# #recalculate metrics
# data_gl <- gl.recalc.metrics(data_gl, v = 3) # recalculate loci metrics
# save(data_gl, file = file.path("./Rdata", "2022_platy_gl.RData"))
load("./Rdata/2022_platy_gl.RData") # data_gl

# choose 'basic' for basic filtering, suitable for: XXXXX. This should result in 63 individuals and ~786 loci.
# or
# choose 'null' for additional null allele filtering, suitable for heterozygosity analyses. This should result in XXX individuals and XXX loci.
source("./R_scripts/1b_2022heron_seq_load_filt.R")

# basic filtering ---------------------------------------------------------
intial <- filter_data(data_gl, filter_type = "basic")
data_gl_filtered <- intial$data_gl_filtered
data_genind <- intial$data_genind
# data_genind_adult <- intial$data_genind_adult
# data_gl_filtered_adult <- intial$data_gl_filtered_adult
# data_genind_progeny <-intial$data_genind_progeny
# data_gl_adult_unique <- intial$data_gl_adult_unique
# data_genind_adult_unique <- intial$data_genind_adult_unique


# null filtering ----------------------------------------------------------
null <- filter_plus_null(data_genind = data_genind)

data_gl_filtered <- null$data_gl_filtered
# data_genind <- null$data_genind
# data_genind_adult <- null$data_genind_adult
# #data_genind_adult_unique <- null$data_genind_adult_unique
# data_genind_progeny <- null$data_genind_progeny

## quick load of null allels filtering (may not be current) 
save(data_genind, file = file.path("./Rdata", "2022_Heron_null_filt.RData"))
load("./Rdata/2022_Heron_null_filt.RData") # data_genind
save(data_gl_filtered, file = file.path("./Rdata", "2022_data_gl_filtered_null.RData"))
load("./Rdata/2022_data_gl_filtered_null.RData") # data_gl_filtered


# various individual subsetting -------------------------------------------

#unique adults based on best callrate
genotype_matrix <- data_genind_adult@tab
(callrate <- rowMeans(!is.na(genotype_matrix)))
ind_names <- indNames(data_genind_adult)
genotypes <- data_genind_adult@other$ind.metrics$genotype  # Adjust if necessary
geno_df <- data.frame(individual = ind_names, genotype = genotypes, callrate = callrate, stringsAsFactors = FALSE)
best_geno_df <- geno_df %>% group_by(genotype) %>% slice_max(order_by = callrate, n = 1, with_ties = FALSE) %>% ungroup()
best_ind_names <- best_geno_df$individual
best_indices <- match(best_ind_names, indNames(data_genind_adult))
data_genind_adult_unique <- data_genind_adult[best_indices, ]
data_gl_adult_unique = gi2gl(data_genind_adult_unique, parallel = FALSE, verbose = NULL)


