



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


load("./Rdata/2022_platy_gl.RData") # data_gl

source("./R_scripts/1b_2022heron_seq_load_filt.R")

intial <- filter_data(data_gl, filter_type = "basic")
data_gl_filtered <- intial$data_gl_filtered
data_genind <- intial$data_genind


null <- filter_plus_null(data_genind = data_genind)

data_gl_filtered <- null$data_gl_filtered
data_genind <- null$data_genind

load("./Rdata/2022_Heron_null_filt.RData") # data_genind
load("./Rdata/2022_data_gl_filtered_null.RData") # data_gl_filtered




data_genind_adult <- data_genind[grep("a", indNames(data_genind)), ] # subset by matching 'a' in individual names


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


