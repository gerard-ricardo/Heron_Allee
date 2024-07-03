
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
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code


# import data -------------------------------------------------------------

##  Dartseq platy------------------------------------------------------------

### Preliminary analysis DarT seq dataset

# packages
# install.packages('dartR')

# platy data
#meta_platy <- read.csv("./data/meta_platy.csv", head = T) # make sure samples are in same order as in data_gl
# meta_platy$id
# meta_platy_order <- read.csv("./data/reordered_platy.csv", head = T) # This is the order of samples in the DarT file (SNP mapping)
# meta_platy_order$id
# meta_platy_final <- left_join(meta_platy_order, meta_platy, by = "id") # joining and keeping left
# meta_platy_final$id


# platy
#data_gl <- gl.read.dart(filename = "./data/Report_DPlatyg23-7805_SNP_2 - Copy corrected.csv", ind.metafile = "./data/meta_platy_ordered.csv", topskip = 6)

# data_gl <- data_gl[data_gl@ind.names$stage != "egg", ]
# data_gl$other

#data_gl$n.loc
#data_gl$ind.names

# Assign population to dart genlight object (can be any variable that you want cluster information from)
#data_gl <- gl.reassign.pop(data_gl, as.pop = "stage")
#data_gl$pop
# C:/Users/gerar/OneDrive/1 Work/4 Writing/1 Allee effects/allee experiments
# Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Users/gerar/Desktop/plink_win64_20231018", sep = ";"))
# gl2vcf(data_gl, plink_path = 'C:/Users/gerar/Desktop/plink_win64_20231018', outfile = "platy_vcf", outpath = 'C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/allee_experiments/data')

# recalculate metrics
#data_gl <- gl.recalc.metrics(data_gl, v = 3) # recalculate loci metrics
#save(data_gl, file = file.path("./Rdata", "2022_platy_gl.RData"))
load("./Rdata/2022_platy_gl.RData")  #data_gl


# calculate coverage metrics - mean number of reads that cover reference (30 good). Inc depth/reads will beter this.
# summary(data_gl$other$loc.metrics$coverage)
data_gl$other$loc.metrics
data_gl$other$loc.metrics$coverage <- data_gl$other$loc.metrics$AvgCountRef + data_gl$other$loc.metrics$AvgCountSnp
mean(data_gl$other$loc.metrics$coverage) # 61.59835. PD = 15.6, AH = 35.8
min((data_gl$other$loc.metrics$coverage)) # 5.375. PD = 10, AH = 5
max((data_gl$other$loc.metrics$coverage)) # 974.6452. PD = 29, AH = 372
sd(data_gl$other$loc.metrics$coverage) / sqrt(1996) # 1.858626, AH = 0.66
# PD has relatively consistent coverage, but on the low-sde
# ACROS has decent coverage but high variable. Unsure how good at each species


# data filtering ----------------------------------------------------------
data_gl_filtered <- data_gl
#following the dartr suggested order (see tut5)
# gl <-gl.filter.secondaries(gl)
# gl <- gl.filter.rdepth(gl)
# gl <- gl.filter.reproducibility(gl)
# gl <- gl.filter.callrate(gl, method=”loc”)
# gl <- gl.filter.callrate(gl, method=”ind”)
# gl <- gl.filter.monomorphs(gl)

## pre-filtering
#platy = 79 ind, 29377 loc

##secondaries
gl.report.secondaries(data_gl_filtered)
data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method="random", verbose = 3) #remove loci fragment that shared SNPs. Only keep 1
#platy = 79 ind, 22952  loc

#rdepth
gl.report.rdepth(data_gl_filtered)
#platy = 5.9. Generally 10 is considered min
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered,  lower = 10, v = 3) # filter by loci callrate
#platy = 79 ind, 2170  loc

##reproducibility 
gl.report.reproducibility(data_gl_filtered )
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t=0.95,v=3) #filter out loci with limited reproducibility
#Platy at 95%: ind = 79, loci = 1991  

# callrate loci (non missing data)
gl.report.callrate(data_gl_filtered, method = "loc") 
#PLaty = 71% , acro 30%.   , 49% missing data, but ~57%, AH = 0.6 (40% missing data)
#so lost some reps, only coral to loose was pd4a, rest larvae or eggs.
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
##Platy at 70%: ind = 59, loci = 1002
#AH = 202  ind, 21844 loci


## call rate ind (non missing data). low could indicate poor extract or reference genome or contamination.
#individauls
gl.report.callrate(data_gl_filtered, method = "ind") 
# platy = 96%
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.7, v = 3) # filter by ind callrate
#Platy at 70%:  ind = 61, loci = 1002 
#AH = 202  ind, 50405 loci


data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3) # recalculate loci metrics


##others - not sure if needed
# data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v=3) #remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
# not sure if I need HWE filter because remove inbreeding
# data_gl_filtered <- gl.filter.hwe(data_gl_filtered, alpha_val = 0.05, subset = "each", multi_comp_method = 'bonferroni',v=3) #filter out loci that depart from H-W proportions
# list.match <- data_gl_filtered$loc.names[which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.05 & data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef > 0.05 & data_gl_filtered$other$loc.metrics$coverage > 5)] #remove loci based on minor allele frequency and low data coverage
# data_gl_filtered <- data_gl_filtered[,match(list.match, data_gl_filtered$loc.names)]#keep only loci in the list above



# population filtering and objects ----------------------------------------


# look into genotype as population
data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop = "genotype")
data_gl_filtered


# Convert GENIND OBJECT all indiv
data_genind <- gl2gi(data_gl_filtered)

# Filter out eggs and larvae to keep only adults
adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adults")
data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)

# Convert genind adults only
data_genind_adult <- gl2gi(data_gl_filtered_adult)



