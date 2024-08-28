
##Todo
# add filtering function e.g basic, null etc




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



filter_data <- function(data, filter_type = "basic") {
  # Initialize the objects that may be returned
  data_gl_filtered <- data
  data_genind <- NULL
  data_genind_adult <- NULL
  data_gl_filtered_adult <- NULL
    


  
  #vcf
  # C:/Users/gerar/OneDrive/1 Work/4 Writing/1 Allee effects/allee experiments
  # Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Users/gerar/Desktop/plink_win64_20231018", sep = ";"))
  # gl2vcf(data_gl, plink_path = 'C:/Users/gerar/Desktop/plink_win64_20231018', outfile = "platy_vcf", outpath = 'C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/allee_experiments/data')
  
  
  # calculate coverage metrics - mean number of reads that cover reference (30 good). Inc depth/reads will beter this.
  # summary(data_gl$other$loc.metrics$coverage)
  data_gl$other$loc.metrics
  data_gl$other$loc.metrics$coverage <- data_gl$other$loc.metrics$AvgCountRef + data_gl$other$loc.metrics$AvgCountSnp
  median(data_gl$other$loc.metrics$coverage) #  PD = 11.54839
  min((data_gl$other$loc.metrics$coverage)) # PD = 5
  max((data_gl$other$loc.metrics$coverage)) #PD = 224.7
  sd(data_gl$other$loc.metrics$coverage) / sqrt(1996) # 0.134746
  # PD has relatively consistent coverage, but on the low-sde
  
  # data filtering ----------------------------------------------------------
  if (filter_type == "basic" || filter_type == "medium") {
    
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
  #platy = 79 ind, 2172    loc
  
  ##reproducibility 
  gl.report.reproducibility(data_gl_filtered )
  data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t=0.95, v=3) #filter out loci with limited reproducibility
  #Platy at 95%: ind = 79, loci = 1996    
  
  # callrate loci (non missing data)
  gl.report.callrate(data_gl_filtered, method = "loc") 
  #PLaty = 71%
  data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
  ##Platy at 70%: ind = 79, loci = 1004 
  #so lost some reps, only coral to loose was pd4a, rest larvae or eggs.
  
  #Minor Allele Frequency (MAF) and Coverage Filter:
  list.match <- data_gl_filtered$loc.names[
    which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.01 & 
            data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.99 & 
            data_gl_filtered$other$loc.metrics$OneRatioRef < 0.99 & 
            data_gl_filtered$other$loc.metrics$OneRatioRef > 0.01 & 
            data_gl_filtered$other$loc.metrics$coverage > 4)
  ]
  data_gl_filtered <- data_gl_filtered[, match(list.match, data_gl_filtered$loc.names)]
  #ind = 79, loci = 801  
  
  data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v=3) #remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
  #ind = 79, loci = 801  
  
  ## call rate ind (non missing data). low could indicate poor extract or reference genome or contamination.
  #individuals
  gl.report.callrate(data_gl_filtered, method = "ind") 
  #note that pd2.a.2, pd11.a.2, pd5.l.14.1. 'pd2.a.2' is the only pd2 so kinda important - however could add as unkown. 
  #'pd11.a.2' not important as we have rep 1. 'pd5.l.14.1' is the only pd5 larvae but assignment has been poor anyway, maybe incorrectly labelled pd5. Edit: Was incorrectl labelled in 13 from the tray sheet. 
  # platy = 89%
  pre_filt_ind <- data_gl_filtered@ind.names
  data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.58, v = 3) # filter by ind callrate
  filt_ind <- data_gl_filtered@ind.names
  (lost_ind <- setdiff(pre_filt_ind, filt_ind))
  #Platy at 58%:  ind = 63, loci = 786 : Used this theshold to allow for pd2.a.2 otherwise would not have a rep.
  #Note that if i filter at 85%, I only loose 3 individuals, so maybe more robust (can run with and without)
  length(filt_ind[grep('.l.', filt_ind )])  #count the larvae
  
  data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3) # recalculate loci metrics
  
  ##others - not sure if needed
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
  #genind object are 2-col (ref/var) loci format, where counts =  numbers of each allele i.e 2/0 means two reference. 
  
  # Filter out eggs and larvae to keep only adults
  adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adults")
  data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
  data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)
  
  #unique adults
  ind_names <- indNames(data_gl_filtered_adult)
  genotypes <- data_gl_filtered_adult@other$ind.metrics$genotype
  geno_df <- data.frame(individual = ind_names, genotype = genotypes, stringsAsFactors = FALSE)
  (unique_geno_df <- geno_df %>% distinct(genotype, .keep_all = TRUE))
  unique_ind_names <- unique_geno_df$individual
  unique_indices <- match(unique_ind_names, indNames(data_gl_filtered_adult))
  data_gl_adult_unique <- data_gl_filtered_adult[unique_indices, ]
  
  # Convert genind adults only
  data_genind_adult <- gl2gi(data_gl_filtered_adult)
  data_genind_adult_unique <- gl2gi(data_gl_adult_unique)
  
  #create 0_1 coded df
  mat_0_1_2_coded = data_genind_adult$tab
  mat_0_1_2_coded_char <- as.character(mat_0_1_2_coded)
  mat_0_1_2_coded_char[grepl("^2$", mat_0_1_2_coded_char)] <- "1"
  mat_0_1_coded <- matrix(as.numeric(mat_0_1_2_coded_char), nrow = nrow(mat_0_1_2_coded), ncol = ncol(mat_0_1_2_coded))
  
  }
  
  return(list(data_gl_filtered = data_gl_filtered, data_gl_filtered_adult = data_gl_filtered_adult, data_genind = data_genind, 
              data_genind_adult = data_genind_adult, data_genind_adult_unique = data_genind_adult_unique))
}


# filter likely null alleles (working)------------------------------------------------------------
filter_plus_null <- function(data_genind = data_genind, data_genind_adult = data_genind_adult, data_genind_adult_unique = data_genind_adult_unique, filter_type = "null") {
    ## all individual
  num_loci <- nLoc(data_genind) # Get the number of loci in the genind object
  sampled_loci_indices <- sample(num_loci, num_loci) # Randomly sample x loci (max popgenreport can report)
  
  # Subset the genind object to include only the sampled loci
  sampled_genind_obj <- data_genind[, sampled_loci_indices]
  pop(sampled_genind_obj) <- factor(rep("Combined_Population", nInd(sampled_genind_obj)))
  #table(pop(sampled_genind_obj))
  
  # Test for null alleles
  report1 = popgenreport(sampled_genind_obj, mk.null.all=TRUE, mk.pdf=FALSE)
  null_alleles_rep = report1$counts$nallelesbyloc
  null_alleles = colnames(null_alleles_rep)
  length(null_alleles)
  all_loci <- locNames(data_genind)
  length(all_loci)
  loci_to_keep <- setdiff(all_loci, null_alleles)
  data_genind <- data_genind[loc = loci_to_keep]
  
  
  ## all adults
  data_genind_adult
  num_loci <- nLoc(data_genind_adult) # Get the number of loci in the genind object
  sampled_loci_indices <- sample(num_loci, num_loci) # Randomly sample x loci (max popgenreport can report)
  
  # Subset the genind object to include only the sampled loci
  sampled_genind_obj <- data_genind_adult[, sampled_loci_indices]
  pop(sampled_genind_obj) <- factor(rep("Combined_Population", nInd(sampled_genind_obj)))
  #table(pop(sampled_genind_obj))
  
  # Test for null alleles
  report1 = popgenreport(sampled_genind_obj, mk.null.all=TRUE, mk.pdf=FALSE)
  null_alleles_rep = report1$counts$nallelesbyloc
  null_alleles = colnames(null_alleles_rep)
  length(null_alleles)
  all_loci <- locNames(data_genind_adult)
  length(all_loci)
  loci_to_keep <- setdiff(all_loci, null_alleles)
  data_genind_adult <- data_genind_adult[loc = loci_to_keep]
  
  ##  unique adults
  data_genind_adult_unique
  num_loci <- nLoc(data_genind_adult_unique) # Get the number of loci in the genind object
  sampled_loci_indices <- sample(num_loci, num_loci) # Randomly sample x loci (max popgenreport can report)
  
  # Subset the genind object to include only the sampled loci
  sampled_genind_obj <- data_genind_adult_unique[, sampled_loci_indices]
  pop(sampled_genind_obj) <- factor(rep("Combined_Population", nInd(sampled_genind_obj)))
  #table(pop(sampled_genind_obj))
  
  # Test for null alleles
  report1 = popgenreport(sampled_genind_obj, mk.null.all=TRUE, mk.pdf=FALSE)
  null_alleles_rep = report1$counts$nallelesbyloc
  null_alleles = colnames(null_alleles_rep)
  length(null_alleles)
  all_loci <- locNames(data_genind_adult_unique)
  length(all_loci)
  loci_to_keep <- setdiff(all_loci, null_alleles)
  data_genind_adult_unique <- data_genind_adult_unique[loc = loci_to_keep]

  return(list(data_genind = data_genind, data_genind_adult = data_genind_adult, data_genind_adult_unique = data_genind_adult_unique))

}



#save(data_genind_adult, file = file.path("./Rdata", "2022_Heron_null_filt_adult.RData"))



