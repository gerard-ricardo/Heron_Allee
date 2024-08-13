## Filtering function


# choose 'basic' for basic filtering, suitable for: XXXXX. This should result in 63 individuals and ~786 loci. 
# or
# choose 'null' for additional null allele filtering, suitable for heterozygosity analyses. This should result in XXX individuals and XXX loci. 
source("./R_scripts/1b_2022heron_seq_load_filt.R")

intial = filter_data(data_gl,  filter_type = 'basic')
data_gl_filtered = intial$data_gl_filtered
data_genind = intial$data_genind
data_genind_adult = intial$data_genind_adult

null =filter_plus_null(data_genind,  data_genind_adult, filter_type = 'null')
data_gl_filtered = null$data_gl_filtered
data_genind = null$data_genind
data_genind_adult = null$data_genind_adult




# quick load (may not be current) -----------------------------------------

load("./Rdata/2022_Heron_null_filt_adult.RData") #data_genind_adult
