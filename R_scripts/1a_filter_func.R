## Filtering function


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

intial <- filter_data(data_gl, filter_type = "basic")
data_gl_filtered <- intial$data_gl_filtered
data_genind <- intial$data_genind
data_genind_adult <- intial$data_genind_adult
data_gl_filtered_adult <- intial$data_gl_filtered_adult
data_gl_adult_unique <- intial$data_gl_adult_unique
data_genind_adult_unique <- intial$data_genind_adult_unique

null <- filter_plus_null(data_genind, data_genind_adult, data_genind_adult_unique, filter_type = "null")
data_gl_filtered <- null$data_gl_filtered
data_genind <- null$data_genind
data_genind_adult <- null$data_genind_adult
data_genind_adult_unique <- null$data_genind_adult_unique


# quick load (may not be current) -----------------------------------------

load("./Rdata/2022_Heron_null_filt_adult.RData") # data_genind_adult
