#VCF

library(vcfR)


#vcf
##C:/Users/gerar/OneDrive/1 Work/4 Writing/1 Allee effects/allee experiments
#Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Users/gerar/Desktop/plink_win64_20231018", sep = ";"))
#gl2vcf(data_gl, plink_path = 'C:/Users/gerar/Desktop/plink_win64_20231018', outfile = "platy_vcf", outpath = 'C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/allee_experiments/data')
#Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Users/gerar/Desktop/plink_win64_20240818", sep = ";"))
#gl2vcf(data_gl_filtered, plink.bin.path = 'C:/Users/gerar/Desktop/plink_win64_20240818', outfile = "platy_vcf", outpath = './data')


# Load VCF file
vcf <- read.vcfR("./data/platy_vcf.vcf")   #based on basic filtering for genlight_filtered
tidy_vcf <- vcfR2tidy(vcf)
vcf_fixed <- tidy_vcf$fix
vcf_genotypes <- tidy_vcf$gt


vcf_combined <- left_join(vcf_fixed, vcf_genotypes, by = c("ChromKey", "POS"))


parent_df <- vcf_combined %>% filter(grepl("\\.a\\.", Indiv))
progeny_df <- vcf_combined %>% filter(!grepl("\\.a\\.", Indiv))


pd1.a <- vcf_combined %>% filter(grepl("pd1\\.a\\.", Indiv))
