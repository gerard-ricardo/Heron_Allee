#VCF

library(vcfR)

# Load VCF file
vcf <- read.vcfR("./data/platy_vcf.vcf")   #based on basic filtering for genlight_filtered
tidy_vcf <- vcfR2tidy(vcf)
vcf_fixed <- tidy_vcf$fix
vcf_genotypes <- tidy_vcf$gt
vcf_combined <- left_join(vcf_fixed, vcf_genotypes, by = c("ChromKey", "POS"))


parent_df <- vcf_combined %>% filter(grepl("\\.a\\.", Indiv))
progeny_df <- vcf_combined %>% filter(!grepl("\\.a\\.", Indiv))


pd1.a <- vcf_combined %>% filter(grepl("pd1\\.a\\.", Indiv))
