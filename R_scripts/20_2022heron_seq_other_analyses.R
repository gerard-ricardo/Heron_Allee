# 2023 Sequencing

# notes
# - call rate (NAs) for tenuis is low ~30%, for spath it is better, platy low too. Could be the reference genome or other issues
# - at4 likely spath
# - pca show similarity but not parental assignment.
# 3 larvae of tenuis attached well to adults using PCA (arne data), but some where way off (likely I used second batch).
# at1 and at2 from the leeward side were genetically distinct
# colonies 13 and 15 showed clear clustering to adults, but could be clones



# null alleles (not working)------------------------------------------------------------
# Get the number of loci in the genind object
num_loci <- nLoc(data_genind)

# Randomly sample 1000 loci (max popgenreport can report)
sampled_loci_indices <- sample(num_loci, 100)

# Subset the genind object to include only the sampled loci
sampled_genind_obj <- data_genind[, sampled_loci_indices]

# Test for null alleles
#null_allele_results <- null.all(data_genind)
popgenreport(sampled_genind_obj, mk.null.all=TRUE, mk.pdf=FALSE)





  







# loci under selection ----------------------------------------------------

outflank_results  = gl.outflank(data_genind_adult, plot = TRUE, RightTrimFraction = 0.12)
OutlierFlag = outflank_results$outflank$results$OutlierFlag
summary(OutlierFlag)
#no true outliers flagged







# pedigree and parentage --------------------------------------------------

gl.report.parent.offspring(data_gl_filtered)  #nothing picked up with this basic analysis


# identity by decent (IBD)------------------------------------------------------

res = gl.grm(data_gl_filtered) # relatedness matrix
res2 <- gl.grm.network(res, data_gl_filtered, relatedness_factor = 0.125)


# testing -----------------------------------------------------------------




# Acro Dartseq ------------------------------------------------------------

# meta_acro<-read.csv("./data/meta_acro.csv", head=T) #make sure samples are in same order as in data_gl
# meta_acro$id
# meta_acro_order<-read.csv("./data/reordered_acro.csv", head=T) #This is the order of samples in the DarT file (SNP mapping)
# meta_acro_order$id
# meta_acro_final = left_join(meta_acro_order, meta_acro, by = 'id')  #joining and keeping left
# meta_acro_final$id


# meta_acro_final<-meta_acro[match(meta_order$id,meta_acro_order$id),]
# write.csv(meta_platy_final,file="./data/meta_platy_ordered.csv")

## read Dart file with meta data TURN ONE OFF

data_gl <- gl.read.dart(filename = "./data/Report_DAc23-7804_SNP_2 - Copy.csv", ind.metafile = "./data/meta_acro_ordered2.csv", topskip = 6)


data_gl$n.loc
data_gl$ind.names

# gl2vcf(data_gl, plink_path = 'C:/Users/gerar/Desktop/plink_win64_20231018', outfile = "acro_vcf", outpath = 'C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/allee_experiments/data')


# For the Acro dataset, keep only the species of interest as analysing different species together could influence the analyses
data_gl <- gl.drop.pop(data_gl, "spath", as.pop = "species") # DROP population. tenuis or spath


# Assign population to dart genlight object (can be any variable that you want cluster information from)
data_gl <- gl.reassign.pop(data_gl, as.pop = "stage")


# data_gl$other$ind.metrics

# If you are subsetting a genlight object you need to adjust the metadata file and recalculate the genlight object metrics such as callrate,....
# meta_acro_spat_final<-meta_acro_final[match(data_gl$ind.names, meta_acro_final$id),]
# nrow(meta_acro_spat_final)# check how many individuals remain after removing tenuis samples

# recalculate metrics
data_gl <- gl.recalc.metrics(data_gl, v = 3) # recalculate loci metrics

data_gl_filtered <- data_gl


str(data_gl_filtered$other$loc.metrics)

# calculate Depth metrics - mean number of reads that cover reference (30 good). Inc depth/reads will beter this.
data_gl_filtered$other$loc.metrics$coverage <- data_gl_filtered$other$loc.metrics$AvgCountRef + data_gl_filtered$other$loc.metrics$AvgCountSnp

plot(density(data_gl_filtered$other$loc.metrics$coverage))
quantile(data_gl_filtered$other$loc.metrics$coverage)
data_gl_filtered$other$loc.metrics$coverage <- data_gl_filtered$other$loc.metrics$AvgCountRef + data_gl_filtered$other$loc.metrics$AvgCountSnp
mean(data_gl_filtered$other$loc.metrics$coverage) # 61.59835. PD = 15.6
median(data_gl_filtered$other$loc.metrics$coverage) # 28.60925. PD = 15.6

min((data_gl_filtered$other$loc.metrics$coverage)) # 5.375. PD = 10
max((data_gl_filtered$other$loc.metrics$coverage)) # 974.6452. PD = 29
sd(data_gl_filtered$other$loc.metrics$coverage) / sqrt(1996) # 1.858626
# PD has relatively consistent coverage, but on the low-sde
# ACROS has decent coverage but high variable. Unsure how good at each species   (28.60925 +1.858626)
# Calculate the proportion of coverage values less than 30
coverage_values <- data_gl_filtered$other$loc.metrics$coverage
less_than_30 <- coverage_values < 30
(proportion_less_than_30 <- sum(less_than_30) / length(coverage_values))



####

# call rate (missing data). Could indicate poor extract or reference genome or contamination.
# visualise data such as call rate for loci and individual. 95% considered good
plot(density(data_gl_filtered$other$loc.metrics$CallRate))
quantile(data_gl_filtered$other$loc.metrics$CallRate)
length(data_gl_filtered$other$loc.metrics$CallRate)
data_gl_filtered1 <- gl.report.callrate(data_gl_filtered, method = "loc") # acro 30%.   PLaty = 52%, 49% missing data, but ~57%
data_gl_filtered2 <- gl.report.callrate(data_gl_filtered, method = "ind") # so adults generally around 60% but variable, larvae about 60 but constant

# start with 1996 loci
# Filtering: this process is arbitrary and depends on the distribution and quality of your data. The bottom line is to end up with a dataset that has sufficient good quality input data for genetic comparison analyses
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.5, v = 3) # filter by ind callrate
# 29337
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate



data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
# 552
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.7, v = 3) # filter by ind callrate
# 552
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t = 0.7, v = 3) # filter out loci with limited reproducibility
# 552
data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v = 3) # remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
# 249
data_gl_filtered <- gl.filter.hwe(data_gl_filtered, alpha_val = 0.05, subset = "each", multi_comp_method = "bonferroni", v = 3) # filter out loci that depart from H-W proportions
# 161
data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method = "random", verbose = 3) # remove loci fragment that shared SNPs. Only keep 1
# 151
list.match <- data_gl_filtered$loc.names[which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.05 & data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef > 0.05 & data_gl_filtered$other$loc.metrics$coverage > 10)] # remove loci based on minor allele frequency and low data coverage
data_gl_filtered <- data_gl_filtered[, match(list.match, data_gl_filtered$loc.names)] # keep only loci in the list above
# 141 loci left


meta_acro_spat_final_filtered <- meta_acro_spat_final[match(data_gl_filtered$ind.names, meta_acro_spat_final$id), ] # match metadata file with genlight object

dim(meta_acro_spat_final_filtered)

## droplevels to avoid errors in dataframe size
meta_acro_spat_final_filtered$stage <- as.factor(as.character(meta_acro_spat_final_filtered$stage))
meta_acro_spat_final_filtered$genotype <- as.factor(as.character(meta_acro_spat_final_filtered$genotype))
meta_acro_spat_final_filtered$stage <- droplevels(meta_acro_spat_final_filtered$stage)
meta_acro_spat_final_filtered$genotype <- droplevels(meta_acro_spat_final_filtered$genotype)

data_gl_filtered$pop <- meta_allsamples$stage
data_gl_filtered$pop <- meta_allsamples$genotype

## PCA
# look into genotype as population
data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop = "genotype")
data_gl_filtered

# PCA params
pca_data <- tab(data_gl_filtered, freq = TRUE, NA.method = "mean")
pca <- dudi.pca(pca_data, center = T, scale = F, nf = 2, scannf = FALSE)

pca_complete <- data.frame(pca$li, pop = data_gl_filtered$pop)
# Assuming your data.frame is named pca_complete


pca_complete <- pca_complete %>%
  mutate(
    Stage = ifelse(str_detect(row.names(pca_complete), "\\.a\\."), "Adult", "Larva"),
    MumID = str_extract(row.names(pca_complete), "(?<=pd)\\d+"),
    RepID = str_extract(row.names(pca_complete), "(?<=\\.)\\d+$"),
    NewID = paste0(Stage, MumID, "_", RepID)
  )

data1 <- dplyr::arrange(pca_complete, Axis1) # dplyr - use this. Allows multiple sorts i.e site then orient

# plot pca
# spath
t <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = pop), shape = 21, size = 3) +
  geom_text(aes(label = rownames(pca_complete)), hjust = 0, vjust = 0) +
  scale_fill_manual(values = c("green", "red", "purple", "pink", "blue", "orange", "brown", "black")) # number of colors representing the number of groups in the variable
# scale_fill_manual(values=c("chartreuse2","green","darkgreen","tomato3","firebrick2","red","brown4","gray1","gray23","slategray1","steelblue1","cyan3","royalblue1","blue","purple","magenta4","yellow","ivory3","lavenderblush4","snow4"))
t

# tenuis
t1 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = pop), shape = 21, size = 3) +
  geom_text(aes(label = rownames(pca_complete)), hjust = 0, vjust = 0) +
  scale_fill_manual(values = c("green", "red", "purple", "pink", "blue", "orange", "brown", "steelblue", "gold", "yellow", "grey", "black")) # number of colors representing the number of groups in the variable
# scale_fill_manual(values=c("chartreuse2","green","darkgreen","tomato3","firebrick2","red","brown4","gray1","gray23","slategray1","steelblue1","cyan3","royalblue1","blue","purple","magenta4","yellow","ivory3","lavenderblush4","snow4"))
t1

# maybe I sequenced culture 2 (parents 14, 9, 10) instead of drop 2. Check remaining samples in freezer.

# platy
my_palette <- c(
  "dodgerblue", "salmon", "mediumseagreen", "orchid", "gold", "lightcoral",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato"
)




# Calculate the total variance
total_variance <- sum(pca$eig)
# Calculate the proportion of variance explained by the first principal component
var_explained_axis1 <- pca$eig[1] / total_variance
# Calculate the proportion of variance explained by the second principal component
var_explained_axis2 <- pca$eig[2] / total_variance
# Print the results
cat("Variance explained by Axis1: ", var_explained_axis1 * 100, "%\n")
cat("Variance explained by Axis2: ", var_explained_axis2 * 100, "%\n")


# Scree plot
# Calculate the proportion of variance explained by each component
pca$eig # first 4 >1
pve <- pca$eig / sum(pca$eig)

# Generate the scree plot
plot(pve,
  type = "b",
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained",
  main = "Scree Plot",
  ylim = c(0, max(pve)),
  yaxs = "i"
)

# Add a cumulative sum line
lines(cumsum(pve), type = "b", pch = 19, col = "red")
legend("topright",
  legend = c("PVE", "Cumulative PVE"),
  col = c("black", "red"), lty = 1:1, cex = 0.8
)




# admixture - not working#############
# Install and load the 'LEA' package
library(LEA)
library(adegenet)
# Convert dartR to genlight
data_genlight <- glDart(data_gl_filtered)
# Convert genlight to genind
X.snp <- df2genind(data_genlight, pop = data_genlight@pop)
# Run the admixture analysis
# Note: Set the K value according to your specific study
res.admixture <- snmf(X.snp, K = 3, entropy = TRUE, repetitions = 10)

# Plot the admixture results
barplot(res.admixture$Q)



# relatedness (not working) -------------------------------------------------------------
library(Relatedness)

# Prepare your genotype matrix and frequency matrix
# Assuming your genind object is called data_genind_filtered

# Convert genind object to a bi-allelic genotype matrix
genind_matrix <- as.matrix(data_genind_filtered@tab)

# Compute allele frequencies for each locus
allele_frequencies <- colMeans(genind_matrix, na.rm = TRUE) / 2

# Create a frequency matrix with the same number of columns as genind_matrix
# and one row representing the allele frequencies for each locus
frequency_matrix <- matrix(allele_frequencies, nrow = 1, ncol = length(allele_frequencies), byrow = TRUE)

# Ensure the dimensions match
if (ncol(genind_matrix) != ncol(frequency_matrix)) {
  stop("The number of columns in the genotype matrix and frequency matrix must match.")
}

# Run RelCoef function
relatedness_results <- RelCoef(IndividualGenom = genind_matrix, Freq = frequency_matrix)
# View relatedness results
print(relatedness_results$R)

####

## Usign distance matrix
dist_matrix <- bitwise.dist(data_gl_filtered)
dist_matrix_mat <- as.matrix(dist_matrix)
# Initialize a vector to store the closest individual
# Initialize vectors to store the closest individuals
closest_individual1 <- character(nrow(dist_matrix_mat))
closest_individual2 <- character(nrow(dist_matrix_mat))

# Loop through each individual
for (i in 1:nrow(dist_matrix_mat)) {
  # Exclude the diagonal element by setting it to a large value (Inf)
  dist_matrix_mat[i, i] <- Inf
  # Find the indices of the two minimum values in the row
  min_indices <- order(dist_matrix_mat[i, ])[1:2]
  # Get the names of the closest individuals
  closest_individual1[i] <- rownames(dist_matrix_mat)[min_indices[1]]
  closest_individual2[i] <- rownames(dist_matrix_mat)[min_indices[2]]
}

# Create a dataframe with the results
results <- data.frame(
  Individual = rownames(dist_matrix_mat),
  ClosestIndividual1 = closest_individual1,
  ClosestIndividual2 = closest_individual2
)

# Display the results
print(results)






########DAPC####


##### DAPC --> another clustering method
data_genind_filtered <- gi.reassign.pop(data_genind_filtered, as.pop = "genotype")

data_genind_filtered@pop <- meta_allsamples$stage
data_genind_alldata_neutral_050520$pop <- meta_neutral$pop
data_genind_alldata_neutral_050520$pop <- meta_neutral$Location

meta_allsamples$Location <- droplevels(meta_allsamples$Location)
met
clus <- find.clusters(data_gl_filtered, max.n.clust = 3)
100
3

dapc <- dapc(data_gl_filtered, k = 3)
100
6

scatter.dapc(dapc, col = c(1:3), clab = 0, solid = 1, scree.da = F, posi.leg = "bottomright", leg = F)

#  dartr base code----------------------------------------------------------------
# dartR
# snp = read.csv("./data/Report_DAc23-7804_SNP_2 - Copy.csv")
# meta = read.csv("./data/meta_acro.csv")

gl1 <- gl.read.dart(filename = "./data/Report_DAc23-7804_SNP_2 - Copy.csv", ind.metafile = "./data/meta_acro.csv") # note,
spp <- gl1@other$ind.metrics$species
# gl1_split <- split(gl1, spp)
# seppop(gl1, ~spp)
# seppop(gl1)
# names(gl1)


# file formats needs to be correct before import. Download as csv.
# gl.save(gl1,file="tmp.Rdata")
# gl1 <- gl.load("tmp.Rdata")
gl1@other$loc.metrics$RepAvg[1:10] ## Print the first 10 values of the RepAvg column in the loc.metrics object of gl1

gl1@ind.names
gl1@loc.names ## Print the names of the loci (SNPs) in gl1
gl1@gen[[1]] ## Print the genotype data for the first individual in gl1
# 45 genotypes,  1,996 SNPs
# missing data: 50605 (=56.34 %) scored as NA
nLoc(gl1) # number of loci = 1996

# filtering and diagnostics
library(HardyWeinberg)
# gl.diagnostics.hwe(gl1)
gl2 <- gl.filter.callrate(gl1, threshold = 0.8)
gl.diagnostics.hwe(
  x = gl.filter.allna(gl1[, 1:50]),
  stdErr = FALSE, n.cores = 1
)

gl1 <- gl.filter.monomorphs(gl1) # Filter out monomorphic sites (loci where all individuals have the same genotype)
# No change. 45 genotypes,  1,996 SNPs, missing data: 50605 (=56.34 %) scored as NA

# Identify individuals with missing genotypes
missing_inds <- which(is.na(gl1@gen))

gl.dist.pop(gl1[1:45, 1:100]) # distance matrix  , genetic distances between individuals.
# This matrix contains genetic distances between individuals, which can be used to construct a phylogenetic tree or perform other analyses
str(gl1)

pca1 <- gl.pcoa(gl2, verbose = 2)
gl.pcoa.plot(pca1, gl2) # pcoa  Fig. 1a. # The gl.pcoa.plot function creates a scatterplot of the first two principal coordinates,
# where each point represents an individual

gl.tree.nj(gl1) ## Construct a neighbor-joining tree from the genotype data in gl1




# cervus SNP_2 (not working) ------------------------------------------------------------

# load libraries --

# 1 Import data ---
# read.excel <- function(header=TRUE,...) {read.table("clipboard",sep="\t",header=header, na.strings=c("","-","na"),...)}
# data1 <- read.excel() #read clipboard from excel
# save(data1, file = file.path("./Rdata", "cleaned_acro1.RData"))
# load('./Rdata/cleaned_acro1.RData')
# head(data1)

# C:\Users\gerar\OneDrive\1 Work\3 Results\11 Allee effects\3 field experiments\2022_12 Heron\genetics\Report-DAc23-7804
# Read the data, skipping the first 6 rows and using the 7th row as header
data1 <- read.csv(
  skip = 6, header = FALSE,
  file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Report-DAc23-7804", "Report_DAc23-7804_SNP_2.csv")
)
head(data1)
colnames(data1) <- as.character(data1[1, ]) # Set the names of data1 to be the values from the sixth row
data1 <- data1[-1, ] # Remove the row that was used as header
data1 <- type.convert(data1, as.is = TRUE) # Convert data frame columns to the correct data type

# Check data structure to ensure types are correct
str(data1)

# Convert CloneID to a factor
data1$CloneID <- as.factor(as.character(data1$CloneID))

# Select relevant columns that start with 'sp' or 'at', along with CloneID, SNP, CallRate
data2 <- dplyr::select(data1, c(CloneID, SNP, CallRate), starts_with(c("sp", "at")))

# Group by CloneID and summarize to get the last SNP value and extract reference and variant alleles
data3 <- data2 %>%
  group_by(CloneID) %>%
  summarise(
    Last_SNP = last(SNP),
    Ref = substr(Last_SNP, start = nchar(Last_SNP) - 2, stop = nchar(Last_SNP) - 2),
    Var = substr(Last_SNP, start = nchar(Last_SNP), stop = nchar(Last_SNP))
  )

# Pivot the data to longer format, collapsing 'Ref' and 'Var' into a single column 'base'
data4 <- data3 %>% pivot_longer(cols = -c(CloneID, Last_SNP), names_to = "allele", values_to = "base")

# Add RowID to keep track of original row numbers
data4$RowID <- rep(c(1, 2), nrow(data4) / 2)
data2$RowID <- rep(c(1, 2), nrow(data2) / 2)

# Join data2 and data4 by CloneID and RowID
data5 <- left_join(data2, data4, by = c("CloneID", "RowID"))

# Create a new column 'LocusID' from 'CloneID' and select relevant columns
data6 <- data5 %>%
  mutate(LocusID = CloneID) %>%
  dplyr::select(., c(LocusID, allele, base), starts_with(c("sp", "at")))

# Transpose data6 to switch rows and columns, keeping only 'sp' and 'at' columns
data7 <- dplyr::select(data6, starts_with(c("sp", "at"))) %>% t()

# Group data6 by LocusID and create a modified LocusID with 'a' and 'b' appended
data8 <- data6 %>%
  group_by(LocusID) %>%
  mutate(LocusID_mod = ifelse(row_number() == 1, paste0(LocusID, "a"), paste0(LocusID, "b"))) %>%
  ungroup() %>%
  data.frame() %>%
  select(., LocusID_mod)

# Set the column names of data7 to be the modified LocusIDs with 'a' and 'b'
colnames(data7) <- data8$LocusID_mod
data7$

  # subset to work on managable version
  data7[1:10, 1:10]




#### scrap code?#############################

# data1.long = data1 %>% tidyr::pivot_longer(cols = 21:ncol(.),  names_to = "sampleID" ,values_to = "allele") %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)
# str(data1.long)
# nrow(data1.long)
# unique(data1.long$sampleID)
#
# #write.csv(data1.long, file = file.path("./data", "cleaned_acro1_long.csv"))
#
#
# ###
#
# #convert to id by loci
#
# library(dplyr)
#
# # Assuming your dataframe is named 'df'
# data2 <- data1 %>% select(1, 21:65)
# # Assuming your dataframe is named 'data2'
# # Transpose the dataframe
# # Saving the original column names (excluding CloneID)
# individual_ids <- colnames(data2)[-1]
# CloneID  = data2$CloneID
#
# # Transposing the dataframe excluding the CloneID column
# transposed_data <- t(data2[-1])
#
# # Converting the transposed matrix to a dataframe
# transposed_df <- as.data.frame(transposed_data)
#
#
# # Setting the saved column names as row names
# data3 = cbind(individual_ids, transposed_df)
# #head(data3)
# colnames(data3)[2:length(data3)] <- CloneID
# rownames(data3) <- NULL
#
# data3[data3 == 1] <- 2
# data3[data3 == 0] <- 1
#
# data3[is.na(data3)] <- '*'
#
# # Duplicate each column for each SNP (except the first column which contains individual IDs)
# snp_columns <- data3[, -1] # Excluding the first column which contains individual IDs
# duplicated_data <- lapply(snp_columns, function(col) data.frame(col1 = col, col2 = col))
# combined_data <- cbind(data3[1], do.call(cbind, duplicated_data))
# colnames(combined_data)[-1] <- paste0(rep(colnames(snp_columns), each = 2), rep(c('a', 'b'), length(colnames(snp_columns))))
#
# subset1 = combined_data[1:10]


# write.csv(data3, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "acro id by loci.csv"))


# Acro sequina ------------------------------------------------------------
library(sequoia)

# read.excel <- function(header=TRUE,...) {read.table("clipboard",sep="\t",header=header, na.strings=c("","-","na"),...)}
# data1.1 = read.excel() #read clipboard from excel. load ae7:bw3999
# save(data1.1, file = file.path("./Rdata", "2022acro_snp.RData"))
load("./Rdata/2022acro_snp.RData") # loaded individual as rows


data1.1[1:50, 1:42] # check subset to know it is okay - full df doesn't display right, or use f2
library(tidyr)
data1.1[is.na(data1.1)] <- -9 # turn na to -9
data1 <- data.frame(t(data1.1)) # rows ID, col = SNPs
data1$spp <- substr(rownames(data1), 1, 2) # extract values from character string
str(data1)
data1_spp <- split(data1, data1$spp)
data1_at <- data1_spp$at %>%
  dplyr::select(., -c(spp)) %>%
  as.matrix()
data1_sp <- data1_spp$sp %>%
  dplyr::select(., -c(spp)) %>%
  as.matrix()

# write.csv(data1_at,'./data/snp_ten.csv')

# mat1 = as.matrix(data1) #use f2 to view

# 1b Import meta data
# read.excel <- function(header=TRUE,...) {read.table("clipboard",sep="\t",header=header, na.strings=c("","-","na"),...)}
# data2=read.excel() #read clipboard from excel
# save(data2, file = file.path("./Rdata", "2022acro_meta.RData"))
load("./Rdata/2022acro_meta.RData")
data2
# write.csv(data2,'./Rdata/data2.csv')
data2$BirthYear <- 2022
data2 <- dplyr::select(data2, -c(age)) # remove column. Make sure have package on front
str(data2)
data2$ID <- as.factor(as.character(data2$ID))
data2$spp <- substr(data2$ID, 1, 2) # extract values from character string
data2_spp <- split(data2, data2$spp)
data2_at <- data2_spp$at
data2_sp <- data2_spp$sp



# checking
CheckGeno(data1_at, Plot = T)
# Warning:  There are 458 SNPs scored for <5% of individuals, these will be excluded
# Warning:  There are 988 monomorphic (fixed) SNPs, these will be excluded
# Warning:  In addition, there are 1149 SNPs scored for <50% of individuals
# After exclusion, There are  17  individuals and  2546  SNPs.

CheckGeno(data1_sp, Plot = T)
# Warning:  There are 878 SNPs scored for <5% of individuals, these will be excluded
# Warning:  There are 843 monomorphic (fixed) SNPs, these will be excluded
# Warning:  In addition, there are 704 SNPs scored for <50% of individuals
# After exclusion, There are  28  individuals and  2271  SNPs.


# args.AP1 <- list(type = "flat", ages = c(0, 1))


SeqOUT1_at <- sequoia(
  GenoM = data1_at,
  LifeHistData = data2_at,
  Err = 0.005, # genotyping error rate
  Module = "par",
  # args.AP = args.AP1,
  quiet = "verbose",
  Plot = T
)
# assigned 0 dams and 0 sires to 17 individuals


SeqOUT1_sp <- sequoia(
  GenoM = data1_sp,
  LifeHistData = data2_sp,
  Err = 0.005, # genotyping error rate
  Module = "par",
  # args.AP = args.AP1,
  quiet = "verbose",
  Plot = T
)
# assigned 0 dams and 0 sires to 28 individuals


# For each pair of candidate relatives, the likelihoods are calculated of them being parent-offspring (PO), full siblings (FS), half siblings (HS),
# grandparent-grandoffspring (GG), full avuncular (niece/nephew - aunt/uncle; FA), half avuncular/great-grandparental/cousins (HA), or unrelated (U).
# Assignments are made if the likelihood ratio (LLR) between the focal relationship and the most likely alternative exceed the threshold Tassign.

# Warning:  There are 2752 SNPs scored for <50% of individuals
# After exclusion, There are  45  individuals and  3991  SNPs.
# it's possible that the SNP is rare in the population you are studying, or it might be located in a region of the genome that is difficult to
# sequence accurately.#it could also be clones.

# Platys sequoia ------------------------------------------------------------------

library(sequoia)
#
# read.excel <- function(header=TRUE,...) {read.table("clipboard",sep="\t",header=header, na.strings=c("","-","na"),...)}
# data1.1 = read.excel() #read clipboard from excel. load ae7:bw3999
# save(data1.1, file = file.path("./Rdata", "2022platy_snp.RData"))
load("./Rdata/2022platy_snp.RData") # loaded individual as rows
data1.1[1:50, 1:42] # check subset to know it is okay - full df doesn't display right, or use f2
library(tidyr)
data1.1[is.na(data1.1)] <- -9 # turn na to -9
data1 <- t(data1.1)
data1[1:50, 1:42]
# mat1 = as.matrix(data1) #use f2 to view


# 1b Import meta data --
# read.excel <- function(header=TRUE,...) {read.table("clipboard",sep="\t",header=header, na.strings=c("","-","na"),...)}
# data2 = read.excel() #read clipboard from excel
# save(data2, file = file.path("./Rdata", "2022platy_meta.RData"))
load("./Rdata/2022platy_meta.RData")
data2
# write.csv(data2,'./Rdata/data2.csv')
# data2$BirthYear = 2022
data2 <- dplyr::select(data2, -c(age)) # remove column. Make sure have package on front
str(data2)
data2$ID <- as.factor(as.character(data2$ID))


# checking
CheckGeno(data1, Plot = T, Return = "excl")
# Warning:  There are 27516 SNPs scored for <50% of individuals, it is strongly advised to exclude those
#
# Warning:  *********** There are 9 individuals scored for <5% of SNPs, these WILL BE IGNORED ***********
#
#   Warning:  There are 20 individuals scored for <50% of SNPs, it is strongly advised to exclude those
#
# After exclusion, There are  70  individuals and  58754  SNPs.

# args.AP1 <- list(type = "flat", ages = c(0, 1))


SeqOUT1 <- sequoia(
  GenoM = data1,
  LifeHistData = data2,
  Err = 0.005, # genotyping error rate
  Module = "par",
  # args.AP = args.AP1,
  quiet = "verbose",
  Plot = T
)

GetMaybeRel(
  GenoM = data1,
  LifeHistData = data2,
  Err = 0.005, # genotyping error rate
  Module = "par",
  Herm = "B"
)

# For each pair of candidate relatives, the likelihoods are calculated of them being parent-offspring (PO), full siblings (FS), half siblings (HS),
# grandparent-grandoffspring (GG), full avuncular (niece/nephew - aunt/uncle; FA), half avuncular/great-grandparental/cousins (HA), or unrelated (U).
# Assignments are made if the likelihood ratio (LLR) between the focal relationship and the most likely alternative exceed the threshold Tassign.

# Warning:  There are 2752 SNPs scored for <50% of individuals
# After exclusion, There are  45  individuals and  3991  SNPs.
# it's possible that the SNP is rare in the population you are studying, or it might be located in a region of the genome that is difficult to
# sequence accurately.#it could also be clones.

SeqOUT1$DupGenotype
# this might indicate self fert
# Not all pairs flagged as potential duplicate genotypes are necessarily actual duplicates: inbred individuals
# may be nearly indistinguishable from their parent(s), especially when the number of SNPs is limited.


#  platy hiphop------------------------------------------------------------

load("./Rdata/2022platy_snp.RData") # loaded individual as rows
data1.1[1:50, 1:42] # check subset to know it is okay - full df doesn't display right, or use f2
library(tidyr)
data1.1[is.na(data1.1)] <- -9 # turn na to -9
data1 <- t(data1.1)
data1[1:50, 1:42]
# mat1 = as.matrix(data1) #use f2 to view


load("./Rdata/2022platy_meta.RData")
data2
# write.csv(data2,'./Rdata/data2.csv')
# data2$BirthYear = 2022
str(data2)
data2$ID <- as.factor(as.character(data2$ID))



# clonal assessemnt (nneeds finishing)-------------------------------------------------------
library(tidyr)
library(adegenet)
library(poppr)

# load SNP data (2 row locus; 0 and 1s)
# load('./Rdata/2022platy_snp.RData') #loaded individual as cols data1.1, loci as rows

# load 0,1,2 data (single row locus)
data1.1 <- read.csv(file = "./data/Report_DPlatyg23-7805_SNP_mapping_2_biallelic.csv", na = "-")
data1.1[data1.1 == "-"] <- NA
head(data1.1)



# data1.1[1:50, 1:42]   #check subset to know it is okay - full df doesn't display right, or use f2
missing_percent <- sapply(data1.1, function(x) sum(is.na(x)) / length(x) * 100)
cols_to_remove <- names(missing_percent[missing_percent > 70])
data1.1 <- data1.1[, !(names(data1.1) %in% cols_to_remove)]

data1.2 <- data1.1[, 2:ncol(data1.1)]
data1.2 <- t(data1.2) %>% data.frame()
data1 <- data1.2
# data1.1[is.na(data1.1)] <- -9  #turn na to -9
# data1 = t(data1.1) %>%  data.frame() #rows ID, cols allelles
str(data1)

# Convert SNP data to a genind object
(my_genind <- new("genind", ind.names = rownames(data1), loc.names = colnames(data1), tab = as.matrix(data1), ploidy = 2))


(pc <- as.genclone(my_genind))
# (mlg_table <- mlg(pc, quiet = FALSE))  #multilocus genotypes
# (mlg_tabl2 <- mlg.table(pc, bar = TRUE, quiet = FALSE))  #multilocus genotypes
# mlg.crosspop(pc, quiet = FALSE)

## Use mlg.filter to create a distance threshold to define multilocus genotypes.
mlg.filter(pc, threshold = 0.01, missing = "ignore") # 0.1 allows 10% dissimilarity
pc
