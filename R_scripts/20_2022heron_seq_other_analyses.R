










outflank_results  = gl.outflank(data_genind_adult, plot = TRUE, RightTrimFraction = 0.12)
OutlierFlag = outflank_results$outflank$results$OutlierFlag
summary(OutlierFlag)








gl.report.parent.offspring(data_gl_filtered)  #nothing picked up with this basic analysis



res = gl.grm(data_gl_filtered) # relatedness matrix
res2 <- gl.grm.network(res, data_gl_filtered, relatedness_factor = 0.125)











data_gl <- gl.read.dart(filename = "./data/Report_DAc23-7804_SNP_2 - Copy.csv", ind.metafile = "./data/meta_acro_ordered2.csv", topskip = 6)


data_gl$n.loc
data_gl$ind.names



data_gl <- gl.drop.pop(data_gl, "spath", as.pop = "species") # DROP population. tenuis or spath


data_gl <- gl.reassign.pop(data_gl, as.pop = "stage")




data_gl <- gl.recalc.metrics(data_gl, v = 3) # recalculate loci metrics

data_gl_filtered <- data_gl


str(data_gl_filtered$other$loc.metrics)

data_gl_filtered$other$loc.metrics$coverage <- data_gl_filtered$other$loc.metrics$AvgCountRef + data_gl_filtered$other$loc.metrics$AvgCountSnp

plot(density(data_gl_filtered$other$loc.metrics$coverage))
quantile(data_gl_filtered$other$loc.metrics$coverage)
data_gl_filtered$other$loc.metrics$coverage <- data_gl_filtered$other$loc.metrics$AvgCountRef + data_gl_filtered$other$loc.metrics$AvgCountSnp
mean(data_gl_filtered$other$loc.metrics$coverage) # 61.59835. PD = 15.6
median(data_gl_filtered$other$loc.metrics$coverage) # 28.60925. PD = 15.6

min((data_gl_filtered$other$loc.metrics$coverage)) # 5.375. PD = 10
max((data_gl_filtered$other$loc.metrics$coverage)) # 974.6452. PD = 29
sd(data_gl_filtered$other$loc.metrics$coverage) / sqrt(1996) # 1.858626
coverage_values <- data_gl_filtered$other$loc.metrics$coverage
less_than_30 <- coverage_values < 30
(proportion_less_than_30 <- sum(less_than_30) / length(coverage_values))




plot(density(data_gl_filtered$other$loc.metrics$CallRate))
quantile(data_gl_filtered$other$loc.metrics$CallRate)
length(data_gl_filtered$other$loc.metrics$CallRate)
data_gl_filtered1 <- gl.report.callrate(data_gl_filtered, method = "loc") # acro 30%.   PLaty = 52%, 49% missing data, but ~57%
data_gl_filtered2 <- gl.report.callrate(data_gl_filtered, method = "ind") # so adults generally around 60% but variable, larvae about 60 but constant

data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.5, v = 3) # filter by ind callrate
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate



data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.7, v = 3) # filter by ind callrate
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t = 0.7, v = 3) # filter out loci with limited reproducibility
data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v = 3) # remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
data_gl_filtered <- gl.filter.hwe(data_gl_filtered, alpha_val = 0.05, subset = "each", multi_comp_method = "bonferroni", v = 3) # filter out loci that depart from H-W proportions
data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method = "random", verbose = 3) # remove loci fragment that shared SNPs. Only keep 1
list.match <- data_gl_filtered$loc.names[which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.05 & data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef > 0.05 & data_gl_filtered$other$loc.metrics$coverage > 10)] # remove loci based on minor allele frequency and low data coverage
data_gl_filtered <- data_gl_filtered[, match(list.match, data_gl_filtered$loc.names)] # keep only loci in the list above


meta_acro_spat_final_filtered <- meta_acro_spat_final[match(data_gl_filtered$ind.names, meta_acro_spat_final$id), ] # match metadata file with genlight object

dim(meta_acro_spat_final_filtered)

meta_acro_spat_final_filtered$stage <- as.factor(as.character(meta_acro_spat_final_filtered$stage))
meta_acro_spat_final_filtered$genotype <- as.factor(as.character(meta_acro_spat_final_filtered$genotype))
meta_acro_spat_final_filtered$stage <- droplevels(meta_acro_spat_final_filtered$stage)
meta_acro_spat_final_filtered$genotype <- droplevels(meta_acro_spat_final_filtered$genotype)

data_gl_filtered$pop <- meta_allsamples$stage
data_gl_filtered$pop <- meta_allsamples$genotype

data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop = "genotype")
data_gl_filtered

pca_data <- tab(data_gl_filtered, freq = TRUE, NA.method = "mean")
pca <- dudi.pca(pca_data, center = T, scale = F, nf = 2, scannf = FALSE)

pca_complete <- data.frame(pca$li, pop = data_gl_filtered$pop)


pca_complete <- pca_complete %>%
  mutate(
    Stage = ifelse(str_detect(row.names(pca_complete), "\\.a\\."), "Adult", "Larva"),
    MumID = str_extract(row.names(pca_complete), "(?<=pd)\\d+"),
    RepID = str_extract(row.names(pca_complete), "(?<=\\.)\\d+$"),
    NewID = paste0(Stage, MumID, "_", RepID)
  )

data1 <- dplyr::arrange(pca_complete, Axis1) # dplyr - use this. Allows multiple sorts i.e site then orient

t <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = pop), shape = 21, size = 3) +
  geom_text(aes(label = rownames(pca_complete)), hjust = 0, vjust = 0) +
  scale_fill_manual(values = c("green", "red", "purple", "pink", "blue", "orange", "brown", "black")) # number of colors representing the number of groups in the variable
t

t1 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = pop), shape = 21, size = 3) +
  geom_text(aes(label = rownames(pca_complete)), hjust = 0, vjust = 0) +
  scale_fill_manual(values = c("green", "red", "purple", "pink", "blue", "orange", "brown", "steelblue", "gold", "yellow", "grey", "black")) # number of colors representing the number of groups in the variable
t1


my_palette <- c(
  "dodgerblue", "salmon", "mediumseagreen", "orchid", "gold", "lightcoral",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato"
)




total_variance <- sum(pca$eig)
var_explained_axis1 <- pca$eig[1] / total_variance
var_explained_axis2 <- pca$eig[2] / total_variance
cat("Variance explained by Axis1: ", var_explained_axis1 * 100, "%\n")
cat("Variance explained by Axis2: ", var_explained_axis2 * 100, "%\n")


pca$eig # first 4 >1
pve <- pca$eig / sum(pca$eig)

plot(pve,
  type = "b",
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained",
  main = "Scree Plot",
  ylim = c(0, max(pve)),
  yaxs = "i"
)

lines(cumsum(pve), type = "b", pch = 19, col = "red")
legend("topright",
  legend = c("PVE", "Cumulative PVE"),
  col = c("black", "red"), lty = 1:1, cex = 0.8
)




library(LEA)
library(adegenet)
data_genlight <- glDart(data_gl_filtered)
X.snp <- df2genind(data_genlight, pop = data_genlight@pop)
res.admixture <- snmf(X.snp, K = 3, entropy = TRUE, repetitions = 10)

barplot(res.admixture$Q)



library(Relatedness)


genind_matrix <- as.matrix(data_genind_filtered@tab)

allele_frequencies <- colMeans(genind_matrix, na.rm = TRUE) / 2

frequency_matrix <- matrix(allele_frequencies, nrow = 1, ncol = length(allele_frequencies), byrow = TRUE)

if (ncol(genind_matrix) != ncol(frequency_matrix)) {
  stop("The number of columns in the genotype matrix and frequency matrix must match.")
}

relatedness_results <- RelCoef(IndividualGenom = genind_matrix, Freq = frequency_matrix)
print(relatedness_results$R)


dist_matrix <- bitwise.dist(data_gl_filtered)
dist_matrix_mat <- as.matrix(dist_matrix)
closest_individual1 <- character(nrow(dist_matrix_mat))
closest_individual2 <- character(nrow(dist_matrix_mat))

for (i in 1:nrow(dist_matrix_mat)) {
  dist_matrix_mat[i, i] <- Inf
  min_indices <- order(dist_matrix_mat[i, ])[1:2]
  closest_individual1[i] <- rownames(dist_matrix_mat)[min_indices[1]]
  closest_individual2[i] <- rownames(dist_matrix_mat)[min_indices[2]]
}

results <- data.frame(
  Individual = rownames(dist_matrix_mat),
  ClosestIndividual1 = closest_individual1,
  ClosestIndividual2 = closest_individual2
)

print(results)


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


gl1 <- gl.read.dart(filename = "./data/Report_DAc23-7804_SNP_2 - Copy.csv", ind.metafile = "./data/meta_acro.csv") # note,
spp <- gl1@other$ind.metrics$species


gl1@other$loc.metrics$RepAvg[1:10] ## Print the first 10 values of the RepAvg column in the loc.metrics object of gl1

gl1@ind.names
gl1@loc.names ## Print the names of the loci (SNPs) in gl1
gl1@gen[[1]] ## Print the genotype data for the first individual in gl1
nLoc(gl1) # number of loci = 1996

library(HardyWeinberg)
gl2 <- gl.filter.callrate(gl1, threshold = 0.8)
gl.diagnostics.hwe(
  x = gl.filter.allna(gl1[, 1:50]),
  stdErr = FALSE, n.cores = 1
)

gl1 <- gl.filter.monomorphs(gl1) # Filter out monomorphic sites (loci where all individuals have the same genotype)

missing_inds <- which(is.na(gl1@gen))

gl.dist.pop(gl1[1:45, 1:100]) # distance matrix  , genetic distances between individuals.
str(gl1)

pca1 <- gl.pcoa(gl2, verbose = 2)
gl.pcoa.plot(pca1, gl2) # pcoa  Fig. 1a. # The gl.pcoa.plot function creates a scatterplot of the first two principal coordinates,

gl.tree.nj(gl1) ## Construct a neighbor-joining tree from the genotype data in gl1







data1 <- read.csv(
  skip = 6, header = FALSE,
  file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Report-DAc23-7804", "Report_DAc23-7804_SNP_2.csv")
)
head(data1)
colnames(data1) <- as.character(data1[1, ]) # Set the names of data1 to be the values from the sixth row
data1 <- data1[-1, ] # Remove the row that was used as header
data1 <- type.convert(data1, as.is = TRUE) # Convert data frame columns to the correct data type

str(data1)

data1$CloneID <- as.factor(as.character(data1$CloneID))

data2 <- dplyr::select(data1, c(CloneID, SNP, CallRate), starts_with(c("sp", "at")))

data3 <- data2 %>%
  group_by(CloneID) %>%
  summarise(
    Last_SNP = last(SNP),
    Ref = substr(Last_SNP, start = nchar(Last_SNP) - 2, stop = nchar(Last_SNP) - 2),
    Var = substr(Last_SNP, start = nchar(Last_SNP), stop = nchar(Last_SNP))
  )

data4 <- data3 %>% pivot_longer(cols = -c(CloneID, Last_SNP), names_to = "allele", values_to = "base")

data4$RowID <- rep(c(1, 2), nrow(data4) / 2)
data2$RowID <- rep(c(1, 2), nrow(data2) / 2)

data5 <- left_join(data2, data4, by = c("CloneID", "RowID"))

data6 <- data5 %>%
  mutate(LocusID = CloneID) %>%
  dplyr::select(., c(LocusID, allele, base), starts_with(c("sp", "at")))

data7 <- dplyr::select(data6, starts_with(c("sp", "at"))) %>% t()

data8 <- data6 %>%
  group_by(LocusID) %>%
  mutate(LocusID_mod = ifelse(row_number() == 1, paste0(LocusID, "a"), paste0(LocusID, "b"))) %>%
  ungroup() %>%
  data.frame() %>%
  select(., LocusID_mod)

colnames(data7) <- data8$LocusID_mod
data7$

  data7[1:10, 1:10]









library(sequoia)

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



load("./Rdata/2022acro_meta.RData")
data2
data2$BirthYear <- 2022
data2 <- dplyr::select(data2, -c(age)) # remove column. Make sure have package on front
str(data2)
data2$ID <- as.factor(as.character(data2$ID))
data2$spp <- substr(data2$ID, 1, 2) # extract values from character string
data2_spp <- split(data2, data2$spp)
data2_at <- data2_spp$at
data2_sp <- data2_spp$sp



CheckGeno(data1_at, Plot = T)

CheckGeno(data1_sp, Plot = T)




SeqOUT1_at <- sequoia(
  GenoM = data1_at,
  LifeHistData = data2_at,
  Err = 0.005, # genotyping error rate
  Module = "par",
  quiet = "verbose",
  Plot = T
)


SeqOUT1_sp <- sequoia(
  GenoM = data1_sp,
  LifeHistData = data2_sp,
  Err = 0.005, # genotyping error rate
  Module = "par",
  quiet = "verbose",
  Plot = T
)





library(sequoia)
load("./Rdata/2022platy_snp.RData") # loaded individual as rows
data1.1[1:50, 1:42] # check subset to know it is okay - full df doesn't display right, or use f2
library(tidyr)
data1.1[is.na(data1.1)] <- -9 # turn na to -9
data1 <- t(data1.1)
data1[1:50, 1:42]


load("./Rdata/2022platy_meta.RData")
data2
data2 <- dplyr::select(data2, -c(age)) # remove column. Make sure have package on front
str(data2)
data2$ID <- as.factor(as.character(data2$ID))


CheckGeno(data1, Plot = T, Return = "excl")



SeqOUT1 <- sequoia(
  GenoM = data1,
  LifeHistData = data2,
  Err = 0.005, # genotyping error rate
  Module = "par",
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



SeqOUT1$DupGenotype



load("./Rdata/2022platy_snp.RData") # loaded individual as rows
data1.1[1:50, 1:42] # check subset to know it is okay - full df doesn't display right, or use f2
library(tidyr)
data1.1[is.na(data1.1)] <- -9 # turn na to -9
data1 <- t(data1.1)
data1[1:50, 1:42]


load("./Rdata/2022platy_meta.RData")
data2
str(data2)
data2$ID <- as.factor(as.character(data2$ID))



library(tidyr)
library(adegenet)
library(poppr)


data1.1 <- read.csv(file = "./data/Report_DPlatyg23-7805_SNP_mapping_2_biallelic.csv", na = "-")
data1.1[data1.1 == "-"] <- NA
head(data1.1)



missing_percent <- sapply(data1.1, function(x) sum(is.na(x)) / length(x) * 100)
cols_to_remove <- names(missing_percent[missing_percent > 70])
data1.1 <- data1.1[, !(names(data1.1) %in% cols_to_remove)]

data1.2 <- data1.1[, 2:ncol(data1.1)]
data1.2 <- t(data1.2) %>% data.frame()
data1 <- data1.2
str(data1)

(my_genind <- new("genind", ind.names = rownames(data1), loc.names = colnames(data1), tab = as.matrix(data1), ploidy = 2))


(pc <- as.genclone(my_genind))

mlg.filter(pc, threshold = 0.01, missing = "ignore") # 0.1 allows 10% dissimilarity
pc
