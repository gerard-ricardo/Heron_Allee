#DAPC


# tutorial https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html----------------------------------------------------------------


#To use model-free methods K-means clustering based on genetic distance or DAPC are more convenient approaches for populations that 
#are clonal or partially clonal.

# DAPC requires the adegenet package. Let's load this package:
library("adegenet")
library("poppr")
library("magrittr")
data(H3N2) # load the H3N2 influenza data. Type ?H3N2 for more info.
pop(H3N2) <- H3N2$other$epid
dapc.H3N2 <- dapc(H3N2, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(H3N2) - 1)
scatter(dapc.H3N2, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

set.seed(4)
contrib <- loadingplot(dapc.H3N2$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)

temp    <- seploc(H3N2)       # seploc {adegenet} creates a list of individual loci.
snp906  <- tab(temp[["906"]]) # tab {adegenet} returns a matrix of genotypes
snp399  <- tab(temp[["399"]])

# The following two commands find the average allele frequencies per population
(freq906 <- apply(snp906, 2, function(e) tapply(e, pop(H3N2), mean, na.rm = TRUE)))
(freq399 <- apply(snp399, 2, function(e) tapply(e, pop(H3N2), mean, na.rm = TRUE)))

# First, set the plotting parameters
# mfrow = number of columns, rows
# mar   = plot margin size
# las   = axis label style (3: always vertical)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 0) + 0.1, las = 3)

matplot(freq906,  pch = c("c", "t"), type = "b",
        xlab = "year", ylab = "allele frequency", main = "SNP # 906",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:6, lab = 2001:2006)

matplot(freq399, pch = c("c", "t"), type = "b",
        xlab = "year", ylab = "allele frequency", main = "SNP #399",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:6, lab = 2001:2006)


###crossval to determine number of dimensions
library("poppr")
data("Pram")
set.seed(999)
pramx <- xvalDapc(tab(Pram, NA.method = "mean"), pop(Pram))
set.seed(999)
system.time(pramx <- xvalDapc(tab(Pram, NA.method = "mean"), pop(Pram),
                              n.pca = 10:20, n.rep = 1000,
                              parallel = "snow", ncpus = 4L))
names(pramx) # The first element are all the samples
pramx[-1]
scatter(pramx$DAPC, col = other(Pram)$comparePal, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

# my data -----------------------------------------------------------------

# Perform DAPC
dapc_genind_adult <- dapc(data_genind_adult, var.contrib = TRUE, scale = FALSE, n.pca = 3, n.da = nPop(data_genind_adult) - 1)
# Scatter plot of the DAPC results
scatter(dapc_genind_adult, cell = 0, pch = 18:23, cstar = 0, mstree = F, lwd = 5, lty = 5)



set.seed(4)
contrib <- loadingplot(dapc_genind_adult$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)



### Perform cross-validation to determine the optimal number of principal components  (not working)
set.seed(999)
xval_result <- xvalDapc(tab(data_genind_adult, NA.method = "mean"), pop(data_genind_adult))

# If parallel computation is desired
set.seed(999)
system.time(xval_result <- xvalDapc(tab(data_genind_adult, NA.method = "mean"), pop(data_genind_adult),
                                    n.pca = 10:20, n.rep = 1000,
                                    parallel = "snow", ncpus = 4L))

# Examine the cross-validation result
names(xval_result) # The first element contains all the samples
xval_result[-1]

# Scatter plot of the DAPC results
scatter(xval_result$DAPC, col = rainbow(nPop(data_genind_adult)), cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
