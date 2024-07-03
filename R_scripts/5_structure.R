# structure analysis (dart - working) ----------------------------------------------------
# tutorial https://green-striped-gecko.github.io/kioloa/session12.html
#test run
struct_adult = gl.run.structure(data_gl_filtered_adult, verbose = 3, burnin = 1000, numreps = 1000, k.range = 2:5, num.k.rep = 2, 
                                seed = 1, noadmix=FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
#linux
# tic("Run Structure") # Start the timer
# struct_adult = gl.run.structure(data_gl_filtered_adult, verbose = 3, burnin = 1000, numreps = 1000, k.range = 2:5, num.k.rep = 2, 
#                                 seed = 1, noadmix=FALSE, exec = "/home/gricardo/structure/structure.exe")
# toc() # End the timer (170.7 sec   = is slower)

# formal run
struct_adult = gl.run.structure(data_gl_filtered_adult, verbose = 3, burnin = 10000, numreps = 100000, k.range = 2:4, num.k.rep = 10, 
                                seed = 1, noadmix = FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
#seems to work for k.range = 2:4, num.k.rep = 2 (but not other settings)
#save(struct_adult, file = file.path("./Rdata", "struct_adult_1_10.RData"))
load("./Rdata/struct_adult_1_10.RData")
#seems to vary each time
str(struct_adult)
ev <- gl.evanno(struct_adult, plot.out = TRUE)

qmat <- gl.plot.structure(struct_adult, K = 3)
head(qmat)

#create a map showing groupings
gl.map.structure(qmat = qmat, x = data_gl_filtered_adult, K = 3, scalex = 1, scaley = 0.5)