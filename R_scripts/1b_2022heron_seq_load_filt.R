


filter_data <- function(data, filter_type = "basic") {
  data_gl_filtered <- data
  data_genind <- NULL
  data_genind_adult <- NULL
  data_gl_filtered_adult <- NULL
  data_genind_progeny <-  NULL 

  data_gl$other$loc.metrics
  data_gl$other$loc.metrics$coverage <- data_gl$other$loc.metrics$AvgCountRef + data_gl$other$loc.metrics$AvgCountSnp
  median(data_gl$other$loc.metrics$coverage) #  PD = 11.54839
  min((data_gl$other$loc.metrics$coverage)) # PD = 5
  max((data_gl$other$loc.metrics$coverage)) #PD = 224.7
  sd(data_gl$other$loc.metrics$coverage) / sqrt(1996) # 0.134746
  
  if (filter_type == "basic" || filter_type == "medium") {
    
  data_gl_filtered <- data_gl
  
  
  gl.report.secondaries(data_gl_filtered)
  data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method="random", verbose = 3) #remove loci fragment that shared SNPs. Only keep 1
  
  gl.report.rdepth(data_gl_filtered)
  data_gl_filtered <- gl.filter.rdepth(data_gl_filtered,  lower = 10, v = 3) # filter by loci callrate
  
  gl.report.reproducibility(data_gl_filtered )
  data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t=0.95, v=3) #filter out loci with limited reproducibility
  
  gl.report.callrate(data_gl_filtered, method = "loc") 
  data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
  
  list.match <- data_gl_filtered$loc.names[
    which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.01 & 
            data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.99 & 
            data_gl_filtered$other$loc.metrics$OneRatioRef < 0.99 & 
            data_gl_filtered$other$loc.metrics$OneRatioRef > 0.01 & 
            data_gl_filtered$other$loc.metrics$coverage > 4)
  ]
  data_gl_filtered <- data_gl_filtered[, match(list.match, data_gl_filtered$loc.names)]
  
  data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v=3) #remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
  
  
  gl.report.callrate(data_gl_filtered, method = "ind") 
  pre_filt_ind <- data_gl_filtered@ind.names
  data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.58, v = 3) # filter by ind callrate
  filt_ind <- data_gl_filtered@ind.names
  (lost_ind <- setdiff(pre_filt_ind, filt_ind))
  length(filt_ind[grep('.l.', filt_ind )])  #count the larvae
  
  data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3) # recalculate loci metrics
  

  
  
  
  
  data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop = "genotype")
  data_gl_filtered
  
  data_genind <- gl2gi(data_gl_filtered)
  
  

  
  
  
  
  }
  
  return(list(data_gl_filtered = data_gl_filtered, data_gl_filtered_adult = data_gl_filtered_adult, data_genind = data_genind, 
              data_genind_adult = data_genind_adult, data_genind_progeny = data_genind_progeny))
}


filter_plus_null <- function(data_genind = data_genind) {
  
  
  
  num_loci <- nLoc(data_genind) # Get the number of loci in the genind object
  sampled_loci_indices <- sample(num_loci, num_loci) # Randomly sample x loci (max popgenreport can report)
  sampled_genind_obj <- data_genind[, sampled_loci_indices]
  pop(sampled_genind_obj) <- factor(rep("Combined_Population", nInd(sampled_genind_obj)))
  report1 = popgenreport(sampled_genind_obj, mk.null.all=TRUE, mk.pdf=FALSE)
  
  null_alleles_rep = report1$counts$nallelesbyloc
  null_alleles = colnames(null_alleles_rep)
  length(null_alleles)
  all_loci <- locNames(data_genind)
  length(all_loci)
  loci_to_keep <- setdiff(all_loci, null_alleles)
  data_genind <- data_genind[loc = loci_to_keep]
  data_genind@loc.n.all
  
  data_genind@other <- NULL
  data_genind <- new("genind",
                     tab = data_genind@tab,
                     pop = data_genind@pop,
                     ploidy = data_genind@ploidy,
                     loc.names = locNames(data_genind),
                     ind.names = indNames(data_genind),
                     strata = strata(data_genind))  # Include if you have stratification

  n_loci <- nLoc(data_genind)  # Should be 366 after subsetting
  
  
  data_gl_filtered = gi2gl(data_genind, parallel = FALSE, verbose = NULL)
  
  
  
  

  return(list(data_genind = data_genind,  data_gl_filtered = data_gl_filtered))

}




