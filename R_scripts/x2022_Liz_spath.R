library(dartR)
library(adegenet)
library(HardyWeinberg)
library(ggplot2)






data_gl <- gl.read.dart(filename = "./data/Report_DAc23-7804_SNP_2 - Copy.csv", ind.metafile="./data/meta_acro_ordered2.csv",topskip=6)
meta_acro_final<-read.csv("./data/meta_acro_final.csv", head=T) #make sure samples are in same order as in data_gl


data_gl$n.loc
data_gl <-gl.drop.pop(data_gl,"tenuis", as.pop="species") #drop population. tenuis or spath
data_gl$ind.names

data_gl <- gl.reassign.pop(data_gl, as.pop="stage") 

meta_acro_spat_final<-meta_acro_final[match(data_gl$ind.names, meta_acro_final$id),]
nrow(meta_acro_spat_final)# check how many individuals remain after removing tenuis samples

data_gl<-gl.recalc.metrics(data_gl, v=3) #recalculate loci metrics
data_gl_filtered=data_gl

data_gl_filtered$other$loc.metrics$coverage<-data_gl_filtered$other$loc.metrics$AvgCountRef + data_gl_filtered$other$loc.metrics$AvgCountSnp
mean(data_gl_filtered$other$loc.metrics$coverage) 
min((data_gl_filtered$other$loc.metrics$coverage))
max((data_gl_filtered$other$loc.metrics$coverage))


data_gl_filtered1 <- gl.report.callrate(data_gl_filtered,method = "loc")
data_gl_filtered2 <- gl.report.callrate(data_gl_filtered,method = "ind")

data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v=3) #filter by loci callrate
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method="ind", threshold = 0.7, v=3) #filter by ind callrate
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t=0.7,v=3) #filter out loci with limited reproducibility
data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered,v=3) #remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
data_gl_filtered <- gl.filter.hwe(data_gl_filtered, alpha_val = 0.05, subset = "each", multi_comp_method = 'bonferroni',v=3) #filter out loci that depart from H-W proportions
data_gl_filtered<- gl.filter.secondaries(data_gl_filtered,method="random", verbose = 3) #remove loci fragment that shared SNPs. Only keep 1
list.match <- data_gl_filtered$loc.names[which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.05 & data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef > 0.05 & data_gl_filtered$other$loc.metrics$coverage > 10)] #remove loci based on minor allele frequency and low data coverage
data_gl_filtered<-data_gl_filtered[,match(list.match, data_gl_filtered$loc.names)]#keep only loci in the list above
meta_acro_spat_final_filtered<-meta_acro_spat_final[match(data_gl_filtered$ind.names, meta_acro_spat_final$id),] #match metadata file with genlight object

dim(meta_acro_spat_final_filtered)

meta_acro_spat_final_filtered$stage <- as.factor(as.character(meta_acro_spat_final_filtered$stage))
meta_acro_spat_final_filtered$genotype <- as.factor(as.character(meta_acro_spat_final_filtered$genotype))
meta_acro_spat_final_filtered$stage<-droplevels(meta_acro_spat_final_filtered$stage)
meta_acro_spat_final_filtered$genotype<-droplevels(meta_acro_spat_final_filtered$genotype)


data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop="genotype") 
data_gl_filtered

pca.data <- tab(data_gl_filtered, freq=TRUE, NA.method="mean")
pca <- dudi.pca(pca.data, center=T, scale=F, nf=2, scannf=FALSE)

pca_complete <-data.frame(pca$li, pop = data_gl_filtered$pop)
data1 = dplyr::arrange(pca_complete, Axis1) #dplyr - use this. Allows multiple sorts i.e site then orient

t <-ggplot(pca_complete,aes(x=Axis1,y=Axis2))+geom_point(aes(fill=pop),shape=21,size=3)+
  geom_text(aes(label = rownames(pca_complete)), hjust = 0, vjust = 0)+
  scale_fill_manual(values = c("green","red","purple", "pink", "blue", "orange", "brown", "black")) #number of colors representing the number of groups in the variable
t
