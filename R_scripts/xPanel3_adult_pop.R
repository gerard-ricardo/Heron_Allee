#Panel Fig 3 genetics adult

# # PCA adult
# load("./Rdata/heron_adult_pca.RData")  #t2
# # Structure
# load("./Rdata/structure_plot.RData") #p3
# # Spatial auto
# load("./Rdata/spatial_auto_plot.RData") #p4
# #Depth by clustering
# load("./Rdata/bath_cluster.RData")  #p5


library(cowplot)
p_combined1 = plot_grid(t2, struc_plot, p4, bathy_plot, ncol = 2)
p_combined1

ggsave(p_combined1, filename = 'fig3.pdf',  path = "./plots", device = 'pdf',  width = 8, height = 6)  #


# library(gridExtra)
# gs = list(t2, p3, p4, p5)
# plots1 = grid.arrange(grobs = gs)
