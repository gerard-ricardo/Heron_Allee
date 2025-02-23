struct_adult = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 2000, numreps = 3000, k.range = 1:5, num.k.rep = 2, 
                                seed = 1, noadmix=FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")

struct_all = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 2000, numreps = 2000, k.range = 1:5, num.k.rep = 2, 
                                seed = 1, noadmix=FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")

tic("Running structure analysis") # start the timer with a message
struct_adult = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 8000, numreps = 20000, k.range = 1:5, num.k.rep = 2, 
                                seed = 1, noadmix = FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
toc() 

tic("Running structure analysis") # start the timer with a message
struct_adult = gl.run.structure(data_gl_adult_unique, verbose = 3, burnin = 10000, numreps = 50000, k.range = 2:4, num.k.rep = 3, 
                                seed = 1, noadmix = FALSE, exec = "C:/Users/gerar/OneDrive/Documents/structure/structure.exe")
toc() 

load("./Rdata/struct_adult_1_10.RData")
str(struct_adult)
ev <- gl.evanno(struct_adult, plot.out = TRUE)



qmat <- dartR::gl.plot.structure(struct_adult, K = 3, colors_clusters = list("grey", "grey", "grey", 'grey'), clumpak = T, save2tmp = T)


Q_melt <- do.call("rbind", lapply(qmat, reshape2::melt, id.vars = c("Label", "K", "orig.pop", "ord"), variable.name = "Cluster" ))

Q_melt$orig.pop <- factor(Q_melt$orig.pop, levels = unique(struct_adult[[1]]$q.mat$orig.pop))
Q_melt$Cluster <- factor(Q_melt$Cluster, levels = c("cluster1", "cluster2", "cluster3"))  # Set consistent factor levels
levels(Q_melt$Cluster)  # Check factor levels


struc_plot3 <- ggplot(Q_melt, aes(x = factor(ord), y = value, fill = Cluster)) +
  geom_col(color = "black", size = 0.25, width = 1) +
  facet_grid(K ~ orig.pop, scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(
    breaks = unique(Q_melt$ord), labels = unique(Q_melt$orig.pop), expand = c(0, 0)
  ) +
  scale_fill_manual(values = c(cluster1 = "mediumseagreen", cluster2 = "salmon", cluster3 = "dodgerblue" )) +  # Colour palette
  theme_sleek2() +
  theme(
    panel.spacing = unit(0, "lines"), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),  # Remove x-axis facet labels (orig.pop)
    strip.text.y = element_blank(),  # Remove y-axis facet labels (K)
    axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = 'Individuals', y = 'K=3')
struc_plot3


qmat2 <- dartR::gl.plot.structure(struct_adult, K = 2, colors_clusters = list("grey", "grey", "grey", 'grey'), clumpak = T, save2tmp = T)
Q_melt2 <- do.call("rbind", lapply(qmat2, reshape2::melt, id.vars = c("Label", "K", "orig.pop", "ord"), variable.name = "Cluster" ))

Q_melt2$orig.pop <- factor(Q_melt2$orig.pop, levels = unique(struct_adult[[1]]$q.mat$orig.pop))
Q_melt2$Cluster <- factor(Q_melt2$Cluster, levels = c("cluster1", "cluster2", "cluster3"))  # Set consistent factor levels
levels(Q_melt2$Cluster)  # Check factor levels


struc_plot2 <- ggplot(Q_melt2, aes(x = factor(ord), y = value, fill = Cluster)) +
  geom_col(color = "black", size = 0.25, width = 1) +
  facet_grid(K ~ orig.pop, scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(
    breaks = unique(Q_melt2$ord), labels = unique(Q_melt2$orig.pop), expand = c(0, 0)
  ) +
  scale_fill_manual(values = c(cluster1 = "dodgerblue", cluster2 = "salmon", cluster3 = "dodgerblue" )) +  # Colour palette
  theme_sleek2() +
  theme(
    panel.spacing = unit(0, "lines"), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_blank(),
    strip.text.x = element_blank(),  # Remove x-axis facet labels (orig.pop)
    strip.text.y = element_blank(),  # Remove y-axis facet labels (K)
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = 'K=2')
struc_plot2

library(gridExtra)

struc_plot = grid.arrange(struc_plot2, struc_plot3, ncol = 1, heights = c(1, 1.2))  # Arrange the two plots side by side



