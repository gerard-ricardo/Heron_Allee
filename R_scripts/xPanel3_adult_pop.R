


library(cowplot)
p_combined1 = plot_grid(t2, struc_plot, p4, bathy_plot, ncol = 2)
p_combined1

ggsave(p_combined1, filename = 'fig3.pdf',  path = "./plots", device = 'pdf',  width = 8, height = 6)  #


