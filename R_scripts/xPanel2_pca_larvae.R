#Panel 2 - larval assignments


load("./Rdata/scree_larvaue.RData")  #p0
load("./Rdata/pca_larvae.RData")  #t2





# library(gridExtra)
# p4 = grid.arrange(arrangeGrob(p1, p2, ncol = 2))


##create an inset
library(cowplot)
# p_combined_top <- ggdraw() +
#   draw_plot(p1) +  # The main plot
#   draw_plot(p2, x = 0.55, y = 0.45, width = 0.5, height = 0.5)  # The inset plot
# 
# p_final <- plot_grid(p_combined_top, p3, ncol = 1, rel_heights = c(1, 2))
# p_final


# Create the main plot with p3
p_combined <- ggdraw() +
  draw_plot(t2) +  # The main plot (p3)
  draw_plot(p0, x = 0.1, y = 0.65, width = 0.33, height = 0.33)   # Lower left inset (p1)
p_combined

ggsave(p_combined, filename = 'fig2.tiff',  path = "./plots", device = 'tiff',  width = 8, height = 7)  #
