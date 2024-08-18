#Panel 1


load("./Rdata/heron_intercol_all.RData")  #p1
load("./Rdata/heron_density.RData")  #p2
load("./Rdata/heron_adult_site.RData") #p3




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
  draw_plot(p3) +  # The main plot (p3)
  draw_plot(p1, x = 0.15, y = 0.15, width = 0.35, height = 0.35) +  # Lower left inset (p1)
  draw_plot(p2, x = 0.66, y = 0.55, width = 0.35, height = 0.35)  # Upper right inset (p2)
p_combined

#ggsave(filename = "p_combined.png", plot = p_combined, path = "./plots", device = "png", width = 7.1, height = 7.5)


#final_plot <- plot_grid(p_combined_fixed, p4, ncol = 1, rel_heights = c(1, 1))


ggsave(p_combined, filename = 'fig1.tiff',  path = "./plots", device = 'tiff',  width = 7.1, height = 7.5)  #

ggsave(p_combined, filename = 'fig1.pdf',  path = "./plots", device = 'pdf',  width = 7.1, height = 7.5)  #


