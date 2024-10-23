#Panel 2 - larval assignments

library(cowplot)
library(gridExtra)

##create an inset

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


panel2 = grid.arrange(
  arrangeGrob(p_combined, ncol = 1),
  arrangeGrob(t3, t4, ncol = 2),
  nrow = 2,
  heights = c(2, 1)  # Top row is double the height of the bottom row
)

#ggsave(p_combined, filename = 'fig2.tiff',  path = "./plots", device = 'tiff',  width = 8, height = 7)  #
ggsave(panel2, filename = 'fig2.pdf',  path = "./plots", device = 'pdf',  width = 9, height = 9.5)  #
