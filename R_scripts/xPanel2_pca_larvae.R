
library(cowplot)
library(gridExtra)




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

ggsave(panel2, filename = 'fig2.pdf',  path = "./plots", device = 'pdf',  width = 9, height = 9.5)  #
