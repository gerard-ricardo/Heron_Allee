

load("./Rdata/heron_intercol_all.RData")  #p1
load("./Rdata/heron_density.RData")  #p2
load("./Rdata/heron_adult_site.RData") #p3






library(cowplot)


p_combined <- ggdraw() +
  draw_plot(p3) +  # The main plot (p3)
  draw_plot(p1, x = 0.15, y = 0.15, width = 0.35, height = 0.35) +  # Lower left inset (p1)
  draw_plot(p2, x = 0.66, y = 0.55, width = 0.35, height = 0.35)  # Upper right inset (p2)
p_combined





ggsave(p_combined, filename = 'fig1.tiff',  path = "./plots", device = 'tiff',  width = 7.1, height = 7.5)  #

ggsave(p_combined, filename = 'fig1.pdf',  path = "./plots", device = 'pdf',  width = 7.1, height = 7.5)  #


