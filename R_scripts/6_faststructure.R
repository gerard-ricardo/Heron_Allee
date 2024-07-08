# faststructure (works - only via Ubuntu)-----------------------------------------------------------

fast_struc <- gl.run.faststructure(data_gl_filtered_adult, num.k.rep = 2, k.range = 2:5, 
                                   exec = '/home/gricardo/fastStructure', output = '/home/gricardo/fastStructure1')
fast_struc
#something weird happends with only num.k.rep = 1, make sure to have two. Also does not seem to work for k = 1. 
#ggsave(filename = "fast_structure_plot.png", plot = fast_struc$plot, path = getwd())

#mixture plot (not working) - close but the issue is that gl.run.faststructure creates a bogus added value ontop of max k
max_k <- max(as.numeric(names(fast_struc$q_list)))
fast_struc$q_list <- fast_struc$q_list[names(fast_struc$q_list) != as.character(max_k)]
qmat <- gl.plot.faststructure(fast_struc, k.range = 2:3)
gl.map.structure(qmat, K = 2, data_gl_filtered_adult, scalex = 1, scaley = 0.5)


# Extract Q-values for K = 2, Replicate 1
# q_values_k2_rep1 <- fast_struc$q_list$`2`$`1`
# # Extract Q-values for K = 3, Replicate 1
# q_values_k3_rep1 <- fast_struc$q_list$`3`$`1`
# ggsave(filename = "fast_structure_plot2.png", plot = p_k2, path = getwd())
# ggsave(filename = "fast_structure_plot2.png", plot = p_k3, path = getwd())