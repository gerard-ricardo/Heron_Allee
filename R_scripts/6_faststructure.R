
fast_struc <- gl.run.faststructure(data_gl_filtered_adult, num.k.rep = 2, k.range = 2:5, 
                                   exec = '/home/gricardo/fastStructure', output = '/home/gricardo/fastStructure1')
fast_struc

max_k <- max(as.numeric(names(fast_struc$q_list)))
fast_struc$q_list <- fast_struc$q_list[names(fast_struc$q_list) != as.character(max_k)]
qmat <- gl.plot.faststructure(fast_struc, k.range = 2:3)
gl.map.structure(qmat, K = 2, data_gl_filtered_adult, scalex = 1, scaley = 0.5)


