##heron gamete output

#platy
colony_diam = 31   #n= 5
col_area <- (4/3 *pi * (colony_diam / 2)^3) * 0.8   #used volume minus 20%
polyp_den = 4   #very rough estimate from colony 13 setting
bundles_col = col_area * polyp_den
fecun_zone = 0.7
fecund <-  bundles_col * fecun_zone

#sperm and eggs
sc = 5*10^6  #per bundle
sperm_col = sc * fecund   #total sperm
eggs = 53
egg_col = eggs * fecund
