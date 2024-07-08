#STRUCTURE

#alpha = degree of admixture. 0 = no mixing. > 1  = substantial mixing. 
#Fst  = between group variability (over total variability).  0 =  none, 1 = high, 0.25 = substantial



# 2000_10000 --------------------------------------------------------------

#note that probably need a burnin over 6000

# Load necessary libraries
library(ggplot2)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code


# Data from the provided table
data <- data.frame(
  K = 1:8,
  Ln_P_D = c(-543765.2, -543712.3, -543705.6, -544103.1, -543663.4, -543808.1, -543780.7, -543736.7),
  Var_Ln_P_D = c(8247.3, 8119.4, 8091.3, 8689.1, 8402.3, 8279.8, 8247.9, 8139.1),
  Alpha = c(NA, 6.6448, 8.7197, 9.028, 8.5596, 9.0107, 8.9091, 8.8645),
  Fst_1 = c(0.0068, 0.0094, 0.0102, 0.0093, 0.0099, 0.0099, 0.0092, 0.0098),
  Fst_2 = c(NA, 0.0089, 0.0103, 0.0095, 0.0099, 0.0099, 0.0093, 0.0087),
  Fst_3 = c(NA, NA, 0.0092, 0.0092, 0.0107, 0.0107, 0.0093, 0.0097),
  Fst_4 = c(NA, NA, NA, 0.0087, 0.0091, 0.0091, 0.0101, 0.0101),
  Fst_5 = c(NA, NA, NA, NA, 0.0097, 0.0097, 0.0104, 0.0105),
  Fst_6 = c(NA, NA, NA, NA, 0.0100, 0.0100, 0.0100, 0.0103),
  Fst_7 = c(NA, NA, NA, NA, NA, 0.0107, 0.0100, 0.0090),
  Fst_8 = c(NA, NA, NA, NA, NA, NA, NA, 0.0105)
)

# Plot Ln P(D) with K
p1 = ggplot(data, aes(x = K, y = Ln_P_D)) +
  geom_point() +
  geom_line() +
  labs(title = "Log Probability of Data (Ln P(D)) for Different K",
       x = "K (Number of Clusters)",
       y = "Ln P(D)")

# Plot Variance of Ln P(D) with K
p2 = ggplot(data, aes(x = K, y = Var_Ln_P_D)) +
  geom_point(color = "red") +
  geom_line(color = "red") +
  labs(title = "Variance of Ln P(D) for Different K",
       x = "K (Number of Clusters)",
       y = "Variance of Ln P(D)") 

# Plot Alpha with K
p3 = ggplot(data, aes(x = K, y = Alpha)) +
  geom_point(color = "green") +
  geom_line(color = "green") +
  labs(title = "Alpha for Different K",
       x = "K (Number of Clusters)",
       y = "Alpha") 

# Plot Fst values for different K
p4 = ggplot(data) +
  geom_point(aes(x = K, y = Fst_1, color = "Fst_1")) +
  geom_line(aes(x = K, y = Fst_1, color = "Fst_1")) +
  geom_point(aes(x = K, y = Fst_2, color = "Fst_2")) +
  geom_line(aes(x = K, y = Fst_2, color = "Fst_2")) +
  geom_point(aes(x = K, y = Fst_3, color = "Fst_3")) +
  geom_line(aes(x = K, y = Fst_3, color = "Fst_3")) +
  geom_point(aes(x = K, y = Fst_4, color = "Fst_4")) +
  geom_line(aes(x = K, y = Fst_4, color = "Fst_4")) +
  geom_point(aes(x = K, y = Fst_5, color = "Fst_5")) +
  geom_line(aes(x = K, y = Fst_5, color = "Fst_5")) +
  geom_point(aes(x = K, y = Fst_6, color = "Fst_6")) +
  geom_line(aes(x = K, y = Fst_6, color = "Fst_6")) +
  geom_point(aes(x = K, y = Fst_7, color = "Fst_7")) +
  geom_line(aes(x = K, y = Fst_7, color = "Fst_7")) +
  geom_point(aes(x = K, y = Fst_8, color = "Fst_8")) +
  geom_line(aes(x = K, y = Fst_8, color = "Fst_8")) +
  labs(title = "Fst Values for Different K",
       x = "K (Number of Clusters)",
       y = "Fst Values") +
  scale_color_discrete(name = "Fst Values") +
  theme(legend.position = "bottom")

#5 clusters best overall?

library(gridExtra)
grid.arrange(arrangeGrob(p1, p2, ncol = 2), arrangeGrob(p3, p4, ncol = 2), nrow = 2)
