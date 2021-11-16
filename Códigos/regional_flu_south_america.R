library(phylodyn)
library(ggplot2)

data("regional_flu")
parameters <- c(2.544882, 1.203437)
results_matching_gamma <- BNPR(data = regional_flu$SouthAmerica, lengthout = 100,
                prec_alpha = parameters[1],
                prec_beta = 1/parameters[2])

results_pc_prior <- BNPR(data = regional_flu$SouthAmerica, lengthout = 100, 
     pc_prior = TRUE)

results_gamma_flat <- BNPR(data = regional_flu$SouthAmerica, lengthout =100, prec_alpha = 0.001, prec_beta = 0.001)

plot_BNPR(results_pc_prior, col = rgb(0.829, 0.680, 0.306),
          heatmap_labels_side = "left")

all_results <- data.frame(
  effpop975 <- c(results_gamma_flat$effpop975, results_matching_gamma$effpop975, results_pc_prior$effpop975),
  effpop025 <- c(results_gamma_flat$effpop025, results_matching_gamma$effpop025, results_pc_prior$effpop025), 
  effpop <- c(results_gamma_flat$effpop, results_matching_gamma$effpop, results_pc_prior$effpop),
  Prior = rep(c("Gamma Flat", "Matching Gamma", "PC prior"), each = 100),
  time <- c(results_gamma_flat$x, results_matching_gamma$x, results_pc_prior$x),
  effpopmean <- c(results_gamma_flat$effpopmean, results_matching_gamma$effpopmean, results_pc_prior$effpopmean)
)

options(repr.plot.width = 10, repr.plot.height = 7)
peffpop <- ggplot(data = all_results) +
  scale_x_continuous("Time", expand = c(0,0)) +
  scale_y_continuous("Effective population size", expand = c(0,0)) +
  geom_line(aes(x = time, y = effpop, col = Prior)) + 
  theme_bw() +
  scale_x_reverse() +
  geom_ribbon(data = all_results, aes(x = time, ymin= effpop025, ymax= effpop975,
                                      fill = Prior),alpha = 0.1) +
  geom_hline(yintercept = 1, linetype = "dotted")
  

peffpop
