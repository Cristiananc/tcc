library(phylodyn)
library(ggplot2)

data("regional_flu")

#Considering Pr(tau > 10^-5) = 10^-4 
parameters <- c(0.8007825, 1.0802954)

results <- BNPR(data = regional_flu$SouthAmerica, lengthout = 100,
                prec_alpha = parameters[1],
                prec_beta = 1/parameters[2])

plot_BNPR(results, col = rgb(0.829, 0.680, 0.306),
          heatmap_labels_side = "left")

#Cyclic trajectory
set.seed(9)
trajectory = cyclic_traj
samp_times = sort(runif(100, 0, 20))
n = rep(1,100)

#Simulating a genealogy based on the sample
gene = coalsim(samp_times = samp_times, n_sampled = n,
               traj = trajectory, lower_bound = 1/10)

#Gamma flat
results1 <- BNPR(data = gene, lengthout =100, prec_alpha = 0.001, prec_beta = 0.001)

#Gamma 2
par_ab <- c(2.544882, 1.203437)
results2 <- BNPR(data = gene, lengthout =100, 
                 prec_alpha = par_ab[1], prec_beta = 1/par_ab[2])

#PC prior
results3 <- BNPR(data = gene, lengthout = 100, 
         pc_prior = TRUE)

#Plot
require(gridExtra)
options(repr.plot.width = 10, repr.plot.height = 5)

par(mfrow=c(1, 3))
plot_BNPR(results1, col = rgb(0.829, 0.680, 0.306),
          heatmap_labels_side = "left", traj = cyclic_traj)

plot_BNPR(results2, col = rgb(0.829, 0.680, 0.306),
          heatmap_labels_side = "left", traj = cyclic_traj)

plot_BNPR(results3, col = rgb(0.829, 0.680, 0.306),
          heatmap_labels_side = "left", traj = cyclic_traj)
