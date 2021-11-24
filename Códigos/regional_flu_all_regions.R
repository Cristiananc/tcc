library(phylodyn)
library(ggplot2)
library(ggpubr)

data("regional_flu")
set.seed(29)

parameters <- c(1.13340147e+02, 1.00071471e-02) #c(2.544882, 1.203437)

zero_dates <- list(USACanada = 2012.301, Europe = 2011.044, NorthChina = 2011.285,
                   JapanKorea = 2012.29, India = 2010.814, SouthChina = 2011.282,
                   SouthAmerica = 2011.518, SoutheastAsia = 2011.995, Oceania = 2010.964)

years = list(USACanada = 12, Europe = 11, NorthChina = 10,
             JapanKorea = 12, India = 10, SouthChina = 11,
             SouthAmerica = 11, SoutheastAsia = 11, Oceania = 10)

start <- 8/12

axlabs <- list(x = seq(1, 0, by=-1/12),
               labs = c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar",
                        "Apr", "May", "Jun", "Jul", "Aug", "Sep"))

#Verificar como plotar pra incluir os dados de meses

#Resultados para cada região
regions <- list("USACanada", "Europe", "NorthChina", "JapanKorea", "India", "SouthChina", "SouthAmerica", "SoutheastAsia", "Oceania")
results <- list()
for (region in regions){
  results[[region]] <- list(bnpr_matching_gamma = BNPR(data = regional_flu[[region]], lengthout = 100,
                                                       prec_alpha = parameters[1],
                                                       prec_beta = 1/parameters[2]),
                            bnpr_pc_prior = BNPR(data = regional_flu[[region]], lengthout = 100, 
                                                 pc_prior = TRUE),
                            bnpr_gamma_flat = BNPR(data = regional_flu[[region]], lengthout =100, 
                                                   prec_alpha = 0.001, 
                                                   prec_beta = 0.001))
}

results_plot <- list()

for (region in regions){
  results_plot[[region]] <- data.frame(
    effpop975 <- c(results[[region]]$bnpr_gamma_flat$effpop975, 
                   results[[region]]$bnpr_matching_gamma$effpop975, 
                   results[[region]]$bnpr_pc_prior$effpop975),
    
    effpop025 <- c(results[[region]]$bnpr_gamma_flat$effpop025, 
                   results[[region]]$bnpr_matching_gamma$effpop025, 
                   results[[region]]$bnpr_pc_prior$effpop025), 
    
    effpop <- c(results[[region]]$bnpr_gamma_flat$effpop, 
                results[[region]]$bnpr_matching_gamma$effpop, 
                results[[region]]$bnpr_pc_prior$effpop),
    
    Prior = rep(c("Gamma Flat", "Matching Gamma", "PC prior"), each = 100),
    
    Time <- c(results[[region]]$bnpr_gamma_flat$x, 
               results[[region]]$bnpr_matching_gamma$x, 
               results[[region]]$bnpr_pc_prior$x),
    
    effpopmean <- c(results[[region]]$bnpr_gamma_flat$effpopmean, 
                    results[[region]]$bnpr_matching_gamma$effpopmean, 
                    results[[region]]$bnpr_pc_prior$effpopmean),
    
    region_name <- rep(region, 300)
  )
}

peffpop <- rbind(results_plot[[1]], results_plot[[2]], results_plot[[3]], 
      results_plot[[4]], results_plot[[5]], results_plot[[6]],
      results_plot[[7]], results_plot[[8]], results_plot[[9]])

colnames(peffpop) <- c("effpop975", "effpop025", "effpop", "Prior", "Time", "effpopmean", "region_name")

options(repr.plot.width = 10, repr.plot.height = 7)
ggplot(data = peffpop) +
  scale_x_continuous("Time", expand = c(0,0)) +
  scale_y_continuous("Effective population size", expand = c(0,0)) +
  geom_line(aes(x = Time, y = effpop, col = Prior)) + 
  theme_bw() +
  scale_x_reverse() +
  geom_ribbon(data = peffpop, aes(x = Time, ymin= effpop025, ymax= effpop975,
                                        fill = Prior), alpha = 0.1) +
  facet_wrap(~region_name, scale= "free_y")
