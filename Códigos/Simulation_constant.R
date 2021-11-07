#Simulation with constant trajectory
library(ggplot2)
library(phylodyn)

date_min <- 1990
date_max <- 2021

const10 <- function(x){
  return(rep(10, length(x)))
}

nrep <- 500

#Sample times - uniform distribution
#20 taxa
stimes20 <- sort(runif(20, date_min, date_max))

#Genealogies
#(taxa — pl.) Any named group of organisms (e.g., the reptiles, Felidae, beetles, Homo sapiens).

genealogies_20taxa <- lapply(1:nrep, function(i){
  phylodyn::coalsim(samp_times = stimes20,
                    n_sampled = rep(1,20),
                    traj = const10,
                    lower_bound = 1)
})


#200 taxa
stimes200 <- sort(runif(200, date_min, date_max))

genealogies_200taxa <- lapply(1:nrep, function(i){
  phylodyn::coalsim(samp_times = stimes200,
                    n_sampled = rep(1,200),
                    traj = const10,
                    lower_bound = 1)
})

#Phylodynamic reconstruction of effective population size
#20 taxa
#Gamma flat
system.time(
  reconstructions_20taxa_gamma <- lapply(genealogies_20taxa, function(x)
    BNPR(data = x, lengthout =100, prec_alpha = 0.001, prec_beta = 0.001))
)

#Gamma 2
par_ab <- c(2.544882, 1.203437)
system.time(
  reconstructions_20taxa_gamma2 <- lapply(genealogies_20taxa, function(x)
    BNPR(data = x, lengthout =100, prec_alpha = par_ab[1], prec_beta = 1/par_ab[2]))
)

#PC Prior
system.time(
  reconstructions_20taxa_gumbel <- lapply(genealogies_20taxa, function(x)
    BNPR(data = x, lengthout = 100, 
         pc_prior = TRUE))
)

#Phylodynamic reconstruction of effective population size - 200 taxa
system.time(
  reconstructions_200taxa_gamma <- lapply(genealogies_200taxa, function(x)
    BNPR(data = x, lengthout = 100))
)

#Gamma 2
par_ab <- c(2.544882, 1.203437)
system.time(
  reconstructions_200taxa_gamma2 <- lapply(genealogies_200taxa, function(x)
    BNPR(data = x, lengthout =100, prec_alpha = par_ab[1], prec_beta = 1/par_ab[2]))
)

#PC prior
system.time(
  reconstructions_200taxa_gumbel <- lapply(genealogies_200taxa, function(x)
    BNPR(data = x, lengthout = 100, 
         pc_prior = TRUE))
)

