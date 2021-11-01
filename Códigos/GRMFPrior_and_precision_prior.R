#Simulating with a fixed value for tau
library(phylodyn)
library(ggplot2)
library(latex2exp)
library(INLA)
library(gridExtra)
library(VaRES)

#Parameters for the precision prior
gamma2 <- c(2.544882, 1.203437)
pc_prior <- c(0.5, 2.302585)

set.seed(9)
trajectory = exp_traj
samp_times = 0
n = 100

#Simulating a genealogy based on the sample
gene = coalsim(samp_times = samp_times, n_sampled = n,
               traj = trajectory, lower_bound = 1/10)

k <- n - 1

midpoint_distance <- (gene$intercoal_times[1:k-1] + gene$intercoal_times[2:k])/2

midpoint_distance_x <- rep(0, k)
for(i in 2:k){
  midpoint_distance_x[i] = midpoint_distance[i-1] + midpoint_distance_x[i-1]
}

#Draw a random value for tau until it's different of zero
random_tau <- function(tau_prior) {
  success <- FALSE
  
  if (tau_prior == 1){
    x <- rgamma(1, shape = 0.001, rate = 0.001)
    }
  if (tau_prior == 2){
      x <- rgamma(1, shape = gamma2[1], scale = gamma2[2])
    }    
  if(tau_prior == 3){
      x <- dgumbel2(1, a = pc_prior[1], b = pc_prior[2])
    }    
  return(x)
}

#Run the simulation a given tau and gamma0 priors
X_simulation <- function(gamma0_prior, tau_prior){
  u0 <- rcauchy(1, location = 0, scale = 1)
  tau <- random_tau(tau_prior)
  prior <- list(0, rnorm(1, 0, 1/sqrt(tau)), runif(1, -10, 10), rnorm(1, 0, 10^(6)), rnorm(1, u0, 10^(6)))
  
  #Preallocate arrays
  dX <- rep(0, k-1)
  X <- rep(0, k)
  var <- 0
  X[1] <- prior[[gamma0_prior]]

  for (i in 2:k){
    var <- midpoint_distance[i-1]/tau
    if (var)
    dX[i-1] <- rnorm(1, mean = 0, sd = sqrt(var))
    #ERROR nans produced when var = Inf
    if (var == Inf){
      var <- .Machine$double.xmax
    }
    X[i] <- X[i-1] + dX[i-1]
  }
  return(X)
}

#Simulation for each type of precision prior and gamma0 prior
nrep <- 1000
df_simulations <- data.frame()
prior_names <- c('0', 'N(0, 1/sqrt(tau))', 'Uniform',
                 'N(0, 10^6)', 
                 'N(u, 10^6) u~Cauchy(0, 1)')
prior_tau_names <- c('Gamma usual', 'Gamma2', 'PC prior')

simulation_for_each_prior <- function (nrep){
  for(l in 1:5){
    for (j in 1:3){
      simulations <- do.call(rbind, lapply(1:nrep, function(i){
        data.frame(x_value = midpoint_distance_x, 
                   value = X_simulation(gamma0_prior = l, tau_prior = j), 
                   replicate = i, 
                   group = c(1:99),
                   prior = prior_names[l],
                   tau = prior_tau_names[j])
      }))
      df_simulations <- rbind(df_simulations, simulations)
    }
  }
  return(df_simulations)
}

#Calling simulation function for each prior
df_simulations <- simulation_for_each_prior(nrep = nrep)

#Plotting
require(gridExtra)
options(repr.plot.width = 10, repr.plot.height = 7)

#Escalas individuais no eixo y
ggplot(data = df_simulations, aes(x = x_value, y = value, 
                                  group = replicate)) +
  scale_x_continuous("Times", expand = c(0,0)) +
  scale_y_continuous("", expand = c(0,0)) +
  geom_line(alpha = .5) + 
  guides(col = "none") +
  theme_bw() +
  #ggtitle(TeX("$\\tau =$")) +
  facet_grid(prior ~ tau, scales = "free_y")

#Pseudolog scale
ggplot(data = df_simulations, aes(x = x_value, y = value, 
                                  group = replicate)) +
  scale_x_continuous("Times", expand = c(0,0)) +
  scale_y_continuous("", expand = c(0,0)) +
  geom_line(alpha = .5) + 
  guides(col = "none") +
  theme_bw() +
  #ggtitle(TeX("$\\tau =$")) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  facet_grid(prior ~ tau)

#Removing the gamma(shape = 0.001, rate = 0.001 )
df_simulations_without_gamma_usual <- df_simulations[df_simulations$tau != "Gamma usual", ]

#Escalas individuais no eixo y
ggplot(data = df_simulations_without_gamma_usual, aes(x = x_value, y = value, 
                                  group = replicate)) +
  scale_x_continuous("Times", expand = c(0,0)) +
  scale_y_continuous("", expand = c(0,0)) +
  geom_line(alpha = .5) + 
  guides(col = "none") +
  theme_bw() +
  #ggtitle(TeX("$\\tau =$")) +
  facet_grid(prior ~ tau, scales = "free_y")
