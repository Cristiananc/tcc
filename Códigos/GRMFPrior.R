#Simulating with a fixed value for tau
library(phylodyn)
library(ggplot2)
library(latex2exp)
library(Matrix)
library(cowplot)
library(INLA)

tau1 <- 1

trajectory = exp_traj
samp_times = 0
n = 100

#Simulating a genealogy based on the sample
gene = coalsim(samp_times = samp_times, n_sampled = n,
               traj = trajectory, lower_bound = 1/10)
k <- n - 1

midpoint_distance <- (gene$intercoal_times[1:k-1] + gene$intercoal_times[2:k])/2

#Precision Matrix
gene_intercoal1 <- replicate(k-1, 1)/midpoint_distance

midpoint_distance_x <- rep(0, k)
for(i in 2:k){
  midpoint_distance_x[i] = midpoint_distance[i-1] + midpoint_distance_x[i-1]
}

diag <- c(gene_intercoal1[1], gene_intercoal1[1:(k-2)] + gene_intercoal1[2:(k-1)], gene_intercoal1[k-1])
Mdiag <- Diagonal(x = diag)
Q1D <- Mdiag - sparseMatrix(i= seq(from = 2, to = k, by = 1), j = seq(from = 1, to = (k-1), by = 1), x= gene_intercoal1[1:k-1], dims=c(k, k), symmetric=TRUE)


data_plot <- data.frame(midpoint_distance_x)


X_simulation <- function(tau){
  #theta0 <- runif(1, -10, 10) #diffuse or non-informative prior
  theta1 <- rnorm(1, 0, 1/sqrt(tau1))
  #Preallocate arrays
  dX <- rep(0, k-1)
  X <- rep(0, k)
  var <- 0
  X[1] <- theta1
  
  for (i in 2:k){
    var <- midpoint_distance[i-1]/tau
    dX[i-1] <- rnorm(1, mean = 0, sd = sqrt(var))
    X[i] <- X[i-1] + dX[i-1]
  }
  return(X)
}

#Simulation
nrep <- 1000
simulations <- do.call(rbind, lapply(1:nrep, function(i){
  data.frame(x_value = midpoint_distance_x, value = X_simulation(tau = tau1), replicate = i, 
             group = c(1:99))
}))

#Plotting
#for_plot$replicate <- as.factor(for_plot$replicate)
options(repr.plot.width = 10, repr.plot.height = 7)
plot_simulation <- ggplot(data = simulations, aes(x = x_value, y = value, group = replicate)) +
  scale_x_continuous("Times", expand = c(0,0)) +
  scale_y_continuous("", expand = c(0,0)) +
  geom_line(alpha = .5) + 
  guides(col = "none") +
  theme_bw() +
  ggtitle(TeX("$\\tau = 1")) +
  NULL

plot_simulation

#Mean
mean <- aggregate(.~group,data=simulations,FUN=sum)
plot(midpoint_distance_x, mean$value/99, type = 'l')
