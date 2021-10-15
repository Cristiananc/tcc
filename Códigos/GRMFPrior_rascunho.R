#Simulating with a fixed value for tau
library(phylodyn)
library(ggplot2)
library(latex2exp)
library(Matrix)
library(cowplot)
library(INLA)

tau <- 1

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

#
x = inla.qsample(1, Q=tau*Q1D)

#Sampling from IGMRF
ev <- eigen(tau*Q1D)
y <- rep(0,k-1)
for(i in 1:k-1){
  y[i] <- rnorm(1, mean = 0, sd = sqrt(1/ev$values[i+1]))
}
vectors <- ev$vectors[,-1] 
a <- apply(y*vectors,1,sum,na.rm = TRUE)
d <- data.frame(a)
d['x'] <- midpoint_distance_x

#For sampling from GMRF priors
L <- chol(tau*Q1D)
v <- matrix(rnorm(k*(k-1),mean=0,sd=1), k, k-1) 
samps   <- solve(t(L),v)
plot(midpoint_distance_x, samps[,1], type = 'l')

data_plot <- data.frame(midpoint_distance_x)

#Wiener Process
X_simulation <- function(n_sim, tau){
  for(j in 1:n_sim){
  #We start with W(0) = 0
  #Preallocate arrays
  dX <- rep(0, k-1)
  X <- rep(0, k)
  var <- 0
  
  for (i in 2:k){
    var <- midpoint_distance[i-1]/tau
    dX[i-1] <- rnorm(1, mean = 0, sd = sqrt(var))
    X[i] <- X[i-1] + dX[i-1]
  }
  data_plot$y = X
  names(data_plot)[j+1] <- paste("X", j, sep="")
  }
  return(data_plot)
}

#Plotting 5 realizations of the process
data_plot <- X_simulation(10, tau)
g <- ggplot(data = data_plot, aes(x = midpoint_distance_x)) +
 geom_line(aes(y = X1)) +
  geom_line(aes(y = X2)) +
  geom_line(aes(y = X3)) +
  geom_line(aes(y = X4)) +
  geom_line(aes(y = X5)) +
  geom_line(aes(y = X6)) +
  geom_line(aes(y = X7)) +
  geom_line(aes(y = X8)) +
  geom_line(aes(y = X9)) +
  geom_line(aes(y = X10)) +
  theme_bw() +
  ggtitle(TeX("$\\tau = 0.1")) +
  xlab("") + ylab("")

ggsave("wiener_process_01.pdf",
       plot = g)

data_plot1 <- data.frame(midpoint_distance_x)

X_simulation <- function(n_sim, tau){
  for(j in 1:n_sim){
    #We start with W(0) = 0
    #Preallocate arrays
    dX <- rep(0, k-1)
    X <- rep(0, k)
    var <- 0
    
    for (i in 2:k){
      var <- midpoint_distance[i-1]/tau
      dX[i-1] <- rnorm(1, mean = 0, sd = sqrt(var))
      X[i] <- X[i-1] + dX[i-1]
    }
    data_plot1$y = X
    names(data_plot1)[j+1] <- paste("X", j, sep="")
  }
  return(data_plot1)
}

#Plotting for a different value of tau
tau1 <- 1
data_plot1 <- X_simulation(10, tau1)
g1 <- ggplot(data = data_plot1, aes(x = midpoint_distance_x)) +
  geom_line(aes(y = X1)) +
  geom_line(aes(y = X2)) +
  geom_line(aes(y = X3)) +
  geom_line(aes(y = X4)) +
  geom_line(aes(y = X5)) +
  geom_line(aes(y = X6)) +
  geom_line(aes(y = X7)) +
  geom_line(aes(y = X8)) +
  geom_line(aes(y = X9)) +
  geom_line(aes(y = X10)) +
  theme_bw() +
  ggtitle(TeX("$\\tau = 1")) +
  xlab("") + ylab("")

g1

ggsave("wiener_process_1.pdf",
       plot = g1)
g
