#### Elicit a gamma distribution such that Pr(X > S) = p
### Note that Pr(1/sqrt(tau) > S) = p can be written as
#### Pr(tau < 1/S^2) = p
#### here, tau ~ Gamma(shape = a, scale = b)
library(ggplot2)
library(nloptr)
library(VaRES) #Gumbel2 distribution
library(latex2exp)

set.seed(29)

#Considerando 
loss <- function(par, S, p){
  a <- exp(par[1])
  abs(
    pgamma(q = 1/S^2, shape = a, scale = 1/a) - p
  )
}

get_gamma_pars <- function(S, p){
  opt <- optim(par = c(0), fn = loss,
               S = S, p = p, method = "Brent", lower = 0, upper = 100)
  out <- exp(opt$par)
  return(out)
}

SS <- 1
pp <- .1

the.pars <- get_gamma_pars(S = SS, p = pp)
pc_prior <- c(0.5, 2.302585)
matching_gamma_par <- c(2.544882, 1.203437)

#Plots
curve(dgamma(x, shape = the.pars[1],
             scale = the.pars[2]), 0, 10*SS,
      ylab = "Density")

#Density plots
ggplot() +
  xlim(0, 10) +
  geom_function(aes(colour = "Gamma Flat"), fun = dgamma, args = list(shape = 0.001, rate = 0.001)) +
  geom_function(aes(colour = "PC Prior"), fun = dgumbel2, args = list(a = pc_prior[1], b = pc_prior[2])) +
  geom_function(aes(colour = "Matching Gamma"), fun = dgamma, args = list(shape = matching_gamma_par[1], scale = matching_gamma_par[2])) +
  theme_bw() +
  theme(text = element_text(size=15)) +
  xlab(TeX("$\\tau$")) +
  ylab("Density") +
  guides(colour=guide_legend(title="Priors"))

precisions <- rgamma(n = 1E6,
            shape = the.pars[1],
            scale = the.pars[2]
            )

sds <- 1/sqrt(precisions)

mean(sds > SS)

hist(precisions)

the.pars

#Pr(tau > 10^-5) = 10^-4 
aux1 <- 1/sqrt(10^(-5)) 
aux2 <- 10^(-5)
get_gamma_pars(S = aux1, p = aux2)

#Fazendo a otimização para o caso onde queremos que E[tau] = 1 usando nloptr
eq <- function(par, S, p) {
  par[1]*par[2] - 1
}

#Restrição e função objetivo não lineares
#Restrição igualdade
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "local_opts" = local_opts,
              "maxeval"= 1600000,
              "print_level" = 0 )

res <- nloptr ( x0 = c(0.1, 0.1),
               eval_f = loss,
               lb = c(0, 0),
               ub = c(100, 100),
               eval_g_eq = eq,
               opts = opts, S = SS, p = pp
)

solution <- res$solution

#Gamma Flat(shape = 0.001, rate = 0.001)
#Pr(1/sqrt(tau)>1)
pgamma(q = 1, shape = 0.001, rate = 0.001)

