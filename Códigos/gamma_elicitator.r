#### Elicit a gamma distribution such that Pr(X > S) = p
### Note that Pr(1/sqrt(tau) > S) = p can be written as
#### Pr(tau < 1/S^2) = p
#### here, tau ~ Gamma(shape = a, scale = b)

loss <- function(par, S, p){
  a <- exp(par[1])
  b <- exp(par[2])
  abs(
    pgamma(q = 1/S^2, shape = a, scale = b) - p
  )
}

get_gamma_pars <- function(S, p){
  opt <- optim(par = c(0, 0), fn = loss,
               S = S, p = p)
  out <- exp(opt$par)
  return(out)
}

solnl(X = NULL, objfun = NULL)

SS <- 1
pp <- .1

the.pars <- get_gamma_pars(S = SS, p = pp)

curve(dgamma(x, shape = the.pars[1],
             scale = the.pars[2]), 0, 10*SS,
      ylab = "Density")

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


#Fazendo a otimização para o caso onde queremos que E[tau] = 1
library(nloptr)

#Restrição e função objetivo não lineares
#Restrição igualdade
local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
opts <- list( "algorithm"= "NLOPT_GN_ISRES",
              "xtol_rel"= 1.0e-15,
              "maxeval"= 1600000,
              "local_opts" = local_opts,
              "print_level" = 0 )

res <- nloptr ( x0 = c(0.1, 0.1),
               eval_f = loss,
               lb = c(0, 0),
               ub = c(10, 10),
               eval_g_eq = eq,
               opts = opts, S = S, p = p
)

#Perguntar sobre o exp(out)