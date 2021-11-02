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
