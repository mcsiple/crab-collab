# crab_CJS1
library(here)
library("rstan")

# Load data
y <- as.matrix(read.table("crabmat.txt"))
colnames(y) <- NULL
nind <- nrow(y)
n_occasions <- ncol(y)

setwd(here("CJS"))
set.seed(123)

## Initial values
inits <- function() list(mean_phi = runif(1, 0, 1),
                         mean_p = runif(1, 0, 1))

## Parameters monitored
params <- c("mean_phi", "mean_p")

## MCMC settings
ni <- 2000 
nt <- 1
nb <- 1000
nc <- 4

# Fit the model
#MCS: rec'd for running stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit <- stan("cjs_ex.stan",
            data=list(nind=nind, n_occasions=n_occasions, y=y), 
            init=inits,
            pars=params,
            chains = nc, iter = ni, warmup = nb, thin = nt)
print(fit, pars=c("mean_phi", "mean_p"))
print(fit, digits=3)

#mean_phi is the mean survival
#mean_p is the recapture probability
matrix_of_draws <- as.matrix(fit)
shinystan::launch_shinystan(fit) # use shinystan to look at convergence, etc.
