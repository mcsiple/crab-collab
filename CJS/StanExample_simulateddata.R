library(here)
setwd(here("CJS"))
# Set up the true parameter values
a <- c(.8, 1)
b <- c(2, .1)
sigma <- .2
# Simulate data
x <- (1:1000)/100
N <- length(x)
ypred <- a[1]*exp(-b[1]*x) + a[2]*exp(-b[2]*x)
y <- ypred*exp(rnorm(N, 0, sigma))
# Fit the model
library("rstan")

#MCS: rec'd for running stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit <- stan("exponentials.stan", data=list(N=N, x=x, y=y), iter=1000, chains=4)
print(fit, pars=c("a", "b", "sigma"))
