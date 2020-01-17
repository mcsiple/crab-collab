# run-model
library(here)

# Read in all the functions from Dungeness crab model
setwd(here("dcrab_popdy-master","R"))
files.sources <- list.files()
sapply(files.sources, source)
# Now set the working dir back to the project dir
setwd(here())


# Now get to work!
Ninit <- 1000
R0 <- 1e6
nyrs <- 20
SurvJ <- 0.3/12 # juvenile survival (annual survival = 0.3, Grubert ete al 2019)
SurvA <- 0.6/12 # adult survival (I made this one up, might decrease)
InitDepl <- 1 # starting pop size relative to K
K <- 2000 
nages <- 10 # "max" age, treated here as plus group
PropFem <- 0.5 # Proportion of crabs that are female at t=1 - can be a vector or single value
nmonths <- 12 # just for testing


# One vector per sex
Neq <- Ninit <- FemInit <- MaleInit <- vector(length = nages+1)

# Initial nums at age needs to be a stable age distribution
# Simplest case: pop starts at equilibrium
Neq[1] <- 1 # Age 0 (new recruits)
Neq[2] <- SurvJ
for(a in 3:nages){
  Neq[a] <- Neq[a-1] * SurvA
}
Neq[nages+1] <- (SurvJ*SurvA^(nages-1))/(1-SurvA) # plus group
R0 <- K/sum(Neq[2:(nages+1)])    # numerical soln

if(InitDepl == 1){
  Ninit <- Neq
}else{
  Ninit <- Nnoneq
}

FemInit <- Ninit*PropFem
MaleInit <- Ninit*(1-PropFem)

Fem[1,] <- FemInit
Male[1,] <- MaleInit

for(t in 2:nmonths){
  Fem[t,2:(nages-1)] <- calc_abundance(Nprev = N[y-t,],
                         M = M,
                         F_ret = F_ret,
                         F_inc = F_inc)
}

calc_abundance(N_prev = N0,
               M = 0.2,
               F_ret = 0.5, # retained fishing mortality
               F_inc = 0.01) # incidental fishing mortality
