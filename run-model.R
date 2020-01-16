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

N[1] <- Ninit
# N needs to be a stable age distribution

for(y in 2:nyrs){
  N[y] <- calc_abundance(N_prev = N[y-1],
                         M = M,
                         F_ret = F_ret,
                         F_inc = F_inc)
}

calc_abundance(N_prev = N0,
               M = 0.2,
               F_ret = 0.5, # retained fishing mortality
               F_inc = 0.01) # incidental fishing mortality
