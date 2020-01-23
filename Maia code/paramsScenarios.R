## CENTRAL PARAMS FILE FOR KONA CRAB DEMOGRAPHY
## delete env
rm(list = ls())
## FIXED POINT ESTIMATES
## GROWTH ##
Linf = c(115,120)[1] ## Australia, Dichmont & Brown (2010)
K = c(0.26,0.45)[1] ## Australia, Dichmont & Brown (2010)
tZero = -0.1 ## Australia, Dichmont & Brown (2010)
# Linf = 121.7 ## Australia, Kirkwood (2004)
# K = 0.29 ## Australia, Kirkwood (2004)
# tZero = -0.24 ## Australia, Kirkwood (2004)
## longevity parameters
# maxage = 16 ## O'Neill
# longevDraw = maxage ## just for syntax
## Growth Increment for Logit(Pmolt); derived from vonB; Chen & Kennelly (1999)
bchen = exp(-K) - 1
achen = -1 * bchen * Linf
sdmolt = 10 ## assumes continuous variance in molt increment
## FECUNDITY ##
L50 = 60 ## Hawaii, Onizuka minimum age
## slope of length-fecundity equation (Onizuka 1972, 1000 eggs per mm)
beta = 2081.8
## fixed sex proportion of females - Onizuka
SR = 0.45
## MORTALITY ##
## hoenig (1983) "all species" values - for estimating mortality
hoenig.slope = -0.982
hoenig.int = 1.44
## HARVEST RATE POLICIES ##
## These aren't really used, have been replaced by decrement
# harvConst = 0.3
# harvMean = 0.3
# harvSD = 0.1
harvest.breaks = 10 ## number of even breaks from 0 - 0.9 to simulate harvest mortality (can't do 100%, gives -Inf
# for rVal)
hmin = 0
hmax = 0.9
h.vector = c(seq(hmin,hmax,1/harvest.breaks))
## Parameters for BISSECTION METHOD to estimate the stationary harvest, Hstat
# Define the starting h values (for function runBissect)
Uhigh <- .9 # High harvest rate to result in an eigenDom < 1
Ulow <- 0 # Low harvest rate to result in an eigenDom >1
# Define number of iterations in the Bissection method
bissectIters <- 200
# Define convergence level in the Bissection method
bissectConv <- 0.00001
## VECTOR OF TIME AT CAPTURE CUTOFFS
tc.vector = c(0,35,L50,102) ## all ages, smallest observed catch, L50, legal size
## make INITS object to save later
inits = list( 'TIMESTAMP' = as.character(date()),
              'Linf' = Linf,
              'K' = K,
              'achen [growth increment for logit p(molt)]' = achen,
              'bchen [growth increment for logit p(molt)]' = bchen,
              'sdmolt' = sdmolt,
              'tZero' = tZero,
              'L50' = L50,
              'hoenig.slope' = hoenig.slope,
              'hoenig.int' = hoenig.int,
              'beta' = beta,
              'SR [sex ratio]' = SR,
              'harvest.breaks' = harvest.breaks,
              'UHigh' = Uhigh,
              'Ulow' = Ulow,
              'bissectIters' = bissectIters,
              'bissectConv' = bissectConv,
              'tc.vector' = tc.vector)