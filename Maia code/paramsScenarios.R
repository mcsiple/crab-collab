## CENTRAL PARAMS FILE FOR KONA CRAB DEMOGRAPHY
## FIXED POINT ESTIMATES


# GROWTH ------------------------------------------------------------------
Linf = 310  # from S. serrata, Moknes et al. 2014
K = c(0.57,1.38)[1] # from S. serrata, Moknes et al. 2014 (first value)
                    ## S. serrata, Indonesia, La Sara (2010): 1.38
tZero = -0.019       ## from S. serrata, Moknes et al. 2014

## longevity parameters
# longevDraw = maxage ## just for syntax

## Growth Increment for Logit(Pmolt); derived from vonB; Chen & Kennelly (1999)
bchen = exp(-K) - 1
achen = -1 * bchen * Linf
sdmolt = 10 ## assumes continuous variance in molt increment



# FECUNDITY ---------------------------------------------------------------
size.at.maturity = 100 # 91-100 mm (females, Prasad & Neelakantan 1989)
L50 = (91+100)/2 ## Midpoint of size at maturity from Prasad & Neelakantan (1989)

## slope of length-fecundity equation 
## For Samoan crab, 15.55 x 100 eggs increase per mm of carapace width - Sarower et al. (2013), Bangladesh - check w Maia ( 1000 eggs per mm is Kona crab value beta = 2081.8; from Onizuka 1972)
beta = 15550 #1500 
## fixed sex proportion of females - fixed by me (Megsie) because we don't know the sex ratio...
SR = 0.5


# MORTALITY ---------------------------------------------------------------
## hoenig (1983) "all species" values - for estimating mortality
hoenig.slope = -0.982
hoenig.int = 1.44


# HARVEST RATE POLICIES ---------------------------------------------------
harvest.breaks = 10 ## number of even breaks from 0 - 0.9 to simulate harvest mortality (can't do 100%, gives -Inf
# for rVal)
hmin = 0
hmax = 0.9 # max harvest rate
h.vector = c(seq(hmin,hmax,1/harvest.breaks))
## Parameters for BISSECTION METHOD to estimate the stationary harvest, Hstat
# Define the starting h values (for function runBissect)
Uhigh = 0.9 # High harvest rate to result in an eigenDom < 1
Ulow = 0 # Low harvest rate to result in an eigenDom >1

# Define number of iterations in the Bissection method
bissectIters <- 200
# Define convergence level in the Bissection method
bissectConv <- 0.00001

# VECTOR OF SIZE AT CAPTURE CUTOFFS ---------------------------------------
tc.vector = c(0,76,L50,152) ## all ages, smallest observed catch (this study), L50, legal size for Samoan crab from HI state (6 in = 152.4 mm)
HeeiaSel <- c(0, 0, 0, 0, 0.02, 0.14, 0.58, 1, 1) # discretized version of logistic sel curve

# make INITS object to save later -----------------------------------------
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