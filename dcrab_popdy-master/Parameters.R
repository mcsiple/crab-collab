# devtools::install()
library(popdy)

# Assume starting from 1999 season
# Based on CDFW data
Dates <- c("October 1, 1999", "April 15, 2019")
Nages <- 10 # Years, which will be converted to weeks later
Nsex <- 2
Nblocks <- 200
Nweeks <- ceiling(calc_week(Dates[2],Dates[1]))
# Number of weeks from Oct 1
Tspawn <- 1 # "Most apparently spawn during October and early November" (Wild & Tasto 1983). As fishing season is never going to start before Nov anyway, we can just assume spawning starts in Oct 1?
FishingWeeks <- ceiling(calc_week("November 15, 1999",Dates[1])):round(calc_week("June 30, 2000",Dates[1])) + 1
# Trecruit <- # Don't think this matters

# Fishery and population parameters
# Table 1 from Froehlich et al. 2016

# Age-specific Natural Mortality
NatM <- rep(NA, Nages*52)
NatM[1:(2*52)] <- 0.8/52
NatM[(2*52+1):(Nages*52)] <- 0.3/52

# Age-specific Fecundity
Fec <- rep(NA, Nages*52)
Fec[1:(1*52)] <- 0
Fec[(1*52+1):(2*52)] <- 1E6
Fec[(2*52+1):(5*52)] <- 2E6
Fec[(5*52+1):(Nages*52)] <- 1E6

# Age-specific Proportion of Mature Females
Prop_MatFem <- rep(NA, Nages*52)
Prop_MatFem[1:(1*52)] <- 0
Prop_MatFem[(1*52+1):(2*52)] <- .2
Prop_MatFem[(2*52+1):(Nages*52)] <- 1

# Effective Reproductive Ratio
#EffProdRatio <- range(0,1) # Given as a range

# Half-saturation constant for "effective reproductive ratio"
Eta <- .1

# Recruitment Parameters
EggsPerRecruit <- 630222 # Calculated based on M, Fec, and Prop_Mature
Rzero <- 7E6 # ***Specific to Hood Canal?***
Steepness <- .65 # Range given between 0.5 and 0.8

# Age- and sex-specific Probability of being Retained
# Currently assumes no illegal catch
# Can also be modelled as a function of weekly catch rate of legal males
Omega <- matrix(ncol=Nages*52,nrow=Nsex) #By sex, row1=male, row2=female
Omega[1,1:(2*52)] <- 0
Omega[1,(2*52+1):(3*52)] <- 0 # Range given between 0 and 0.8
Omega[1,(3*52+1):(Nages*52)] <- 1
Omega[2,1:(3*52)] <- 0
Omega[2,(3*52+1):(Nages*52)] <- 0 # Range given between 0 and 0.8

# Age-specific Probability of Incidental Death
Sigma <- rep(NA, Nages*52)
Sigma[1:(2*52)] <- 0
Sigma[(2*52+1):(Nages*52)] <- .5 # Range given between 0 and 0.8


# Age- and sex-specific Catchability (in the absence of hypoxia)
qpar <- matrix(ncol=Nages*52,nrow=Nsex) # By sex, row1=male, row2=female
qpar[1,1:(2*52)] <- 0
qpar[1,(2*52+1):(Nages*52)] <- .00035
qpar[2,1:(3*52)] <- 0
qpar[2,(3*52+1):(Nages*52)] <- .00035

# Does effort not change over time??
# Figure out how to divide this over the season based on historic estimates
# Bulk of catch made in first six weeks over the holiday season
# Total weekly crabbing effort for CA coast = 173900; out-of-state permits have 21425 traps?
# Needs to be spatial and temporal, what to do in between seasons? Set to zero
# Include recreational effort?
Effort <- rep(c(rep(0,6),rep(17390,6),rep(.1*17390,(length(FishingWeeks)-6)),rep(0,12)),20)
# q_max <-.005 # Do we need this? Don't think so.
# k <- .00015

# Initialising
# Need to figure out which of these have spatial components
Ncrabs <- array(dim=c(Nsex,Nages*52,Nweeks))
Ncrabs[,,1] <- 7E6

Catch <- c()

# Tune effort to data?

for(t in 1:(Nweeks-1)) {

	upar <- calc_u(qpar,Effort[t])
	# F_at_age for a single week (based on effort specified above)
	F_ret <- calc_F_ret(upar,Omega)
	F_inc <- calc_F_inc(upar,Omega,Sigma)
	
	Ncrabs[,,(t+1)] <- calc_abundance(Ncrabs[,,t],M=NatM,F_ret=F_ret,F_inc=F_inc)
	Catch[t] <- calc_catch(N=Ncrabs[,,(t+1)],M=NatM,F_ret=F_ret,F_inc=F_inc) # In numbers of crab

	# Reproduction and recruitment occur here, after catch has happened for that week
	if((t %% 52) == Tspawn) {
		EggProd <- calc_eggprod(N_female=Ncrabs[2,,(t+1)],N_male=Ncrabs[1,,(t+1)],
							fec=Fec, delta=Prop_MatFem, eta=Eta)
		Ncrabs[,1,(t+1)] <- calc_BevHolt(h=Steepness, R0=Rzero, phi0=EggsPerRecruit, E=EggProd)/2 # Assume .5 probability for each sex
	}
}


