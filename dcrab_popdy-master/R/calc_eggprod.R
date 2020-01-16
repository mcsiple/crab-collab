#' Calculate egg production
#' Assumes females will only mate with males older than them. Assumes oldest females will reproduce with oldest males. Sums weeks into years.
#'
#' @param N_female Vector of female abundance by age
#' @param N_male Vector of male abundance by age
#' @param fec Vector of fecundity-at-age
#' @param delta Vector of proportion of mature females at age a
#' @param eta Half saturation constant
#' @param byweek Logical; default is T; combines numbers-at-age from weeks to years
#' @return Incidental fishing mortality

calc_eggprod <- function(N_female, N_male, fec, delta, eta, byweek=T) {
	if(byweek) {
		N_female_y <- N_male_y <- c()
		for(a in 1:(length(N_female)/52)) {
			N_female_y[a] <- sum(N_female[(((a-1)*52)+1):(a*52)])
			N_male_y[a] <- sum(N_male[(((a-1)*52)+1):(a*52)])
		}
	}

	else {
		N_female_y <- N_female
		N_male_y <- N_male
	}
	N_oldermales <- rev(cumsum(rev(N_male_y[-1])))
	N_oldermales <- c(N_oldermales,rev(N_oldermales)[1])
	rho <- N_oldermales / N_female_y
	phi <- rho/(eta+rho)
	sum(N_female * fec * delta * phi)
}
