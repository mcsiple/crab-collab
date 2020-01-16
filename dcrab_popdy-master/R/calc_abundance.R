#' Calculate abundance
#' Using number of individuals. Can be sex- and/or age-specific.
#'
#' @param N_prev Abundance (number of individuals) from previous time step
#' @param M Natural mortality
#' @param F_ret Retained fishing mortality
#' @param F_inc Incidental fishing mortality
#' @return Abundance (number of individuals) in the next time step

calc_abundance <- function(N_prev,M,F_ret,F_inc) {
  N_next <- N_prev*exp(-(M + F_ret + F_inc))
  return(N_next)
}
