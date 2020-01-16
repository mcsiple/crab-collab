#' Calculate Beverton-Holt recruitment
#' Assumes females will only mate with males older than them. Assumes oldest females will reproduce with oldest males.
#' Follows parametrisation from Brooks and Powers 2007
#'
#' @param h steepness
#' @param R0 virgin recruitment
#' @param phi0 virgin level of eggs per recruit
#' @param E number of eggs; the number at time 0 of stage 1

calc_BevHolt <- function(h, R0, phi0, E) {
  R <- (4*h*R0*E) / ((phi0 * R0 * (1-h)) + ((5*h + 1)*E))
  return(R)
}
