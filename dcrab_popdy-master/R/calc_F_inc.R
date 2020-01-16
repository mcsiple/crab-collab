#' Calculate incidental fishing mortality
#' Can be sex-, age-, and time-specific.
#'
#' @param u Probability of being captured
#' @param omega Probability of being retained
#' @param sigma Probability of incidental death
#' @return Incidental fishing mortality

calc_F_inc <- function(u,omega,sigma) {
  return(u * (1-omega) * sigma)
}
