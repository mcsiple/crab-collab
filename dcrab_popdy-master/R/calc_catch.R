#' Calculate catch
#' Based on Baranov catch equation. Not totally sure what units are; not based on abundance?
#' Can be sex- and/or age-specific.
#'
#' @param N Abundance of exploitable population
#' @param M Natural mortality
#' @param F_ret Retained fishing mortality
#' @param F_inc Incidental fishing mortality
#' @return Catch for that time step. Not totally sure what the units are.

calc_catch <- function(N,M,F_ret,F_inc) {
  catch_t <- N * M * F_ret * (1 - exp(-(M + F_ret + F_inc))) / (M + F_ret + F_inc)
  return(sum(catch_t))
}
