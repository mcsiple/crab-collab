#' Calculate retained fishing mortality
#' Can be sex- and/or age-specific.
#'
#' @param u Probability of being captured
#' @param omega Probability of being retained
#' @return Retained fishing mortality

calc_F_ret <- function(u,omega) {
  return(u * omega)
}
