#' Calculate probability of capture
#' Can be sex-, age-, and time-specific.
#'
#' @param q Catchability
#' @param Effort Total effort at time step
#' @return Probability of capture

calc_u <- function(q, Effort) {
  return(q * Effort)
}
