#' Calculate week duration between two dates
#'
#' @param curdt Date to calculate to in string format, e.g., "August 27 2019"
#' @param stdt Date to calculate from in string format. Default is "November 1, 1999"
#' @return Difference in two dates in weeks

calc_week <- function(curdt,stdt="November 1, 1999") {
  lst <- lubridate::mdy(stdt)
  lcr <- lubridate::mdy(curdt)
  return(as.double(difftime(lcr,lst,units="weeks")))
}
