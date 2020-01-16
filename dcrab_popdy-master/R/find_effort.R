find_effort <- function(effort) {
  upar <- calc_u(qpar,effort)
  F_ret <- calc_F_ret(u, omega)
  predcatch <- calc_catch(N,M,F_ret,F_inc)

  NLL <- (predcatch-truecatch)^2
}
