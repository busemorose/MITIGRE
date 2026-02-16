MC_DM <- function(T, DP, isotope, max_t = 1000) {
  t <- 1:max_t

  a <- (t - 0.5) / T
  C <- (1 / ((t - 0.5) * sqrt(4 * pi * DP * a))) * exp(-((1 - a) ^ 2) / (4 * DP * a))

  MC <- C / sum(C)
  if (isotope ==  "3h") MC <- MC * exp(-0.0562619464740215 * t)
  return(MC)
}
