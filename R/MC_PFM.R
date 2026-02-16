MC_PFM <- function(T, Cin, isotope) {
  T_sup <- ceiling(T)
  sim <- data.table::shift(Cin, T_sup, fill = 0)
  if (isotope ==  "3h") sim <- sim * exp(-0.0562619464740215 * T)
  return(sim)
}
