ui_files <- c("explore.R", "notes.R", "ui_helpers.R")
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) {
  if (!basename(f) %in% ui_files) source(f, encoding = "UTF-8")
}

ok <- function(cond, msg) {
  cat(if (isTRUE(cond)) "[OK]   " else "[ÉCHEC]", msg, "\n")
  invisible(isTRUE(cond))
}
results <- logical(0)

prep <- list(dt = 1 / 12, o = 1, uz = 0, kmax = 1e6, lam_year = 0, isotope = "18o", tol = 1e-7, max_t = 1000)
for (cs in list(list("EMM", 300, NA), list("EPM", 300, 1.5), list("PEM", 300, 1.1),
                list("DM", 300, 0.5), list("DM", 300, 0.01))) {
  k <- tlpm_kernel(cs[[1]], cs[[2]], cs[[3]], prep, decay = FALSE)
  results <- c(results, ok(abs(sum(k$w) - 1) < 1e-6 && abs(kernel_mean_age(k) - 25) < 1e-3,
                           sprintf("TracerLPM %s : masse %.7f, âge moyen %.4f a (attendu 25)",
                                   cs[[1]], sum(k$w), kernel_mean_age(k))))
}

for (n in c(1.1, 2, 5)) {
  a <- tlpm_kernel_w("PEM", 30, n, log(2) / 12.32, 1 / 12, 1, 0, 3000)$w
  b <- tlpm_kernel_w("EPM", 30, 1 + log(n), log(2) / 12.32, 1 / 12, 1, 0, 3000)$w
  results <- c(results, ok(max(abs(a - b)) < 1e-12, sprintf("PEM(n = %.1f) == EPM(1 + ln n)", n)))
}

set.seed(3)
x <- stats::rnorm(600, -7, 1)
d <- data.frame(obs = x, Cin = x, date = NA_real_)
p <- make_prep(d, "18o", 1 / 12, "0", 0L, "1", max_t = 5000, tol = 1e-12, kmax = 1e6, uz = 0)
for (cs in list(list("EMM", 30, NA), list("DM", 30, 0.3), list("PEM", 30, 1.5))) {
  s1 <- mitigre_component(cs[[1]], cs[[2]], cs[[3]], p)
  s2 <- claude_component(cs[[1]], cs[[2]], cs[[3]], p)
  results <- c(results, ok(max(abs(s1 - s2)) < 1e-8,
                           sprintf("MITIGRE == Forward (%s, 18O) : max |écart| = %.2e", cs[[1]], max(abs(s1 - s2)))))
}

vba_EMM <- function(dates, Cin, Tau, S, lam) {
  TauRes <- Tau
  MinAge <- S
  TI <- if (Tau <= 100) 1 / 12 else abs(dates[2] - dates[1])
  MaxAge <- trunc(MinAge) + trunc((MinAge - trunc(MinAge)) / TI) * TI
  End <- S
  MinAge <- S - MinAge
  MaxAge <- S - MaxAge
  j <- length(dates)
  while (j != 0 && dates[j] >= End - MinAge) j <- j - 1
  j <- max(j, 1)
  mult <- 1 / TauRes * (1 / (1 / TauRes + lam))
  h1 <- exp(-MinAge * (1 / TauRes + lam))
  h2 <- exp(-MaxAge * (1 / TauRes + lam))
  out <- Cin[j] * mult * (h1 - h2)
  MinAge <- MaxAge
  for (i in 1:2000000) {
    h1 <- h2
    MaxAge <- MinAge + i * TI
    while (End - MaxAge < dates[j] && j != 1) j <- j - 1
    prev <- out
    h2 <- exp(-MaxAge * (1 / TauRes + lam))
    out <- out + Cin[j] * mult * (h1 - h2)
    if (MaxAge >= Tau && out > 0 && (out - prev) / out < 1e-6) break
  }
  out
}
vba_PFM <- function(dates, Cin, Tau, S, lam) {
  End <- S - Tau
  j <- length(dates)
  while (End < dates[j] && j != 1) j <- j - 1
  Cin[j] * exp(-lam * Tau)
}
dates <- 1900:2020
Cin <- 10 + 400 * exp(-0.5 * ((dates - 1963) / 3)^2)
lam <- log(2) / 12.32
for (o in c(1, 0.5)) {
  e_emm <- e_pfm <- numeric(0)
  for (k in seq(61, 121, by = 10)) {
    S <- dates[k] + o
    xk <- Cin[1:k]
    w <- tlpm_kernel_w("EMM", 25, NA, lam, 1, o, 0, k)
    y <- apply_kernel(xk, w$w, Cin[1], w$total)[k]
    e_emm <- c(e_emm, abs(y / vba_EMM(dates, Cin, 25, S, lam) - 1))
    w <- tlpm_kernel_w("PFM", 25.3, NA, lam, 1, o, 0, k)
    y <- apply_kernel(xk, w$w, Cin[1], w$total)[k]
    e_pfm <- c(e_pfm, abs(y / vba_PFM(dates, Cin, 25.3, S, lam) - 1))
  }
  results <- c(results, ok(max(e_emm) < 5e-4, sprintf("TracerLPM vs EMM_VB (o = %.1f) : écart relatif max %.1e", o, max(e_emm))))
  results <- c(results, ok(max(e_pfm) < 1e-10, sprintf("TracerLPM vs PFM_VB (o = %.1f) : écart relatif max %.1e", o, max(e_pfm))))
}

for (o in c(1, 0.5)) {
  pa <- make_prep(data.frame(obs = Cin, Cin = Cin, date = NA_real_), "3h", 1, "none", 0L,
                  as.character(o), 1000, 1e-7, 1e6, 0)
  ks <- seq(61, 121, by = 10)
  b_emm <- tlpm_back_component("EMM", 25, NA, pa, at = ks)[ks]
  v_emm <- vapply(ks, function(k) vba_EMM(dates, Cin, 25, dates[k] + o, lam), numeric(1))
  b_pfm <- tlpm_back_component("PFM", 25.3, NA, pa, at = ks)[ks]
  v_pfm <- vapply(ks, function(k) vba_PFM(dates, Cin, 25.3, dates[k] + o, lam), numeric(1))
  results <- c(results, ok(max(abs(b_emm / v_emm - 1)) < 1e-10,
                           sprintf("TracerLPM-backward == EMM_VB (o = %.1f) : écart relatif max %.1e", o, max(abs(b_emm / v_emm - 1)))))
  results <- c(results, ok(max(abs(b_pfm / v_pfm - 1)) < 1e-10,
                           sprintf("TracerLPM-backward == PFM_VB (o = %.1f) : écart relatif max %.1e", o, max(abs(b_pfm / v_pfm - 1)))))
}

d3 <- read_dataset("data/exemple_3H_mensuel.csv")
p3 <- make_prep(d3, "3h", infer_dt(d3$date), "none", 0L, "1", 1000, 1e-7, 1e6, 0)
s3 <- tlpm_component("DM", 144, 0.3, p3)
rel3 <- stats::sd((d3$obs - s3) / s3, na.rm = TRUE)
results <- c(results, ok(abs(rel3 - 0.02) < 0.006,
                         sprintf("Exemple 3H : TracerLPM(DM, 12 a, DP 0,3) reproduit obs (bruit relatif %.3f, attendu 0,02)", rel3)))
d18 <- read_dataset("data/exemple_18O_mensuel.csv")
p18 <- make_prep(d18, "18o", infer_dt(d18$date), "u_Cin", 1000L, "1", 1000, 1e-7, 1e6, 0)
s18 <- tlpm_component("EPM", 8, 1.3, p18)
sd18 <- stats::sd(d18$obs - s18)
results <- c(results, ok(abs(sd18 - 0.05) < 0.015,
                         sprintf("Exemple 18O : TracerLPM(EPM, 8 mois, n 1,3) reproduit obs (résidu %.3f, attendu 0,05)", sd18)))

p3b <- make_prep(d3, "3h", 1 / 12, "1_Cin", 1000L, "1", 1000, 1e-7, 1e6, 0)
sb <- mitigre_component("EMM", 120, NA, p3b)
sc <- claude_component("EMM", 120, NA, p3b)
st <- tlpm_component("EMM", 120, NA, p3b)
cat(sprintf("       info : 3H mensuel EMM tau = 10 a, moyenne 2000-2020  MITIGRE / TracerLPM = %.3f ; Forward / TracerLPM = %.5f\n",
            mean(tail(sb / st, 240)), mean(tail(sc / st, 240))))
results <- c(results, ok(abs(mean(tail(sc / st, 240)) - 1) < 1e-3, "Forward ≈ TracerLPM (EMM, 3H mensuel)"))

at3 <- which(!is.na(d3$obs))
for (cs in list(list("EMM", 120, NA), list("EPM", 180, 1.3), list("PEM", 96, 2), list("DM", 144, 0.3))) {
  fw <- tlpm_component(cs[[1]], cs[[2]], cs[[3]], p3)[at3]
  bw <- tlpm_back_component(cs[[1]], cs[[2]], cs[[3]], p3, at = at3)[at3]
  results <- c(results, ok(max(abs(bw / fw - 1)) < 1e-3,
                           sprintf("Forward vs backward (%s, 3H) : écart relatif max %.1e", cs[[1]], max(abs(bw / fw - 1)))))
}

comps <- list(list(type = "EPM", tau = c(2, 20), p2 = c(1, 2), ratio = c(1, 1)))
ev <- which(!is.na(d18$obs))
r <- run_comparison(ENGINE_IDS, comps, p18, d18$obs, ev, "NSE", TRUE, 200, 1)
for (e in ENGINE_IDS) {
  b <- r$best[[e]]$idx
  cat(sprintf("       %s : NSE = %.4f, tau = %.2f pas, n = %.3f (%.2f s)\n", ENGINE_LABELS[[e]],
              r$obj[b, e], r$par$tau[b, 1], r$par$p2[b, 1], r$time[[e]]))
}
d_fb <- max(abs(r$obj[, "tlpm_back"] - r$obj[, "tracerlpm"]), na.rm = TRUE)
results <- c(results, ok(d_fb < 5e-4, sprintf("MC : |NSE backward - NSE forward| max = %.1e sur 200 runs", d_fb)))
results <- c(results, ok(all(r$obj[r$best$tracerlpm$idx, "tracerlpm"] > 0.95), "Calage MC TracerLPM sur l'exemple 18O (NSE > 0,95)"))

cat(sprintf("\n%d / %d tests réussis\n", sum(results), length(results)))
if (!all(results)) quit(status = 1)
