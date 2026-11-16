APP_ENGINES <- c("claude", "tlpm_back")
ENGINE_IDS <- c("mitigre", "claude", "tracerlpm", "tlpm_back")
ENGINE_LABELS <- c(mitigre = "MITIGRE (original)", claude = "Forward",
                   tracerlpm = "TracerLPM forward (exact)", tlpm_back = "Backward")
ENGINE_COLORS <- c("MITIGRE (original)" = "#009E73", "Forward" = "#0072B2",
                   "TracerLPM forward (exact)" = "#E69F00", "Backward" = "#D55E00")
POINTWISE_ENGINES <- "tlpm_back"

MODEL_TYPES <- c("EMM", "EPM", "PEM", "DM", "PFM")
ISOTOPES <- c("18O" = "18o", "2H" = "2h", "3H" = "3h", "14C" = "14c")
HALF_LIFE <- c("18o" = Inf, "2h" = Inf, "3h" = 12.32, "14c" = 5730)

lambda_year <- function(isotope) {
  hl <- HALF_LIFE[[isotope]]
  if (is.infinite(hl)) 0 else log(2) / hl
}

uses_p2 <- function(type) type %in% c("EPM", "PEM", "DM")

p2_name <- function(type) {
  switch(type, EPM = "n", PEM = "n", DM = "DP", "-")
}

decimal_year <- function(d) {
  y <- as.integer(format(d, "%Y"))
  start <- as.Date(paste0(y, "-01-01"))
  end <- as.Date(paste0(y + 1L, "-01-01"))
  y + as.numeric(d - start) / as.numeric(end - start)
}

parse_dates <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  filled <- !is.na(x) & x != ""
  num <- suppressWarnings(as.numeric(x))
  if (any(filled) && all(!is.na(num[filled]))) return(num)
  for (fmt in c("%Y-%m-%d", "%d/%m/%Y", "%Y/%m/%d", "%d-%m-%Y", "%d.%m.%Y")) {
    d <- as.Date(x, format = fmt)
    if (sum(!is.na(d)) >= 0.9 * sum(filled)) return(decimal_year(d))
  }
  rep(NA_real_, length(x))
}

read_dataset <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (grepl(";", first, fixed = TRUE)) ";" else if (grepl("\t", first, fixed = TRUE)) "\t" else ","
  df <- utils::read.table(path, header = TRUE, sep = sep, dec = ".", stringsAsFactors = FALSE,
                          na.strings = c("NA", "", "NaN", "na"), check.names = FALSE,
                          quote = "\"", comment.char = "")
  low <- tolower(trimws(names(df)))
  if (!all(c("obs", "cin") %in% low)) {
    stop("Le fichier doit contenir au moins les colonnes 'obs' et 'Cin'.")
  }
  out <- data.frame(obs = suppressWarnings(as.numeric(df[[which(low == "obs")[1]]])),
                    Cin = suppressWarnings(as.numeric(df[[which(low == "cin")[1]]])))
  out$date <- if ("date" %in% low) parse_dates(df[[which(low == "date")[1]]]) else NA_real_
  out
}

DT_CHOICES <- c("Jour" = "day", "Semaine" = "week", "Quinzaine" = "half_month", "Mois" = "month",
                "Trimestre" = "quarter", "Semestre" = "half_year", "Année" = "year")
DT_YEARS <- c(day = 1 / 365.25, week = 7 / 365.25, half_month = 1 / 24, month = 1 / 12,
              quarter = 0.25, half_year = 0.5, year = 1)

DT_UNITS <- c(day = "jours", week = "semaines", half_month = "quinzaines", month = "mois",
              quarter = "trimestres", half_year = "semestres", year = "années")

nearest_dt_code <- function(dt) names(DT_YEARS)[which.min(abs(log(DT_YEARS / dt)))]

infer_dt <- function(date) {
  if (all(is.na(date))) return(NA_real_)
  d <- stats::median(diff(date[!is.na(date)]))
  if (!is.finite(d) || d <= 0) return(NA_real_)
  usual <- c(1 / 365.25, 7 / 365.25, 1 / 24, 1 / 12, 0.25, 0.5, 1)
  k <- which.min(abs(usual - d) / usual)
  if (abs(usual[k] - d) / usual[k] < 0.1) usual[k] else d
}

build_input <- function(obs, Cin, method, n_warm) {
  if (method == "none") n_warm <- 0L
  val <- switch(method,
                u_Cin = mean(Cin, na.rm = TRUE),
                u_Cobs = mean(obs, na.rm = TRUE),
                "1_Cin" = Cin[1],
                "1_Cobs" = obs[1],
                "0" = 0,
                none = Cin[1])
  if (is.na(val)) val <- Cin[1]
  list(x = c(rep(val, n_warm), Cin), n_warm = as.integer(n_warm), c0 = val)
}

conv_causal <- function(x, w) {
  n <- length(x)
  if (length(w) > n) w <- w[seq_len(n)]
  m <- length(w)
  if (m == 0L) return(numeric(n))
  L <- stats::nextn(n + m - 1L)
  X <- stats::fft(c(x, numeric(L - n)))
  W <- stats::fft(c(w, numeric(L - m)))
  Re(stats::fft(X * W, inverse = TRUE))[seq_len(n)] / L
}

apply_kernel <- function(x, w, c0, total) {
  n <- length(x)
  y <- conv_causal(x, w)
  if (c0 != 0) {
    cw <- if (length(w)) cumsum(w) else 0
    cw <- if (length(cw) >= n) cw[seq_len(n)] else c(cw, rep(cw[length(cw)], n - length(cw)))
    y <- y + c0 * (total - cw)
  }
  y
}

pinvgauss <- function(x, mu, lambda) {
  out <- numeric(length(x))
  pos <- is.finite(x) & x > 0
  if (any(pos)) {
    xs <- x[pos]
    r <- sqrt(lambda / xs)
    a <- stats::pnorm(r * (xs / mu - 1))
    b <- exp(2 * lambda / mu + stats::pnorm(-r * (xs / mu + 1), log.p = TRUE))
    out[pos] <- pmin(a + b, 1)
  }
  out[is.infinite(x) & x > 0] <- 1
  out
}

support_length <- function(type, tau, p2, tol) {
  switch(type,
         EMM = tau * log(1 / tol),
         EPM = tau * (1 - 1 / p2) + tau / p2 * log(1 / tol),
         PEM = {
           taq <- tau / (1 + log(p2))
           taq * log(p2) + taq * log(1 / tol)
         },
         DM = {
           a <- tau
           while (1 - pinvgauss(a, tau, tau / (2 * p2)) > tol && a < 1e9 * tau) a <- a * 2
           a
         },
         PFM = tau)
}

sample_params <- function(comps, n_run, seed) {
  set.seed(seed)
  nc <- length(comps)
  tau <- p2 <- ratio <- matrix(NA_real_, n_run, nc)
  for (c in seq_len(nc)) {
    cp <- comps[[c]]
    tau[, c] <- stats::runif(n_run, cp$tau[1], cp$tau[2])
    if (uses_p2(cp$type)) p2[, c] <- stats::runif(n_run, cp$p2[1], cp$p2[2])
    ratio[, c] <- if (nc > 1) stats::runif(n_run, cp$ratio[1], cp$ratio[2]) else 1
  }
  rs <- rowSums(ratio)
  rs[rs <= 0] <- 1
  list(tau = tau, p2 = p2, ratio = ratio / rs)
}

component_fun <- function(engine) {
  switch(engine, mitigre = mitigre_component, claude = claude_component, tracerlpm = tlpm_component,
         tlpm_back = tlpm_back_component)
}

kernel_fun <- function(engine) {
  switch(engine, mitigre = mitigre_kernel, claude = claude_kernel, tracerlpm = tlpm_kernel,
         tlpm_back = tlpm_back_kernel)
}

simulate_run <- function(engine, types, tau, p2, ratio, prep, at = NULL) {
  f <- component_fun(engine)
  pointwise <- engine %in% POINTWISE_ENGINES
  out <- numeric(prep$n_obs)
  for (c in seq_along(types)) {
    if (ratio[c] > 0) {
      sim_c <- if (pointwise) f(types[c], tau[c], p2[c], prep, at) else f(types[c], tau[c], p2[c], prep)
      out <- out + ratio[c] * sim_c
    }
  }
  if (pointwise && !is.null(at)) out[setdiff(seq_len(prep$n_obs), at)] <- NA_real_
  out
}

mixture_kernel <- function(engine, types, tau, p2, ratio, prep, decay = FALSE) {
  f <- kernel_fun(engine)
  ks <- lapply(seq_along(types), function(c) f(types[c], tau[c], p2[c], prep, decay))
  nmax <- max(vapply(ks, nrow, integer(1)))
  ref <- ks[[which.max(vapply(ks, nrow, integer(1)))]]
  w <- numeric(nmax)
  for (c in seq_along(ks)) {
    k <- ks[[c]]
    w[seq_len(nrow(k))] <- w[seq_len(nrow(k))] + ratio[c] * k$w
  }
  data.frame(lag = ref$lag, lo = ref$lo, hi = ref$hi, w = w)
}

kernel_mean_age <- function(k) {
  if (sum(k$w) <= 0) return(NA_real_)
  sum(k$w * (k$lo + k$hi) / 2) / sum(k$w)
}

make_prep <- function(data, isotope, dt, warm_method, n_warm, obs_pos, max_t, tol, kmax, uz,
                      quad = "mid", tlpm_compat = FALSE) {
  inp <- build_input(data$obs, data$Cin, warm_method, n_warm)
  d0 <- if (!is.null(data$date) && !is.na(data$date[1])) data$date[1] else 2000
  dates <- d0 + (seq_along(inp$x) - 1 - inp$n_warm) * dt
  list(x = inp$x, n_warm = inp$n_warm, c0 = inp$c0, n_obs = nrow(data), dates = dates,
       isotope = isotope, lam_year = lambda_year(isotope), dt = dt,
       o = as.numeric(obs_pos), max_t = max_t, tol = tol, kmax = kmax, uz = uz, quad = quad,
       tlpm_compat = isTRUE(tlpm_compat))
}
