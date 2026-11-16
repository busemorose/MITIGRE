sa_params <- function(r) {
  nc <- length(r$types)
  out <- list()
  for (c in seq_len(nc)) {
    lab <- if (nc > 1) paste0("C", c, " ", r$types[c], " · ") else paste0(r$types[c], " · ")
    out[[length(out) + 1]] <- list(name = paste0(lab, "τ [", r$unit, "]"), value = r$par$tau[, c],
                                   prior = r$comps[[c]]$tau)
    if (uses_p2(r$types[c])) {
      out[[length(out) + 1]] <- list(name = paste0(lab, p2_name(r$types[c])), value = r$par$p2[, c],
                                     prior = r$comps[[c]]$p2)
    }
    if (nc > 1) {
      out[[length(out) + 1]] <- list(name = paste0(lab, "ratio"), value = r$par$ratio[, c],
                                     prior = c(NA_real_, NA_real_))
    }
  }
  out
}

sa_long <- function(r, pct, e) {
  ps <- sa_params(r)
  pnames <- vapply(ps, function(p) p$name, character(1))
  obj <- r$obj[, e]
  behav <- seq_len(r$n_run) %in% top_runs(obj, r$crit, pct)
  do.call(rbind, lapply(ps, function(p) {
    data.frame(run = seq_len(r$n_run), param = factor(p$name, levels = pnames), value = p$value,
               obj = obj, behav = behav)
  }))
}

crit_decile <- function(obj, crit) {
  rk <- rank(if (crit == "RMSE") obj else -obj, na.last = "keep", ties.method = "first")
  n <- sum(!is.na(obj))
  factor(ceiling(10 * rk / n), levels = 1:10)
}

ks_distance <- function(x, behav) {
  xb <- x[behav & !is.na(x)]
  xn <- x[!behav & !is.na(x)]
  if (length(xb) < 2 || length(xn) < 2) return(NA_real_)
  v <- sort(unique(c(xb, xn)))
  max(abs(stats::ecdf(xb)(v) - stats::ecdf(xn)(v)))
}

sa_table <- function(r, pct, e) {
  obj <- r$obj[, e]
  top <- top_runs(obj, r$crit, pct)
  behav <- seq_len(r$n_run) %in% top
  do.call(rbind, lapply(sa_params(r), function(p) {
    q <- stats::quantile(p$value[top], c(0.05, 0.5, 0.95), na.rm = TRUE)
    data.frame(`Paramètre` = p$name,
               `A priori` = if (anyNA(p$prior)) "—" else sprintf("%.4g – %.4g", p$prior[1], p$prior[2]),
               `Médiane` = sprintf("%.4g", q[2]),
               `5–95 %` = sprintf("%.4g – %.4g", q[1], q[3]),
               KS = sprintf("%.2f", ks_distance(p$value, behav)),
               check.names = FALSE)
  }))
}

bound_warnings <- function(r, e) {
  i <- r$best[[e]]$idx
  if (is.na(i)) return(character(0))
  nc <- length(r$types)
  chk <- function(v, rng, name, fac = 1, unit = "") {
    if (anyNA(rng) || diff(rng) <= 0) return(NULL)
    f <- (v - rng[1]) / diff(rng)
    if (f > 0.97) return(sprintf("%s en butée haute (%.4g%s)", name, rng[2] * fac, unit))
    if (f < 0.03) return(sprintf("%s en butée basse (%.4g%s)", name, rng[1] * fac, unit))
    NULL
  }
  out <- character(0)
  for (c in seq_len(nc)) {
    cp <- r$comps[[c]]
    pre <- if (nc > 1) paste0("C", c, " ") else ""
    out <- c(out, chk(r$par$tau[i, c], cp$tau, paste0(pre, "τ"), 1, paste0(" ", r$unit)))
    if (uses_p2(r$types[c])) out <- c(out, chk(r$par$p2[i, c], cp$p2, paste0(pre, p2_name(r$types[c]))))
  }
  out
}
