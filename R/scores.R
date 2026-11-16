KGE <- function(sim, obs, na.rm = TRUE, allow_NA = FALSE) {
  if (allow_NA) {
    na <- is.na(obs)
    obs <- obs[!na]
    sim <- sim[!na]
  }
  mean_s <- mean(sim, na.rm = na.rm)
  mean_o <- mean(obs, na.rm = na.rm)
  sd_s <- stats::sd(sim, na.rm = na.rm)
  sd_o <- stats::sd(obs, na.rm = na.rm)
  alpha <- sd_s / sd_o
  beta <- mean_s / mean_o
  rpearson <- stats::cor(sim, obs, method = "pearson")
  1 - sqrt((rpearson - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}

NSE <- function(sim, obs, na.rm = TRUE, allow_NA = FALSE) {
  if (allow_NA) {
    na <- is.na(obs)
    obs <- obs[!na]
    sim <- sim[!na]
  }
  mean_o <- mean(obs, na.rm = na.rm)
  1 - sum((sim - obs)^2) / sum((obs - mean_o)^2)
}

KGENP <- function(sim, obs, na.rm = TRUE, allow_NA = FALSE) {
  if (allow_NA) {
    na <- is.na(obs)
    obs <- obs[!na]
    sim <- sim[!na]
  }
  mean_s <- mean(sim, na.rm = na.rm)
  mean_o <- mean(obs, na.rm = na.rm)
  fdc_s <- sort(sim / (mean_s * length(sim)))
  fdc_o <- sort(obs / (mean_o * length(obs)))
  alphaNP <- 1 - 0.5 * sum(abs(fdc_s - fdc_o))
  beta <- mean_s / mean_o
  rspearman <- stats::cor(sim, obs, method = "spearman")
  1 - sqrt((rspearman - 1)^2 + (alphaNP - 1)^2 + (beta - 1)^2)
}

KGEABS <- function(sim, obs, na.rm = TRUE, allow_NA = FALSE) {
  if (allow_NA) {
    na <- is.na(obs)
    obs <- obs[!na]
    sim <- sim[!na]
  }
  sd_s <- stats::sd(sim, na.rm = na.rm)
  sd_o <- stats::sd(obs, na.rm = na.rm)
  dif <- sim - obs
  beta_abs <- sum(abs(dif)) / sum(obs)
  alpha <- sd_s / sd_o
  rpearson <- stats::cor(sim, obs, method = "pearson")
  1 - sqrt((rpearson - 1)^2 + (alpha - 1)^2 + beta_abs^2)
}

RMSE <- function(sim, obs, na.rm = TRUE, allow_NA = FALSE) {
  if (allow_NA) {
    na <- is.na(obs)
    obs <- obs[!na]
    sim <- sim[!na]
  }
  sqrt(sum((sim - obs)^2) / length(obs))
}

CRITERIA <- c("KGE", "NSE", "KGENP", "KGEABS", "RMSE")

crit_value <- function(crit, sim, obs, allow_NA) {
  f <- switch(crit, KGE = KGE, NSE = NSE, KGENP = KGENP, KGEABS = KGEABS, RMSE = RMSE)
  suppressWarnings(f(sim, obs, allow_NA = allow_NA))
}

best_index <- function(obj, crit) {
  if (all(is.na(obj))) return(NA_integer_)
  if (crit == "RMSE") which.min(obj) else which.max(obj)
}

top_runs <- function(obj, crit, pct) {
  if (is.null(pct) || length(pct) != 1 || is.na(pct) || pct <= 0) pct <- 100
  pct <- min(pct, 100)
  ok <- which(!is.na(obj))
  if (!length(ok)) return(integer(0))
  ord <- ok[order(obj[ok], decreasing = crit != "RMSE")]
  ord[seq_len(max(1L, ceiling(length(ord) * pct / 100)))]
}

score_table <- function(sim, obs, allow_NA = TRUE) {
  if (allow_NA) {
    na <- is.na(obs)
    obs <- obs[!na]
    sim <- sim[!na]
  }
  mean_s <- mean(sim)
  mean_o <- mean(obs)
  sd_s <- stats::sd(sim)
  sd_o <- stats::sd(obs)
  r <- suppressWarnings(stats::cor(sim, obs))
  alpha <- sd_s / sd_o
  beta <- mean_s / mean_o
  beta_abs <- sum(abs(sim - obs)) / sum(obs)
  rs <- suppressWarnings(stats::cor(sim, obs, method = "spearman"))
  fdc_s <- sort(sim / (mean_s * length(sim)))
  fdc_o <- sort(obs / (mean_o * length(obs)))
  alphaNP <- 1 - 0.5 * sum(abs(fdc_s - fdc_o))
  c(NSE = 1 - sum((sim - obs)^2) / sum((obs - mean_o)^2),
    KGE = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2),
    KGENP = 1 - sqrt((rs - 1)^2 + (alphaNP - 1)^2 + (beta - 1)^2),
    KGE_abs = 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + beta_abs^2),
    RMSE = sqrt(mean((sim - obs)^2)),
    alpha = alpha, beta = beta, r = r)
}
