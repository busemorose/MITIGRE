MITIGRE <- function(...) {
  app_dir <- system.file("app", package = "MITIGRE")
  if (!nzchar(app_dir)) stop("Application not found: try re-installing MITIGRE.", call. = FALSE)
  shiny::runApp(app_dir, display.mode = "normal", ...)
}
