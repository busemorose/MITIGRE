library(shiny)
library(bslib)
library(ggplot2)

app_dir <- getwd()
for (f in sort(list.files(file.path(app_dir, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f, local = FALSE, encoding = "UTF-8")
}

WARMUP_METHODS <- c("Moyenne de Cin" = "u_Cin", "Moyenne de obs" = "u_Cobs",
                    "1re valeur de Cin" = "1_Cin", "1re valeur de obs" = "1_Cobs",
                    "Zéro" = "0")

CRIT_LABELS <- c(NSE = "NSE", KGE = "KGE", KGENP = "KGENP", KGE_abs = "KGE abs", RMSE = "RMSE",
                 alpha = "α", beta = "β", r = "r")

EVAL_COLOR <- "#E69F00"

IS_WINDOWS <- .Platform$OS.type == "windows"
DEFAULT_CORES <- if (IS_WINDOWS) 1L else max(1L, min(8L, as.integer(parallel::detectCores()) - 1L), na.rm = TRUE)

theme_app <- function(legend = "right") {
  theme_bw(base_size = 13) +
    theme(legend.position = legend, legend.title = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#F2F6F9", colour = NA))
}

no_results <- function(r) {
  validate(need(!is.null(r), "Lancez un calage depuis l'onglet Configuration."))
}

gear <- function(..., title = "Affichage") {
  bslib::popover(tags$span(icon("gear"), style = "cursor:pointer; color:#6C7A86"), ...,
                 title = title, placement = "left")
}

or_default <- function(x, default) if (isTruthy(x)) x else default

eval_layers <- function(x, obs, ev, in_legend = FALSE) {
  if (length(ev) < 2) return(list())
  d <- data.frame(x = x[ev], value = obs[ev])
  d <- d[!is.na(d$value), ]
  line <- if (in_legend) {
    geom_line(data = d, aes(x, value, colour = "eval"), linewidth = 2.5, alpha = 0.6)
  } else {
    geom_line(data = d, aes(x, value), colour = EVAL_COLOR, linewidth = 2.5, alpha = 0.6)
  }
  list(line, geom_vline(xintercept = range(x[ev]), colour = EVAL_COLOR, linetype = "dashed", linewidth = 0.6))
}

comp_ui <- function(i, input, n, unit = "mois", dt = 1 / 12, f = 1) {
  id <- function(x) paste0(x, "_", i)
  val <- function(x, default) {
    v <- input[[id(x)]]
    if (is.null(v)) default else v
  }
  tau_val <- function(x, default) {
    v <- input[[id(x)]]
    if (is.null(v)) default else signif(v * f, 4)
  }
  type_id <- id("type")
  bslib::card(
    class = "comp-card mb-0",
    bslib::card_header(paste("Composante", i)),
    bslib::card_body(
      selectInput(type_id, label_info("Modèle", TIP_MODELS), MODEL_TYPES, selected = val("type", "EMM")),
      fluidRow(
        column(6, numericInput(id("tau_min"), label_info(sprintf("τ min [%s]", unit), TIP_TAU),
                               tau_val("tau_min", 1), min = 0, step = 1)),
        column(6, numericInput(id("tau_max"), sprintf("τ max [%s]", unit),
                               tau_val("tau_max", round(20 / dt)), min = 0, step = 1))
      ),
      conditionalPanel(
        sprintf("input['%s'] == 'EPM' || input['%s'] == 'PEM'", type_id, type_id),
        fluidRow(
          column(6, numericInput(id("n_min"), "n min", val("n_min", 1), min = 1, step = 0.05)),
          column(6, numericInput(id("n_max"), "n max", val("n_max", 3), min = 1, step = 0.05))
        )
      ),
      conditionalPanel(
        sprintf("input['%s'] == 'DM'", type_id),
        fluidRow(
          column(6, numericInput(id("dp_min"), "DP min", val("dp_min", 0.01), min = 0, step = 0.01)),
          column(6, numericInput(id("dp_max"), "DP max", val("dp_max", 1), min = 0, step = 0.01))
        )
      ),
      if (n > 1) {
        fluidRow(
          column(6, numericInput(id("ratio_min"), "ratio min", val("ratio_min", 0), min = 0, max = 1, step = 0.05)),
          column(6, numericInput(id("ratio_max"), "ratio max", val("ratio_max", 1), min = 0, max = 1, step = 0.05))
        )
      }
    )
  )
}

config_page <- bslib::layout_columns(
  col_widths = c(3, 9, 4, 4, 4),

  bslib::card(
    bslib::card_header(label_info("Données", TIP_DATA)),
    fileInput("file", NULL, accept = c(".csv", ".txt", "text/plain"),
              buttonLabel = "Parcourir", placeholder = "Importer un fichier"),
    selectInput("isotope", "Traceur", ISOTOPES, selected = "18o"),
    selectInput("dt", "Pas de temps", DT_CHOICES, selected = "month")
  ),

  bslib::card(
    full_screen = TRUE,
    bslib::card_header(
      class = "d-flex justify-content-between align-items-center",
      label_info("Aperçu des données", "Sélectionner une zone pour zoomer (x et y), double-clic pour revenir à la vue entière."),
      tags$span(class = "small text-muted fw-normal", textOutput("data_summary", inline = TRUE))
    ),
    plotOutput("plot_data", height = "300px",
               brush = zoom_brush("data_brush"),
               dblclick = "data_dblclick")
  ),

  bslib::card(
    bslib::card_header("Modèle"),
    radioButtons("n_comp", "Composantes", choices = 1:4, selected = 1, inline = TRUE),
    uiOutput("ui_comps")
  ),

  bslib::card(
    bslib::card_header("Moteurs"),
    checkboxGroupInput("engines", NULL,
                       choiceNames = list(tagList("Forward", info(TIP_FORWARD)),
                                          tagList("Backward", info(TIP_BACKWARD))),
                       choiceValues = APP_ENGINES, selected = APP_ENGINES, inline = TRUE),
    tags$hr(class = "my-2"),
    selectInput("warm_method", label_info("Concentration avant la série", TIP_WARMUP), WARMUP_METHODS, selected = "u_Cin"),
    radioButtons("quad", label_info("Intégration Forward", TIP_QUAD),
                 choices = c("Point milieu" = "mid", "Exacte" = "exact"), selected = "mid", inline = TRUE),
    radioButtons("obs_pos", label_info("Observation datée", TIP_OBSPOS),
                 choices = c("Début" = "0", "Milieu" = "0.5", "Fin du pas" = "1"), selected = "1", inline = TRUE),
    checkboxInput("tlpm_compat", label_info("Backward : compatibilité TracerLPM (DM)", TIP_COMPAT), value = FALSE),
    fluidRow(
      column(6, numericInput("tol", label_info("Tolérance", TIP_TOL), value = 1e-7, min = 1e-15, step = 1e-7)),
      column(6, numericInput("kmax", "K max [pas]", value = 1e6, min = 10, step = 1e5))
    ),
    numericInput("cores", label_info("Cœurs de calcul", TIP_CORES), value = DEFAULT_CORES, min = 1, step = 1)
  ),

  bslib::card(
    bslib::card_header("Calage"),
    fluidRow(
      column(6, numericInput("n_run", "Itérations", value = 1000, min = 1, max = 1e5, step = 100)),
      column(6, numericInput("seed", label_info("Graine aléatoire", TIP_SEED), value = 1, min = 0, step = 1))
    ),
    selectInput("crit", "Critère", CRITERIA, selected = "KGE"),
    checkboxInput("allow_NA", label_info("Ignorer les NA", "Les dates sans observation sont exclues du calcul du critère."), value = TRUE),
    uiOutput("ui_eval_range"),
    bslib::card_footer(
      actionButton("run", "Lancer le calage", icon = icon("play"), class = "btn-primary btn-lg w-100")
    )
  )
)

results_page <- tagList(
  uiOutput("vboxes"),
  bslib::card(
    full_screen = TRUE,
    bslib::card_header(
      class = "d-flex justify-content-between align-items-center",
      tags$span("Simulations", info("Sélectionner une zone pour zoomer (x et y), double-clic pour revenir à la vue entière.")),
      gear(checkboxInput("band", "Enveloppe des meilleurs runs", value = FALSE),
           numericInput("band_pct", "Meilleurs runs [%]", value = 10, min = 0.1, max = 100, step = 1),
           sliderInput("band_q", "Quantiles", min = 0, max = 1, value = c(0.05, 0.95), step = 0.01))
    ),
    plotOutput("plot_sim", height = "420px",
               brush = zoom_brush("sim_brush"),
               dblclick = "sim_dblclick")
  ),
  bslib::layout_columns(
    col_widths = c(5, 7),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(label_info("Performances", "Critères sur la période d'évaluation, au meilleur run de chaque moteur.")),
      bslib::card_body(
        class = "results-table",
        tableOutput("tab_perf"),
        tableOutput("tab_params")
      )
    ),
    bslib::navset_card_underline(
      title = "Fonction de transfert",
      full_screen = TRUE,
      bslib::nav_panel("Densité", plotOutput("plot_kern", height = "340px")),
      bslib::nav_panel("Cumulée", plotOutput("plot_kern_cum", height = "340px")),
      bslib::nav_spacer(),
      bslib::nav_item(gear(checkboxInput("kern_decay", "Inclure la décroissance", value = FALSE),
                           checkboxInput("kern_log", "Âge en échelle log", value = FALSE),
                           numericInput("kern_max", "Âge max (unité du pas)", value = NA, min = 0)))
    )
  )
)

sens_page <- tagList(
  tags$div(
    class = "sa-bar d-flex align-items-center mb-3 px-1",
    radioButtons("sa_engine", NULL, choiceNames = unname(ENGINE_LABELS[APP_ENGINES]),
                 choiceValues = APP_ENGINES, inline = TRUE),
    tags$div(class = "ms-auto",
             gear(sliderInput("sa_pct", label_info("Runs retenus [% meilleurs]",
                                                   "Les X % meilleurs runs selon le critère de calage : utilisés pour les distributions, les interactions et la synthèse."),
                              min = 1, max = 100, value = 10, step = 1),
                  title = "Seuil de comportement"))
  ),
  bslib::card(
    full_screen = TRUE,
    bslib::card_header(label_info("Critère en fonction de chaque paramètre",
                                  "Chaque point est un run Monte Carlo. Une enveloppe supérieure marquée signale un paramètre influent ; un nuage plat, un paramètre peu influent.")),
    plotOutput("sa_dotty", height = "300px")
  ),
  bslib::layout_columns(
    col_widths = c(6, 6),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(label_info("Distributions des runs retenus",
                                    "Histogramme des runs retenus ; tirets : distribution a priori (uniforme). Une distribution resserrée indique un paramètre bien contraint.")),
      plotOutput("sa_hist", height = "300px")
    ),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(label_info("Sensibilité régionale",
                                    "Fonctions de répartition du paramètre par décile du critère (1 = 10 % meilleurs runs). Des courbes écartées indiquent un paramètre influent (Hornberger & Spear, 1981 ; Freer et al., 1996).")),
      plotOutput("sa_rsa", height = "300px")
    )
  ),
  bslib::layout_columns(
    col_widths = c(6, 6),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(
        class = "d-flex justify-content-between align-items-center",
        label_info("Interactions entre paramètres",
                   "Runs colorés par le critère ; cercles noirs : runs retenus. Une structure en bande révèle une compensation entre paramètres (équifinalité)."),
        gear(uiOutput("ui_sa_xy"), title = "Paramètres affichés")
      ),
      plotOutput("sa_xy", height = "320px")
    ),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(label_info("Synthèse",
                                    "Médiane et intervalle 5–95 % des runs retenus. KS : distance de Kolmogorov-Smirnov entre runs retenus et non retenus (0 : sans influence, 1 : très influent).")),
      bslib::card_body(class = "results-table", tableOutput("sa_table"))
    )
  )
)

explore_page <- explore_page_ui()

export_ui <- bslib::popover(
  tags$span(id = "export_icon", class = "banner-icon export-off", title = "Exporter les résultats",
            style = "cursor:pointer", icon("download")),
  tags$div(
    class = "export-menu",
    downloadLink("dl_best", tagList(icon("file-csv"), " Meilleures simulations")),
    downloadLink("dl_runs", tagList(icon("table"), " Tous les runs Monte Carlo")),
    downloadLink("dl_kern", tagList(icon("wave-square"), " Fonctions de transfert")),
    downloadLink("dl_sens", tagList(icon("chart-area"), " Synthèse de sensibilité"))
  ),
  title = "Exporter (CSV ;)", placement = "bottom"
)

ui <- bslib::page_navbar(
  title = "MITIGRE",
  id = "main",
  window_title = "MITIGRE",
  theme = bslib::bs_theme(version = 5, primary = BANNER_COLOR),
  header = tagList(app_css(), app_js()),
  bslib::nav_panel("Configuration", icon = icon("sliders"), config_page),
  bslib::nav_panel("Résultats", icon = icon("chart-line"), results_page),
  bslib::nav_panel("Sensibilité", icon = icon("chart-area"), sens_page),
  bslib::nav_panel("Exploration", icon = icon("flask"), explore_page),
  bslib::nav_spacer(),
  bslib::nav_item(export_ui),
  bslib::nav_item(actionLink("notes_btn", icon("book-open"), class = "banner-icon notes", title = "Notes"))
)

server <- function(input, output, session) {

  res <- reactiveVal(NULL)
  set_results <- function(x) {
    res(x)
    session$sendCustomMessage("export_state", !is.null(x))
  }

  zoom_state <- function(brush_id, dbl_id) {
    z <- reactiveVal(NULL)
    observeEvent(input[[brush_id]], {
      b <- input[[brush_id]]
      z(list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax)))
    })
    observeEvent(input[[dbl_id]], z(NULL))
    z
  }
  zoom_coord <- function(z) if (is.null(z)) coord_cartesian() else coord_cartesian(xlim = z$x, ylim = z$y)
  data_zoom <- zoom_state("data_brush", "data_dblclick")
  ex_zoom <- zoom_state("ex_brush", "ex_dblclick")

  data_path <- reactiveVal(NULL)
  observeEvent(input$file, data_path(input$file$datapath))
  observeEvent(data_path(), {
    set_results(NULL)
    data_zoom(NULL)
    ex_zoom(NULL)
  })

  data_r <- reactive({
    req(data_path())
    d <- tryCatch(read_dataset(data_path()), error = function(e) {
      showNotification(conditionMessage(e), type = "error")
      NULL
    })
    req(d)
    if (anyNA(d$Cin)) {
      ok <- which(!is.na(d$Cin))
      validate(need(length(ok) >= 2, "Cin doit contenir au moins deux valeurs."))
      d$Cin <- stats::approx(ok, d$Cin[ok], seq_len(nrow(d)), rule = 2)$y
      showNotification("Cin contenait des valeurs manquantes : interpolation linéaire.", type = "warning")
    }
    d
  })

  has_dates <- reactive(!anyNA(data_r()$date))

  observeEvent(data_r(), {
    dt <- infer_dt(data_r()$date)
    if (!is.na(dt)) updateSelectInput(session, "dt", selected = nearest_dt_code(dt))
  })
  dt_r <- reactive(DT_YEARS[[input$dt]])
  unit_r <- reactive(DT_UNITS[[input$dt]])
  n_cores <- reactive(if (IS_WINDOWS) 1L else max(1L, as.integer(or_default(input$cores, 1))))

  x_axis <- reactive({
    d <- data_r()
    if (isTRUE(has_dates())) list(x = d$date, lab = "Année") else list(x = seq_len(nrow(d)), lab = "t [pas]")
  })

  n_comp <- reactive(as.integer(input$n_comp))

  last_dt_comps <- DT_YEARS[["month"]]
  output$ui_comps <- renderUI({
    n <- n_comp()
    req(n >= 1, n <= 4)
    dt <- dt_r()
    f <- last_dt_comps / dt
    last_dt_comps <<- dt
    cards <- isolate(lapply(seq_len(n), comp_ui, input = input, n = n, unit = unit_r(), dt = dt, f = f))
    do.call(bslib::layout_column_wrap, c(list(width = 1, gap = "0.5rem"), cards))
  })

  comps_r <- reactive({
    n <- n_comp()
    lapply(seq_len(n), function(i) {
      g <- function(x) input[[paste0(x, "_", i)]]
      type <- g("type")
      req(type)
      list(type = type,
           tau = c(g("tau_min"), g("tau_max")),
           p2 = switch(type, EPM = c(g("n_min"), g("n_max")), PEM = c(g("n_min"), g("n_max")),
                       DM = c(g("dp_min"), g("dp_max")), c(NA_real_, NA_real_)),
           ratio = if (n > 1) c(g("ratio_min"), g("ratio_max")) else c(1, 1))
    })
  })

  output$ui_eval_range <- renderUI({
    d <- data_r()
    ok <- which(!is.na(d$obs))
    req(length(ok) > 0)
    if (isTRUE(has_dates())) {
      lo <- min(d$date); hi <- max(d$date); v <- c(d$date[min(ok)], d$date[max(ok)])
      step <- isolate(dt_r()); lab <- "Période d'évaluation"
    } else {
      lo <- 1; hi <- nrow(d); v <- c(min(ok), max(ok)); step <- 1; lab <- "Période d'évaluation [pas]"
    }
    tagList(
      sliderInput("eval_range", lab, min = lo, max = hi, value = v, step = step, sep = ""),
      fluidRow(
        column(6, numericInput("eval_from", "Début", value = round(v[1], 4), min = lo, max = hi, step = step)),
        column(6, numericInput("eval_to", "Fin", value = round(v[2], 4), min = lo, max = hi, step = step))
      )
    )
  })

  eval_tol <- reactive(if (isTRUE(has_dates())) 0.5 * dt_r() else 0.5)
  observeEvent(input$eval_range, {
    rg <- input$eval_range
    if (!isTruthy(input$eval_from) || abs(input$eval_from - rg[1]) > eval_tol()) {
      updateNumericInput(session, "eval_from", value = round(rg[1], 4))
    }
    if (!isTruthy(input$eval_to) || abs(input$eval_to - rg[2]) > eval_tol()) {
      updateNumericInput(session, "eval_to", value = round(rg[2], 4))
    }
  })
  observeEvent(list(input$eval_from, input$eval_to), {
    req(input$eval_range, isTruthy(input$eval_from), isTruthy(input$eval_to))
    v <- sort(c(input$eval_from, input$eval_to))
    if (any(abs(v - input$eval_range) > eval_tol())) updateSliderInput(session, "eval_range", value = v)
  })

  eval_idx_r <- reactive({
    d <- data_r()
    rg <- input$eval_range
    req(rg)
    x <- x_axis()$x
    tol <- if (isTRUE(has_dates())) 0.5 * infer_dt(d$date) else 0
    which(x >= rg[1] - tol & x <= rg[2] + tol)
  })

  output$plot_data <- renderPlot({
    validate(need(!is.null(data_path()), "Importer un fichier de données."))
    d <- data_r()
    ax <- x_axis()
    ev <- if (is.null(input$eval_range)) integer(0) else eval_idx_r()
    cin <- data.frame(x = ax$x, value = d$Cin)
    obs <- data.frame(x = ax$x, value = d$obs)
    obs_ok <- obs[!is.na(obs$value), ]
    ggplot() +
      geom_line(data = cin, aes(x, value, colour = "cin"), linewidth = 0.4) +
      eval_layers(ax$x, d$obs, ev, in_legend = TRUE) +
      geom_line(data = obs_ok, aes(x, value, colour = "obs"), linewidth = 0.3) +
      geom_point(data = obs_ok, aes(x, value, colour = "obs"), size = 1) +
      scale_colour_manual(values = c(cin = "grey60", obs = "black", eval = EVAL_COLOR),
                          breaks = c("cin", "obs", "eval"),
                          labels = c("Entrée (Cin)", "Observations", "Période d'évaluation")) +
      zoom_coord(data_zoom()) +
      labs(x = ax$lab, y = names(ISOTOPES)[ISOTOPES == input$isotope]) +
      theme_app()
  })

  output$data_summary <- renderText({
    req(data_path())
    d <- data_r()
    per <- if (isTRUE(has_dates())) sprintf(" · %.2f – %.2f", min(d$date), max(d$date)) else ""
    sprintf("%d pas · %d observations%s", nrow(d), sum(!is.na(d$obs)), per)
  })

  prep_r <- reactive({
    req(isTruthy(input$tol), isTruthy(input$kmax))
    make_prep(data_r(), input$isotope, dt_r(), input$warm_method,
              0L, input$obs_pos, quad = or_default(input$quad, "mid"),
              tlpm_compat = isTRUE(input$tlpm_compat),
              max_t = 1000, tol = input$tol, kmax = as.integer(input$kmax), uz = 0)
  })

  observeEvent(input$run, {
    if (is.null(data_path())) {
      showNotification("Importer un fichier de données.", type = "error")
      return()
    }
    d <- data_r()
    engines <- intersect(APP_ENGINES, input$engines)
    comps <- comps_r()
    ev <- eval_idx_r()
    err <- c(
      if (!length(engines)) "Sélectionner au moins un moteur.",
      validate_comps(comps),
      if (length(ev) < 2) "Période d'évaluation trop courte.",
      if (!isTruthy(input$n_run) || input$n_run < 1) "Nombre d'itérations invalide.",
      if (!isTruthy(input$tol) || input$tol <= 0 || input$tol >= 1) "Tolérance invalide (0 < tol < 1).",
      if (!isTruthy(input$kmax) || input$kmax < 10) "K max invalide."
    )
    if (length(err)) {
      showNotification(paste(err, collapse = " "), type = "error")
      return()
    }
    prep <- prep_r()
    n_run <- max(1L, as.integer(input$n_run))
    total <- n_run * length(engines)
    out <- withProgress(message = "Monte Carlo", value = 0, {
      run_comparison(engines, comps, prep, d$obs, ev, input$crit, input$allow_NA, n_run,
                     or_default(input$seed, 1),
                     cores = n_cores(),
                     progress = function(e, r, n_done) {
                       incProgress(n_done / total, detail = sprintf("%s : %d / %d", ENGINE_LABELS[[e]], r, n_run))
                     })
    })
    out$data <- d
    out$x <- x_axis()
    out$isotope <- input$isotope
    out$comps <- comps
    out$unit <- unit_r()
    set_results(out)
    if (length(out$errors)) showNotification(paste(out$errors, collapse = " ; "), type = "warning")
    for (e in engines) {
      w <- bound_warnings(out, e)
      if (length(w)) {
        showNotification(sprintf("%s : %s. Élargir la plage de tirage.", ENGINE_LABELS[[e]], paste(w, collapse = ", ")),
                         type = "warning", duration = 12)
      }
    }
    updateRadioButtons(session, "sa_engine", choiceNames = unname(ENGINE_LABELS[engines]),
                       choiceValues = engines, selected = engines[1], inline = TRUE)
    sim_range(NULL)
    bslib::nav_select("main", "Résultats", session = session)
  })

  best_kernel <- function(r, e, decay = FALSE) {
    i <- r$best[[e]]$idx
    if (is.na(i)) return(NULL)
    mixture_kernel(e, r$types, r$par$tau[i, ], r$par$p2[i, ], r$par$ratio[i, ], r$prep, decay = decay)
  }

  output$vboxes <- renderUI({
    r <- res()
    if (is.null(r)) {
      return(tags$div(class = "text-muted text-center py-3",
                      icon("circle-play"), " Lancez un calage depuis l'onglet Configuration."))
    }
    boxes <- lapply(r$engines, function(e) {
      b <- r$best[[e]]
      lab <- ENGINE_LABELS[[e]]
      col <- ENGINE_COLORS[[lab]]
      if (is.na(b$idx)) {
        return(bslib::value_box(title = lab, value = "—", "Aucun run valide"))
      }
      age <- kernel_mean_age(best_kernel(r, e)) / r$prep$dt
      w <- bound_warnings(r, e)
      bslib::value_box(
        title = paste(lab, "·", r$crit),
        value = sprintf("%.3f", r$obj[b$idx, e]),
        showcase = icon(if (e == "tlpm_back") "backward" else "forward"),
        tags$p(sprintf("Âge moyen %.4g %s · run n°%d · %.1f s", age, r$unit, b$idx, r$time[[e]])),
        if (length(w)) tags$p(class = "bound-warn", icon("triangle-exclamation"), paste(w, collapse = " ; ")),
        style = sprintf("background-color:%s; color:#fff; --bslib-color-bg:%s; --bslib-color-fg:#fff;", col, col)
      )
    })
    do.call(bslib::layout_columns, c(list(col_widths = rep(12 / length(boxes), length(boxes))), boxes))
  })

  sim_range <- reactiveVal(NULL)
  observeEvent(input$sim_brush, {
    b <- input$sim_brush
    sim_range(list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax)))
  })
  observeEvent(input$sim_dblclick, sim_range(NULL))

  output$plot_sim <- renderPlot({
    r <- res()
    no_results(r)
    d <- r$data
    x <- r$x$x
    rg <- sim_range()
    idx <- if (is.null(rg)) seq_along(x) else which(x >= rg$x[1] & x <= rg$x[2])
    validate(need(length(idx) > 1, "Sélection trop courte : double-clic pour revenir à la vue entière."))
    p <- ggplot()
    if (isTRUE(input$band) && isTRUE(r$keep_sims)) {
      q_probs <- or_default(input$band_q, c(0.05, 0.95))
      bands <- do.call(rbind, lapply(r$engines, function(e) {
        top <- top_runs(r$obj[, e], r$crit, or_default(input$band_pct, 10))
        if (length(top) < 2) return(NULL)
        q <- apply(r$sims[[e]][top, idx, drop = FALSE], 2, stats::quantile, probs = q_probs, na.rm = TRUE)
        data.frame(x = x[idx], lo = q[1, ], hi = q[2, ], engine = ENGINE_LABELS[[e]])
      }))
      if (!is.null(bands)) p <- p + geom_ribbon(data = bands, aes(x, ymin = lo, ymax = hi, fill = engine), alpha = 0.2)
    }
    sims <- do.call(rbind, lapply(r$engines, function(e) {
      data.frame(x = x[idx], value = r$best[[e]]$sim[idx], engine = ENGINE_LABELS[[e]])
    }))
    obs <- data.frame(x = x[idx], value = d$obs[idx])
    p +
      eval_layers(x, d$obs, r$eval_idx) +
      geom_point(data = obs, aes(x, value), size = 1.3, colour = "grey15", na.rm = TRUE) +
      geom_line(data = sims, aes(x, value, colour = engine), linewidth = 0.8, na.rm = TRUE) +
      scale_colour_manual(values = ENGINE_COLORS) +
      scale_fill_manual(values = ENGINE_COLORS, guide = "none") +
      coord_cartesian(xlim = range(x[idx]), ylim = if (is.null(rg)) NULL else rg$y) +
      labs(x = r$x$lab, y = names(ISOTOPES)[ISOTOPES == r$isotope]) +
      theme_app()
  })

  output$tab_perf <- renderTable({
    r <- res()
    no_results(r)
    ev <- r$eval_idx
    df <- data.frame(`Critère` = c(unname(CRIT_LABELS), "Temps [s]"), check.names = FALSE)
    for (e in r$engines) {
      b <- r$best[[e]]
      vals <- if (is.na(b$idx)) rep(NA_real_, length(CRIT_LABELS) + 1) else
        c(score_table(b$sim[ev], r$data$obs[ev], allow_NA = TRUE)[names(CRIT_LABELS)], r$time[[e]])
      df[[ENGINE_LABELS[[e]]]] <- formatC(vals, digits = 4, format = "fg")
    }
    df
  }, striped = TRUE, hover = TRUE, spacing = "s", width = "100%")

  output$tab_params <- renderTable({
    r <- res()
    no_results(r)
    ps <- sa_params(r)
    df <- data.frame(`Paramètre` = c(vapply(ps, function(p) p$name, character(1)),
                                     sprintf("Âge moyen du noyau [%s]", r$unit), "Meilleur run"),
                     check.names = FALSE)
    for (e in r$engines) {
      i <- r$best[[e]]$idx
      vals <- if (is.na(i)) rep(NA_real_, nrow(df)) else
        c(vapply(ps, function(p) p$value[i], numeric(1)), kernel_mean_age(best_kernel(r, e)) / r$prep$dt, i)
      df[[ENGINE_LABELS[[e]]]] <- formatC(vals, digits = 4, format = "fg")
    }
    df
  }, striped = TRUE, hover = TRUE, spacing = "s", width = "100%")

  kern_df <- reactive({
    r <- res()
    no_results(r)
    do.call(rbind, lapply(r$engines, function(e) {
      k <- best_kernel(r, e, decay = isTRUE(input$kern_decay))
      if (is.null(k)) return(NULL)
      k$density <- k$w / (k$hi - k$lo) * r$prep$dt
      k$mid <- (k$lo + k$hi) / 2 / r$prep$dt
      k$hi_u <- k$hi / r$prep$dt
      k$cum <- cumsum(k$w)
      k$engine <- ENGINE_LABELS[[e]]
      k$run <- r$best[[e]]$idx
      k
    }))
  })

  kern_xmax <- reactive({
    if (isTruthy(input$kern_max) && input$kern_max > 0) return(input$kern_max)
    k <- kern_df()
    mx <- vapply(split(k, k$engine), function(x) x$hi_u[which(x$cum >= 0.995 * sum(x$w))[1]], numeric(1))
    1.2 * max(mx, na.rm = TRUE)
  })

  kern_plot <- function(k, x, y, ylab) {
    p <- ggplot(k, aes(.data[[x]], .data[[y]], colour = engine)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = ENGINE_COLORS) +
      labs(x = sprintf("Âge [%s]", res()$unit), y = ylab) +
      theme_app()
    if (isTRUE(input$kern_log)) p <- p + scale_x_log10()
    p
  }

  output$plot_kern <- renderPlot({
    k <- kern_df()
    k <- k[k$mid > 0 & k$mid <= kern_xmax(), ]
    kern_plot(k, "mid", "density",
              sprintf(if (isTRUE(input$kern_decay)) "g(t)·e^(−λt) [1/%s]" else "g(t) [1/%s]", res()$unit))
  })

  output$plot_kern_cum <- renderPlot({
    k <- kern_df()
    k <- k[k$hi_u > 0 & k$hi_u <= kern_xmax(), ]
    kern_plot(k, "hi_u", "cum", "Fraction cumulée")
  })

  sa_engine_r <- reactive({
    r <- res()
    no_results(r)
    if (isTruthy(input$sa_engine) && input$sa_engine %in% r$engines) input$sa_engine else r$engines[1]
  })
  sa_pct_r <- reactive(or_default(input$sa_pct, 10))
  sa_col <- reactive(ENGINE_COLORS[[ENGINE_LABELS[[sa_engine_r()]]]])
  sa_df <- reactive({
    r <- res()
    no_results(r)
    sa_long(r, sa_pct_r(), sa_engine_r())
  })

  sa_names <- reactive({
    r <- res()
    no_results(r)
    vapply(sa_params(r), function(p) p$name, character(1))
  })

  output$ui_sa_xy <- renderUI({
    r <- res()
    req(r)
    nm <- sa_names()
    req(length(nm) >= 2)
    tagList(
      selectInput("sa_px", "x", nm, selected = or_default(isolate(input$sa_px), nm[1])),
      selectInput("sa_py", "y", nm, selected = or_default(isolate(input$sa_py), nm[2]))
    )
  })

  output$sa_dotty <- renderPlot({
    df <- sa_df()
    ggplot(df, aes(value, obj)) +
      geom_point(colour = sa_col(), alpha = 0.3, size = 0.9, na.rm = TRUE) +
      facet_wrap(~param, scales = "free_x", nrow = if (nlevels(df$param) <= 6) 1 else NULL) +
      labs(x = NULL, y = res()$crit) +
      theme_app()
  })

  output$sa_hist <- renderPlot({
    df <- sa_df()
    top <- df[df$behav, ]
    validate(need(nrow(top) > 0, "Aucun run retenu."))
    pri <- do.call(rbind, lapply(sa_params(res()), function(p) {
      if (anyNA(p$prior) || diff(p$prior) <= 0) return(NULL)
      data.frame(param = p$name, y = 1 / diff(p$prior))
    }))
    p <- ggplot(top, aes(value)) +
      geom_histogram(aes(y = after_stat(density)), bins = 25, fill = sa_col(), alpha = 0.7, na.rm = TRUE) +
      facet_wrap(~param, scales = "free") +
      labs(x = NULL, y = "Densité") +
      theme_app()
    if (!is.null(pri)) {
      pri$param <- factor(pri$param, levels = levels(df$param))
      p <- p + geom_hline(data = pri, aes(yintercept = y), linetype = "dashed", colour = "grey40")
    }
    p
  })

  output$sa_rsa <- renderPlot({
    r <- res()
    e <- sa_engine_r()
    d <- sa_df()
    d$decile <- crit_decile(r$obj[, e], r$crit)[d$run]
    d <- d[!is.na(d$decile) & !is.na(d$value), ]
    validate(need(nrow(d) > 0, "Aucun run valide."))
    ggplot(d, aes(value, colour = decile)) +
      stat_ecdf(geom = "step") +
      facet_wrap(~param, scales = "free_x") +
      scale_colour_viridis_d(name = "Décile", option = "D", end = 0.9) +
      labs(x = NULL, y = "Fonction de répartition") +
      theme_app() +
      theme(legend.title = element_text(size = 10))
  })

  output$sa_xy <- renderPlot({
    r <- res()
    e <- sa_engine_r()
    ps <- sa_params(r)
    nm <- sa_names()
    validate(need(length(nm) >= 2, "Un seul paramètre : pas d'interaction à représenter."))
    px <- or_default(input$sa_px, nm[1])
    py <- or_default(input$sa_py, nm[2])
    req(px %in% nm, py %in% nm)
    vx <- ps[[match(px, nm)]]$value
    vy <- ps[[match(py, nm)]]$value
    top <- top_runs(r$obj[, e], r$crit, sa_pct_r())
    d <- data.frame(x = vx, y = vy, obj = r$obj[, e])
    d <- d[order(if (r$crit == "RMSE") -d$obj else d$obj), ]
    ggplot(d, aes(x, y)) +
      geom_point(aes(colour = obj), size = 1.2, alpha = 0.8, na.rm = TRUE) +
      geom_point(data = data.frame(x = vx[top], y = vy[top]), shape = 1, size = 2, colour = "black") +
      scale_colour_viridis_c(name = r$crit, direction = if (r$crit == "RMSE") -1 else 1) +
      labs(x = px, y = py) +
      theme_app() +
      theme(legend.title = element_text(size = 10))
  })

  sa_tab <- reactive({
    r <- res()
    no_results(r)
    sa_table(r, sa_pct_r(), sa_engine_r())
  })

  output$sa_table <- renderTable(sa_tab(), striped = TRUE, hover = TRUE, spacing = "s", width = "100%")

  ex_scale <- function(input_id, slider_prefix, min_ok) {
    observeEvent(input[[input_id]], {
      mx <- input[[input_id]]
      req(isTruthy(mx), mx > min_ok)
      for (i in 1:2) updateSliderInput(session, paste0(slider_prefix, i), max = mx)
    })
  }
  ex_scale("ex_taumax", "ex_tau_", 0)

  last_dt_ex <- DT_YEARS[["month"]]
  observeEvent(input$dt, {
    dt <- dt_r()
    f <- last_dt_ex / dt
    last_dt_ex <<- dt
    unit <- unit_r()
    mx <- signif(or_default(input$ex_taumax, 600) * f, 3)
    updateNumericInput(session, "ex_taumax", value = mx, label = sprintf("τ max [%s]", unit))
    for (i in 1:2) {
      id <- paste0("ex_tau_", i)
      updateSliderInput(session, id, label = sprintf("τ [%s]", unit), max = mx,
                        value = signif(or_default(input[[id]], 60) * f, 4))
    }
  }, ignoreInit = TRUE)
  ex_scale("ex_nmax", "ex_n_", 1)
  ex_scale("ex_dpmax", "ex_dp_", 0)

  ex_par <- debounce(reactive({
    n <- as.integer(input$ex_n)
    comp <- lapply(seq_len(n), function(i) {
      g <- function(x) input[[paste0("ex_", x, "_", i)]]
      type <- g("type")
      req(type, g("tau"))
      list(type = type, tau = g("tau"),
           p2 = switch(type, EPM = , PEM = g("n"), DM = g("dp"), NA_real_))
    })
    list(engine = input$ex_engine,
         types = vapply(comp, function(cp) cp$type, character(1)),
         tau = vapply(comp, function(cp) cp$tau, numeric(1)),
         p2 = vapply(comp, function(cp) as.numeric(cp$p2), numeric(1)),
         ratio = if (n == 2) c(input$ex_f1, 1 - input$ex_f1) else 1)
  }), 250)

  ex_sim <- reactive({
    validate(need(!is.null(data_path()), "Importer un fichier de données dans l'onglet Configuration."))
    p <- ex_par()
    prep <- prep_r()
    req(all(p$tau > 0))
    simulate_run(p$engine, p$types, p$tau, p$p2, p$ratio, prep)
  })

  ex_kernel <- function(decay) {
    p <- ex_par()
    prep <- prep_r()
    mixture_kernel(p$engine, p$types, p$tau, p$p2, p$ratio, prep, decay = decay)
  }

  output$ex_vboxes <- renderUI({
    sim <- ex_sim()
    ev <- eval_idx_r()
    sc <- score_table(sim[ev], data_r()$obs[ev], allow_NA = TRUE)
    age <- kernel_mean_age(ex_kernel(FALSE)) / dt_r()
    vb <- function(title, value) bslib::value_box(title = title, value = value, class = "ex-vb")
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      vb("KGE", sprintf("%.3f", sc[["KGE"]])),
      vb("NSE", sprintf("%.3f", sc[["NSE"]])),
      vb("RMSE", sprintf("%.3g", sc[["RMSE"]])),
      vb("Âge moyen", sprintf("%.4g %s", age, unit_r()))
    )
  })

  output$ex_plot_sim <- renderPlot({
    sim <- ex_sim()
    d <- data_r()
    ax <- x_axis()
    e <- ex_par()$engine
    lab <- ENGINE_LABELS[[e]]
    df <- data.frame(x = ax$x, value = sim, serie = lab)
    r <- res()
    if (isTRUE(input$ex_show_best) && !is.null(r) && e %in% r$engines) {
      df <- rbind(data.frame(x = ax$x, value = r$best[[e]]$sim, serie = "Meilleur calage"), df)
    }
    cols <- c(stats::setNames(ENGINE_COLORS[[lab]], lab), "Meilleur calage" = "grey55", "Entrée (Cin)" = "grey80")
    p <- ggplot()
    if (isTRUE(input$ex_show_cin)) {
      p <- p + geom_line(data = data.frame(x = ax$x, value = d$Cin, serie = "Entrée (Cin)"),
                         aes(x, value, colour = serie), linewidth = 0.4)
    }
    p +
      eval_layers(ax$x, d$obs, eval_idx_r()) +
      geom_point(data = data.frame(x = ax$x, value = d$obs), aes(x, value), size = 1.3, colour = "grey15", na.rm = TRUE) +
      geom_line(data = df, aes(x, value, colour = serie), linewidth = 0.8, na.rm = TRUE) +
      scale_colour_manual(values = cols) +
      zoom_coord(ex_zoom()) +
      labs(x = ax$lab, y = names(ISOTOPES)[ISOTOPES == input$isotope]) +
      theme_app()
  })

  ex_kern_df <- reactive({
    k <- ex_kernel(isTRUE(input$ex_decay))
    dt <- dt_r()
    k$density <- k$w / (k$hi - k$lo) * dt
    k$mid <- (k$lo + k$hi) / 2 / dt
    k$hi_u <- k$hi / dt
    k$cum <- cumsum(k$w)
    xmax <- 1.2 * k$hi_u[which(k$cum >= 0.995 * sum(k$w))[1]]
    k[k$mid > 0 & k$mid <= xmax, ]
  })

  ex_kern_plot <- function(x, y, ylab) {
    k <- ex_kern_df()
    p <- ggplot(k, aes(.data[[x]], .data[[y]])) +
      geom_line(colour = ENGINE_COLORS[[ENGINE_LABELS[[ex_par()$engine]]]], linewidth = 0.9) +
      labs(x = sprintf("Âge [%s]", unit_r()), y = ylab) +
      theme_app()
    if (isTRUE(input$ex_log)) p <- p + scale_x_log10()
    p
  }

  output$ex_plot_kern <- renderPlot({
    ex_kern_plot("mid", "density",
                 sprintf(if (isTRUE(input$ex_decay)) "g(t)·e^(−λt) [1/%s]" else "g(t) [1/%s]", unit_r()))
  })
  output$ex_plot_cum <- renderPlot(ex_kern_plot("hi_u", "cum", "Fraction cumulée"))

  ex_info_cols <- function() {
    p <- ex_par()
    ev <- eval_idx_r()
    sc <- score_table(ex_sim()[ev], data_r()$obs[ev], allow_NA = TRUE)
    cols <- list(moteur = ENGINE_LABELS[[p$engine]], traceur = names(ISOTOPES)[ISOTOPES == input$isotope],
                 pas = unit_r())
    for (c in seq_along(p$types)) {
      cols[[paste0("modele_C", c)]] <- p$types[c]
      cols[[paste0("tau_C", c)]] <- p$tau[c]
      if (uses_p2(p$types[c])) cols[[paste0(p2_name(p$types[c]), "_C", c)]] <- p$p2[c]
    }
    if (length(p$types) == 2) cols$f1 <- p$ratio[1]
    c(cols, list(KGE = sc[["KGE"]], NSE = sc[["NSE"]], RMSE = sc[["RMSE"]],
                 age_moyen = kernel_mean_age(ex_kernel(FALSE)) / dt_r()))
  }

  output$ex_dl_sim <- downloadHandler(
    filename = function() "mitigre_exploration_simulation.csv",
    content = function(file) {
      d <- data_r()
      out <- data.frame(t = seq_len(nrow(d)), date = d$date, Cin = d$Cin, obs = d$obs, sim = ex_sim(),
                        evaluation = seq_len(nrow(d)) %in% eval_idx_r())
      write_csv(cbind(out, as.data.frame(ex_info_cols(), check.names = FALSE)), file)
    }
  )

  output$ex_dl_kern <- downloadHandler(
    filename = function() "mitigre_exploration_fonction_de_transfert.csv",
    content = function(file) {
      k <- ex_kernel(isTRUE(input$ex_decay))
      dt <- dt_r()
      out <- data.frame(age_debut = k$lo / dt, age_fin = k$hi / dt, poids = k$w,
                        densite = k$w / (k$hi - k$lo) * dt, cumul = cumsum(k$w))
      write_csv(cbind(out, as.data.frame(ex_info_cols(), check.names = FALSE)), file)
    }
  )

  observeEvent(input$ex_from_best, {
    r <- res()
    if (is.null(r)) {
      showNotification("Aucun calage disponible : lancer d'abord un calage.", type = "warning")
      return()
    }
    nc <- length(r$types)
    if (nc > 2) {
      showNotification("Le calage comporte plus de 2 composantes : non reproductible ici.", type = "warning")
      return()
    }
    e <- if (input$ex_engine %in% r$engines) input$ex_engine else r$engines[1]
    i <- r$best[[e]]$idx
    req(!is.na(i))
    tau_y <- r$par$tau[i, ] * r$prep$dt / dt_r()
    p2v <- r$par$p2[i, ]
    widen <- function(input_id, values, default) {
      mx <- or_default(input[[input_id]], default)
      v <- values[!is.na(values)]
      if (length(v) && max(v) > mx) {
        mx <- signif(1.5 * max(v), 2)
        updateNumericInput(session, input_id, value = mx)
      }
      mx
    }
    mx <- widen("ex_taumax", tau_y, 600)
    nmx <- widen("ex_nmax", p2v[r$types %in% c("EPM", "PEM")], 5)
    dpmx <- widen("ex_dpmax", p2v[r$types == "DM"], 2)
    updateRadioButtons(session, "ex_n", selected = as.character(nc))
    updateRadioButtons(session, "ex_engine", selected = e)
    for (c in seq_len(nc)) {
      updateSelectInput(session, paste0("ex_type_", c), selected = r$types[c])
      updateSliderInput(session, paste0("ex_tau_", c), value = tau_y[c], max = mx)
      p2 <- r$par$p2[i, c]
      if (r$types[c] %in% c("EPM", "PEM")) updateSliderInput(session, paste0("ex_n_", c), value = p2, max = nmx)
      if (r$types[c] == "DM") updateSliderInput(session, paste0("ex_dp_", c), value = p2, max = dpmx)
    }
    if (nc == 2) updateSliderInput(session, "ex_f1", value = r$par$ratio[i, 1])
  })

  write_csv <- function(df, file) utils::write.table(df, file, sep = ";", row.names = FALSE, fileEncoding = "UTF-8")

  output$dl_best <- downloadHandler(
    filename = function() "mitigre_meilleures_simulations.csv",
    content = function(file) {
      r <- res()
      req(r)
      out <- data.frame(t = seq_len(nrow(r$data)), date = r$data$date, Cin = r$data$Cin, obs = r$data$obs)
      for (e in r$engines) out[[paste0("sim_", ENGINE_LABELS[[e]])]] <- r$best[[e]]$sim
      write_csv(out, file)
    }
  )

  output$dl_runs <- downloadHandler(
    filename = function() "mitigre_runs_monte_carlo.csv",
    content = function(file) {
      r <- res()
      req(r)
      out <- data.frame(run = seq_len(r$n_run))
      for (p in sa_params(r)) out[[p$name]] <- p$value
      for (e in r$engines) out[[paste0(r$crit, "_", ENGINE_LABELS[[e]])]] <- r$obj[, e]
      write_csv(out, file)
    }
  )

  output$dl_kern <- downloadHandler(
    filename = function() "mitigre_fonctions_de_transfert.csv",
    content = function(file) {
      k <- kern_df()
      write_csv(k[, c("engine", "run", "lag", "lo", "hi", "w", "density", "cum")], file)
    }
  )

  output$dl_sens <- downloadHandler(
    filename = function() paste0("mitigre_sensibilite_", ENGINE_LABELS[[sa_engine_r()]], ".csv"),
    content = function(file) write_csv(sa_tab(), file)
  )

  observeEvent(input$notes_btn, {
    showModal(modalDialog(title = "Notes", notes_ui(), size = "xl", easyClose = TRUE,
                          footer = modalButton("Fermer")))
  })
}

shinyApp(ui, server)
