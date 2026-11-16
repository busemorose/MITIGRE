explore_comp_ui <- function(i) {
  id <- function(x) paste0("ex_", x, "_", i)
  type_id <- id("type")
  tagList(
    tags$div(class = "fw-semibold mt-2 mb-1", paste("Composante", i)),
    selectInput(type_id, NULL, MODEL_TYPES, selected = if (i == 1) "EMM" else "PFM"),
    sliderInput(id("tau"), "τ [mois]", min = 0.1, max = 600, value = if (i == 1) 60 else 240, step = 0.1),
    conditionalPanel(
      sprintf("input['%s'] == 'EPM' || input['%s'] == 'PEM'", type_id, type_id),
      sliderInput(id("n"), "n", min = 1, max = 5, value = 1.5, step = 0.01)
    ),
    conditionalPanel(
      sprintf("input['%s'] == 'DM'", type_id),
      sliderInput(id("dp"), "DP", min = 0.001, max = 2, value = 0.3, step = 0.001)
    )
  )
}

explore_page_ui <- function() {
  bslib::layout_sidebar(
    fillable = FALSE,
    sidebar = bslib::sidebar(
      width = 330,
      radioButtons("ex_engine", NULL, choiceNames = unname(ENGINE_LABELS[APP_ENGINES]),
                   choiceValues = APP_ENGINES, inline = TRUE),
      radioButtons("ex_n", "Composantes", c("1" = "1", "2 (mélange)" = "2"), inline = TRUE),
      label_info("Modèles", TIP_MODELS),
      explore_comp_ui(1),
      conditionalPanel(
        "input.ex_n == '2'",
        explore_comp_ui(2),
        sliderInput("ex_f1", label_info("Fraction f1", "Part de la composante 1 dans le mélange ; la composante 2 représente 1 − f1."),
                    min = 0, max = 1, value = 0.7, step = 0.01)
      ),
      tags$hr(class = "my-2"),
      label_info(tags$b("Échelle des curseurs"), "Valeur maximale de chaque curseur."),
      fluidRow(
        class = "mt-1",
        column(4, numericInput("ex_taumax", "τ max", value = 600, min = 0.1, step = 10)),
        column(4, numericInput("ex_nmax", "n max", value = 5, min = 1.1, step = 1)),
        column(4, numericInput("ex_dpmax", "DP max", value = 2, min = 0.01, step = 0.5))
      ),
      actionButton("ex_from_best", "Partir du meilleur calage", icon = icon("arrows-rotate"),
                   class = "btn-outline-primary btn-sm w-100")
    ),
    uiOutput("ex_vboxes"),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(
        class = "d-flex justify-content-between align-items-center",
        label_info("Simulation", "Critères calculés sur la période d'évaluation définie dans l'onglet Configuration, avec les mêmes options numériques. Sélectionner une zone pour zoomer, double-clic pour revenir à la vue entière."),
        gear(checkboxInput("ex_show_cin", "Afficher l'entrée Cin", value = FALSE),
             checkboxInput("ex_show_best", "Superposer le meilleur calage", value = FALSE),
             tags$hr(class = "my-2"),
             downloadLink("ex_dl_sim", tagList(icon("download"), " Exporter la simulation (CSV)")))
      ),
      plotOutput("ex_plot_sim", height = "380px",
                 brush = zoom_brush("ex_brush"),
                 dblclick = "ex_dblclick")
    ),
    bslib::card(
      full_screen = TRUE,
      bslib::card_header(
        class = "d-flex justify-content-between align-items-center",
        "Fonction de transfert",
        gear(checkboxInput("ex_decay", "Inclure la décroissance", value = FALSE),
             checkboxInput("ex_log", "Âge en échelle log", value = FALSE),
             tags$hr(class = "my-2"),
             downloadLink("ex_dl_kern", tagList(icon("download"), " Exporter la fonction de transfert (CSV)")))
      ),
      bslib::layout_columns(
        col_widths = c(6, 6),
        plotOutput("ex_plot_kern", height = "280px"),
        plotOutput("ex_plot_cum", height = "280px")
      )
    )
  )
}
