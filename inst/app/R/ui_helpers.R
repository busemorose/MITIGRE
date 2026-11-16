BANNER_COLOR <- "#1F4E6B"

info <- function(..., placement = "right") {
  bslib::tooltip(
    tags$span(class = "info-i", icon("circle-info")),
    ...,
    placement = placement
  )
}

label_info <- function(label, ...) tags$span(label, info(...))

app_css <- function() {
  tags$style(HTML(sprintf("
    .navbar { background-color: %1$s !important; border-bottom: 0; }
    .navbar .navbar-brand { color: #fff !important; font-weight: 700; letter-spacing: .08em; }
    .navbar .nav-link { color: rgba(255,255,255,.75) !important; }
    .navbar .nav-link.active, .navbar .nav-link:hover { color: #fff !important; }
    .navbar .nav-link.active { border-bottom: 2px solid #fff; }
    .navbar-toggler { border-color: rgba(255,255,255,.4); }
    .navbar-toggler-icon { filter: invert(1); }
    .banner-icon { font-size: 1.15rem; padding: .35rem .6rem; display: inline-block; text-decoration: none; }
    .banner-icon.export-off { color: rgba(255,255,255,.35); pointer-events: none; cursor: default; }
    .banner-icon.export-on { color: #5BD67A; }
    .banner-icon.export-on:hover { color: #8BF0A4; }
    .banner-icon.notes { color: rgba(255,255,255,.85); }
    .banner-icon.notes:hover { color: #fff; }
    .info-i { color: #8A9BA8; margin-left: .3rem; cursor: help; font-size: .85em; }
    .card-header { font-weight: 600; }
    .card-header .info-i { font-weight: 400; }
    .form-group { margin-bottom: .55rem; }
    .comp-card .card-header { background: #F2F6F9; padding: .35rem .75rem; }
    .comp-card .card-body { padding: .6rem .75rem .2rem; }
    .export-menu a { display: block; padding: .3rem 0; }
    .results-table table { width: 100%%; }
    .results-table td:not(:first-child), .results-table th:not(:first-child) { text-align: right; }
    .results-table tr.sep td { border-top: 2px solid #CBD5DD; }
    .modal-body h5 { margin-top: 1.1rem; color: %1$s; }
    .sa-bar .form-group { margin-bottom: 0; }
    .sa-bar .shiny-options-group { font-size: 1.05rem; }
    .bound-warn { font-size: .85rem; margin: .2rem 0 0; opacity: .95; }
    .ex-vb .value-box-title { font-size: .9rem; color: #6C7A86; margin-bottom: .1rem; }
    .ex-vb .value-box-value { font-size: 1.6rem; color: %1$s; }
    .ex-vb .value-box-area { padding: .6rem 1rem !important; }
  ", BANNER_COLOR)))
}

app_js <- function() {
  tags$script(HTML("
    Shiny.addCustomMessageHandler('export_state', function(ready) {
      var el = document.getElementById('export_icon');
      if (!el) return;
      el.classList.toggle('export-on', ready);
      el.classList.toggle('export-off', !ready);
    });
  "))
}

TIP_MODELS <- tags$div(
  style = "text-align:left",
  tags$b("EMM"), " : mélange exponentiel (aquifère libre homogène, captage sur toute l'épaisseur).", tags$br(),
  tags$b("EPM"), " : écoulement exponentiel puis piston ; n = ratio EPM + 1.", tags$br(),
  tags$b("PEM"), " : exponentiel partiel (seule la partie basse de l'aquifère est captée) ; n = ratio PEM + 1.", tags$br(),
  tags$b("DM"), " : dispersion (advection-dispersion 1D) ; DP = D/vx.", tags$br(),
  tags$b("PFM"), " : piston, sans mélange."
)
TIP_FORWARD <- "Convolution de toute la série d'entrée en une fois (schéma MITIGRE corrigé) : fonction de transfert évaluée au milieu de chaque pas et normalisée, décroissance ln2/T½ × dt, noyau prolongé jusqu'à convergence."
TIP_BACKWARD <- tags$div(
  style = "text-align:left",
  "Pour chaque date d'observation, remontée dans le passé par incréments d'âge jusqu'à ce que la sortie varie de moins de 10⁻⁶ ; intégrale exacte de g(t)·e^(−λt) sur chaque incrément.",
  tags$br(), tags$br(),
  tags$i("Méthode de TracerLPM : Jurgens, Böhlke & Eberts (2012), USGS Techniques and Methods 4-F3.")
)
TIP_TAU <- "Bornes du tirage Monte Carlo de l'âge moyen de la composante, en années. Si le meilleur run tombe en butée, élargir la plage."
TIP_SEED <- "Initialise le tirage aléatoire du Monte Carlo : une même graine redonne exactement les mêmes jeux de paramètres (résultats reproductibles). La changer pour un autre tirage."
TIP_WARMUP <- "Concentration attribuée à l'eau plus ancienne que le début de la série Cin (prolongée indéfiniment dans le passé). Dans TracerLPM, les cases vides avant la série valent 0 : choisir « Zéro » pour le reproduire."
TIP_QUAD <- "Forward : g évaluée au milieu de chaque pas puis normalisée (schéma MITIGRE), ou intégrée exactement sur chaque pas (comme TracerLPM). Avec l'intégrale exacte, Forward et Backward ne diffèrent plus que par le sens du calcul."
TIP_OBSPOS <- "Date d'une observation dans son pas de temps. Début : convention TracerLPM quand les dates sont au 1er jour du pas (ex. 2005, 2005-01-01) — l'eau d'âge 0 à 1 pas vient alors du pas précédent. Milieu : moyennes d'un pas. Fin : convention MITIGRE."
TIP_COMPAT <- "DM uniquement, pour retrouver un classeur TracerLPM : reproduit sa fonction DM (vérifié à 1e-9). TracerLPM arrondit la date d'échantillon au 1/24 d'an (1/12 si DP < 1), la recule d'un mois si elle tombe sur la grille et saute le premier incrément d'âge : la réponse est décalée d'environ un mois par rapport au Forward. Décoché (recommandé) : date exacte, intégrale complète."
TIP_CORES <- "Nombre de processus utilisés en parallèle pour le Monte Carlo. Non disponible sous Windows (1 cœur)."
TIP_TOL <- "Forward : le noyau est prolongé jusqu'à ce que la masse restante soit inférieure à la tolérance, dans la limite de K max pas."
TIP_DATA <- "Texte délimité (; , ou tabulation), décimale « . », colonnes obs et Cin, colonne date optionnelle (AAAA-MM-JJ, JJ/MM/AAAA ou année décimale)."

zoom_brush <- function(id) {
  brushOpts(id, direction = "xy", resetOnNew = TRUE, delay = 60000, delayType = "debounce")
}
