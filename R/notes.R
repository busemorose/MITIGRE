notes_ui <- function() {
  tags$div(
    h5("Forward et Backward", style = "margin-top:0"),
    p("Les deux moteurs calculent la même intégrale de convolution, ",
      tags$code("Cout(t) = ∫ Cin(t − a) · g(a) · e^(−λa) da"), ", dans deux sens :"),
    tags$ul(
      tags$li(tags$b("Forward"), " : toute la série de sortie est obtenue d'un coup par convolution de l'entrée ",
              "avec la fonction de transfert (schéma de MITIGRE, corrigé). g est évaluée au milieu de chaque pas puis normalisée."),
      tags$li(tags$b("Backward"), " : pour chaque date d'observation, on remonte le temps par incréments d'âge ",
              "(1/12 an ; 1/24 an pour le DM si DP ≥ 1 ; pas des données si τ > 100 a) en additionnant ",
              "Cin × ∫ g(a)e^(−λa) da sur l'incrément, jusqu'à ce que l'âge dépasse τ et que l'incrément relatif ",
              "soit < 10⁻⁶ (10⁻⁷ pour le DM). C'est la boucle des fonctions EMM_VB, EPM_VB, PEM_VB, DM_VB et PFM_VB de TracerLPM.")
    ),
    p("Les deux moteurs évaluent exactement les mêmes jeux de paramètres (Monte Carlo à graine commune). ",
      "À paramètres identiques, les écarts proviennent de :"),
    tags$ul(
      tags$li("la quadrature : point milieu (Forward) contre intégrale exacte par incrément (Backward). Négligeable ",
              "pour EMM et DM dès que τ dépasse quelques pas ; jusqu'à ~1–2 % pour EPM et PEM, dont la distribution ",
              "démarre brutalement au milieu d'un pas ;"),
      tags$li("le critère d'arrêt du Backward, qui néglige la queue de la distribution : 10⁻⁵ à 6·10⁻⁴ en relatif.")
    ),
    p("Le sens du calcul lui-même n'introduit aucun écart : une version forward à intégrale exacte reproduit le ",
      "Backward à ce critère d'arrêt près (tests/test_engines.R)."),

    h5("Corrections apportées à MITIGRE (moteur Forward)"),
    tags$ol(
      tags$li("Décroissance par pas : λ = ln2 / T½ × dt (MITIGRE applique ln2/12,32 par pas, juste seulement en pas annuel)."),
      tags$li("Décroissance du ¹⁴C (T½ = 5730 a), absente de MITIGRE."),
      tags$li("Décroissance évaluée au milieu du pas, comme g (MITIGRE : fin de pas)."),
      tags$li("Seuil d'apparition de l'EPM testé au milieu du pas, comme pour le PEM."),
      tags$li("Noyau prolongé jusqu'à convergence au lieu d'être tronqué à max_t et renormalisé ; âges antérieurs à la série = valeur du warmup (MITIGRE : 0)."),
      tags$li("PFM décalé de ⌈τ − o⌉ pas, cohérent avec la convolution (MITIGRE : ⌈τ⌉, un pas de trop)."),
      tags$li("Plus de dépendance à data.table.")
    ),

    h5("Paramètres et conventions"),
    tags$ul(
      tags$li("τ : âge moyen de la composante, dans l'unité du pas de temps des données (jours, mois…) ; bornes de tirage uniformes saisies dans la carte Modèle. Un meilleur run en butée de plage est signalé."),
      tags$li("EPM, PEM : n = ratio + 1 (n = 1 : EMM). DM : DP = D/vx. PEM paramétré par l'âge moyen de la ",
              "portion captée, comme TracerLPM (τ_aquifère = τ / (1 + ln n))."),
      tags$li("Plusieurs composantes : ratios tirés puis normalisés à 1 ; sortie = somme pondérée."),
      tags$li("Position de l'observation dans son pas : le retard 0 couvre les âges [0, o] pas, le retard i les âges ",
              "[o + i − 1, o + i]. o = 1 (fin, convention MITIGRE) ou 0,5 (milieu, moyennes)."),
      tags$li("Eau plus ancienne que le début de la série Cin : concentration choisie dans « Concentration avant la série », prolongée indéfiniment (Forward et Backward). TracerLPM : cases vides avant la série = 0 (choisir « Zéro »), sinon première valeur prolongée.")
    ),

    h5("Retrouver les valeurs du classeur TracerLPM"),
    tags$ul(
      tags$li("Observation datée : « Début » si les dates d'échantillon de TracerLPM sont au 1er jour du pas (2005, 2005-01-01). En pas annuel, « Fin » décale tout d'un an."),
      tags$li("Concentration avant la série : « Zéro » si les cases d'entrée de TracerLPM avant ta série sont vides."),
      tags$li("DM : cocher « Backward : compatibilité TracerLPM ». La fonction DM du classeur (add-in XLL, identique à DM_VB) arrondit la date d'échantillon au 1/24 d'an (1/12 si DP < 1), la recule d'un mois si elle tombe sur la grille, compte les âges depuis cette date et saute le premier incrément. Vérifié : écart < 1e-9 avec un classeur TracerLPM. EMM, EPM, PEM et PFM utilisent la date exacte et sont reproduits sans option (1e-14)."),
      tags$li("Mêmes paramètres : τ en années dans TracerLPM ; DP identique ; n = ratio + 1.")
    ),

    h5("Écarts assumés au VBA de TracerLPM"),
    tags$ul(
      tags$li("Critère d'arrêt en valeur absolue : le VBA exige une sortie positive, ce qui bloque les traceurs négatifs (δ¹⁸O)."),
      tags$li("DM intégré exactement (loi inverse gaussienne) : le VBA utilise un Simpson adaptatif peu précis, compte les âges depuis une date arrondie à la grille et saute le premier incrément."),
      tags$li("Non portés : ³He tritiogénique, ³H₀, ³H/³H₀, pondération du ¹⁴C par le DIC, ⁴He, zone non saturée.")
    ),

    h5("Validation"),
    p("Backward identique à 10⁻¹⁴ près à un portage littéral de EMM_VB, EPM_VB, PEM_VB et PFM_VB ; ",
      "Forward identique à MITIGRE d'origine (à 10⁻¹¹) quand les corrections sont sans effet ; masses unitaires ",
      "et âges moyens égaux à τ. Script : tests/test_engines.R."),

    h5("Source"),
    p("Jurgens, B.C., Böhlke, J.K., Eberts, S.M. (2012). TracerLPM (Version 1): An Excel® workbook for ",
      "interpreting groundwater age distributions from environmental tracer data. U.S. Geological Survey ",
      "Techniques and Methods 4-F3, 60 p. ",
      tags$a(href = "https://pubs.usgs.gov/publication/tm4F3", target = "_blank", "pubs.usgs.gov/publication/tm4F3"))
  )
}
