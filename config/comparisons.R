# =============================================================================
# Definition of the clinical comparisons analysed in the article
# -----------------------------------------------------------------------------
# Each comparison specifies:
#   keep        function(meta) -> logical: which samples enter the comparison
#   group       function(meta) -> factor: grouping variable (first level = reference)
#   labels      display labels for each level (same order as the factor levels)
#   colours     colours for each level
#   da_covariates  covariates added to the differential abundance models
#
# Level names are kept syntactic (no spaces/accents) because they become
# coefficient names in ANCOM-BC / LinDA.
#
# NOTE: the sample filters reproduce those used in the original per-comparison
# scripts ("tipus_kit" for case-control, "basal_line" for severity/mortality).
# =============================================================================

recode_outcome <- function(x) {
  x <- as.character(x)
  x[x == "Viu"]  <- "Alive"
  x[x == "Mort"] <- "Deceased"
  factor(x, levels = c("Alive", "Deceased"))
}

COMPARISONS <- list(

  # --- COVID-19 vs non-COVID-19 (baseline samples) -------------------------
  case_control = list(
    description = "COVID-19 vs non-COVID-19 cancer controls (baseline samples)",
    keep    = function(m) m$tipus_kit %in% "id_kit_basal",
    group   = function(m) factor(as.character(m$COHORT), levels = c("CONTROL", "CASE")),
    labels  = c("Non-COVID-19", "COVID-19"),
    colours = c("lightblue4", "salmon4"),
    da_covariates = NULL
  ),

  # --- Time from diagnosis (all samples) -----------------------------------
  # Status_.PN.: "No covid" | "Resuelto" (>=90 d) | "En resolución" (8-89 d) | "Activo" (0-7 d)
  timepoint = list(
    description = "Non-COVID-19 vs recovered vs recovering vs active infection (all samples)",
    keep    = function(m) !is.na(m$Status_.PN.),
    group   = function(m) {
      x <- as.character(m$Status_.PN.)
      x <- c("No covid" = "NonCOVID", "Resuelto" = "Recovered",
             "En resolución" = "Recovering", "Activo" = "Active")[x]
      factor(x, levels = c("NonCOVID", "Recovered", "Recovering", "Active"))
    },
    labels  = c("Non-COVID-19", "Recovered", "Recovering", "Active"),
    colours = c("lightblue4", "darkslategrey", "lightcoral", "salmon4"),
    da_covariates = NULL
  ),

  # --- Severity: ambulatory vs hospitalised (COVID-19, baseline) -----------
  severity = list(
    description = "Asymptomatic/mild (ambulatory) vs hospitalised (moderate + severe) COVID-19",
    keep    = function(m) m$basal_line %in% "Yes" & m$COHORT %in% "CASE" &
                          !is.na(m$covid_severity_modif) & m$covid_severity_modif != "NA",
    group   = function(m) {
      x <- ifelse(m$covid_severity_modif %in% c("moderate_disease", "severe_disease"),
                  "hospitalized", "ambulatory")
      factor(x, levels = c("ambulatory", "hospitalized"))
    },
    labels  = c("Asymptomatic/Mild", "Hospitalised"),
    colours = c("lightcoral", "salmon4"),
    da_covariates    = NULL
  ),

  # --- Severity with three levels (exploratory) -----------------------------
  severity_3levels = list(
    description = "Asymptomatic/mild vs moderate vs severe COVID-19",
    keep    = function(m) m$basal_line %in% "Yes" & m$COHORT %in% "CASE" &
                          !is.na(m$covid_severity_modif) & m$covid_severity_modif != "NA",
    group   = function(m) {
      x <- c("asymptomatic/mild disease" = "mild", "moderate_disease" = "moderate",
             "severe_disease" = "severe")[as.character(m$covid_severity_modif)]
      factor(x, levels = c("mild", "moderate", "severe"))
    },
    labels  = c("Asymptomatic/Mild", "Moderate", "Severe"),
    colours = c("lightcoral", "salmon2", "salmon4"),
    da_covariates = NULL
  ),

  # --- All-cause mortality (baseline samples, whole cohort) -----------------
  mortality = list(
    description = "Alive vs deceased (baseline samples, COVID-19 and non-COVID-19)",
    keep    = function(m) m$basal_line %in% "Yes" & !is.na(m$outcome),
    group   = function(m) recode_outcome(m$outcome),
    labels  = c("Alive", "Deceased"),
    colours = c("lightblue4", "salmon4"),
    da_covariates = NULL,
    # Run additive / interaction PERMANOVA with COVID-19 status (Supp. Table 12)
    permanova_covid_interaction = TRUE
  ),

  # --- Mortality stratified by COVID-19 status (4 groups) -------------------
  mortality_4groups = list(
    description = "Non-COVID-19/COVID-19 x alive/deceased (baseline samples)",
    keep    = function(m) m$basal_line %in% "Yes" & !is.na(m$outcome),
    group   = function(m) {
      x <- paste(as.character(m$COHORT), as.character(recode_outcome(m$outcome)), sep = "_")
      factor(x, levels = c("CONTROL_Alive", "CONTROL_Deceased", "CASE_Alive", "CASE_Deceased"))
    },
    labels  = c("Non-COVID-19 alive", "Non-COVID-19 deceased", "COVID-19 alive", "COVID-19 deceased"),
    colours = c("lightblue4", "salmon4", "darkslategrey", "lightcoral"),
    da_covariates = NULL
  ),

  # --- Mortality within COVID-19 patients -----------------------------------
  mortality_covid = list(
    description = "Alive vs deceased among COVID-19 patients (baseline samples)",
    keep    = function(m) m$basal_line %in% "Yes" & m$COHORT %in% "CASE" & !is.na(m$outcome),
    group   = function(m) recode_outcome(m$outcome),
    labels  = c("Alive", "Deceased"),
    colours = c("darkslategrey", "salmon4"),
    da_covariates = c("tumor_stage")
  ),

  # --- Mortality within non-COVID-19 patients -------------------------------
  mortality_noncovid = list(
    description = "Alive vs deceased among non-COVID-19 patients (baseline samples)",
    keep    = function(m) m$basal_line %in% "Yes" & m$COHORT %in% "CONTROL" & !is.na(m$outcome),
    group   = function(m) recode_outcome(m$outcome),
    labels  = c("Alive", "Deceased"),
    colours = c("lightblue4", "salmon4"),
    da_covariates = c("tumor_stage")
  )
)

get_comparison <- function(name) {
  if (!name %in% names(COMPARISONS)) {
    stop("Unknown comparison '", name, "'. Available: ",
         paste(names(COMPARISONS), collapse = ", "), call. = FALSE)
  }
  cmp <- COMPARISONS[[name]]
  cmp$name <- name
  cmp
}
