# =============================================================================
# Global configuration: paths, shared parameters and colour palettes
# -----------------------------------------------------------------------------
# All analysis scripts must be run from the repository root, e.g.
#   Rscript analysis/02_diversity.R severity
# Edit the paths below to point to your local copies of the input data
# (see data/README.md).
# =============================================================================

PATHS <- list(
  # Species-level phyloseq object (UHGG v2.0.2, Bracken, genome-length normalised)
  phyloseq        = "data/phyloseq_obj_F_uhgg_sel.rds",
  # External validation cohort (CRA003945), processed with the same pipeline
  phyloseq_valid  = "data/phyloseq_obj_CRA003945.rds",
  # GCAT Microbiome cohort (healthy reference population)
  phyloseq_gcat   = "data/GCAT_phyloseqobject_BEC.rds",
  metadata_gcat   = "data/metadata_GCAT.rds",
  # Faecal cytokines (main batch + pilot batch)
  cytokines       = "data/cytokine_viral_load_covid.xlsx",
  cytokines_pilot = "data/Cytokines data_pilot study.xlsx",
  # GMHI species lists (Gupta et al., Nat Commun 2020)
  gmhi_MH         = "data/gmhi/MH_species.txt",
  gmhi_MN         = "data/gmhi/MN_species.txt",
  # Output folder (git-ignored)
  results         = "results"
)

# Sample excluded from all analyses: patient returned to country of origin,
# no follow-up available (covid_severity = NA).
EXCLUDED_SAMPLES <- c("F_900")

SEED <- 2024

# Prevalence / abundance filter applied before compositional analyses:
# keep species with relative abundance > 0.1% in at least 5% of samples.
FILTER <- list(relabundance = 0.001, prevalence = 0.05)

PERMUTATIONS <- 9999

# Font sizes shared by all figures
BASE_SIZE <- 15

dir.create(PATHS$results, showWarnings = FALSE, recursive = TRUE)
