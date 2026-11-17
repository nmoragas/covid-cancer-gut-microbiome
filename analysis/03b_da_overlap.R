# =============================================================================
# 03b - Overlap of differentially abundant species across contrasts
# -----------------------------------------------------------------------------
# Compares the robust species (ANCOM-BC & LinDA) of:
#   - COVID-19 vs non-COVID-19            (case_control, CASE)
#   - active infection vs non-COVID-19    (timepoint, Active)
#   - hospitalised vs ambulatory          (severity, hospitalized)
# and defines the candidate permissive / protective taxa (Fig 2g) that are
# followed longitudinally in 06_follow_up.R.
#
# Requires 03_differential_abundance.R to have been run for the three
# comparisons.
#   Rscript analysis/03b_da_overlap.R
# =============================================================================

source("config/config.R")
source("R/utils.R")
suppressPackageStartupMessages({ library(openxlsx); library(ggvenn) })

out_dir <- out_dir_for("03_differential_abundance", "overlap")
da_file <- function(cmp, lv) file.path(PATHS$results, "03_differential_abundance", cmp,
                                       paste0(cmp, "_", lv, "_DA_combined.xlsx"))

da <- list(
  COVID19  = read.xlsx(da_file("case_control", "CASE")),
  Active   = read.xlsx(da_file("timepoint", "Active")),
  Severity = read.xlsx(da_file("severity", "hospitalized"))
)

enriched <- lapply(da, function(x) x$taxon[x$group == "Enriched"])
depleted <- lapply(da, function(x) x$taxon[x$group == "Depleted"])

all_taxa <- unique(unlist(lapply(da, `[[`, "taxon")))
overlap <- data.frame(
  taxon = all_taxa,
  sapply(names(da), function(n) {
    d <- da[[n]]; d$group[match(all_taxa, d$taxon)]
  }),
  check.names = FALSE
)
overlap$n_contrasts <- rowSums(!is.na(overlap[, names(da)]))
overlap <- overlap[order(-overlap$n_contrasts, overlap$taxon), ]
write.xlsx(overlap, file.path(out_dir, "DA_overlap_table.xlsx"))

p_enr <- ggvenn(enriched, fill_color = c("salmon4", "lightcoral", "salmon2"),
                set_name_size = 4, text_size = 5) + ggtitle("Enriched")
p_dep <- ggvenn(depleted, fill_color = c("lightblue4", "darkslategrey", "cadetblue3"),
                set_name_size = 4, text_size = 5) + ggtitle("Depleted")
save_plot(p_enr, file.path(out_dir, "venn_enriched.pdf"), 14, 14)
save_plot(p_dep, file.path(out_dir, "venn_depleted.pdf"), 14, 14)

# Candidate permissive taxa: enriched in all three contrasts, or in
# COVID-19 & active infection, or in COVID-19 & severity
core3      <- Reduce(intersect, enriched)
covid_act  <- setdiff(intersect(enriched$COVID19, enriched$Active), core3)
covid_sev  <- setdiff(intersect(enriched$COVID19, enriched$Severity), core3)

candidates <- list(
  enriched_all_three          = core3,
  enriched_covid_and_active   = covid_act,
  enriched_covid_and_severity = covid_sev,
  depleted_active             = depleted$Active
)
print(candidates)
saveRDS(candidates, file.path(out_dir, "candidate_taxa.rds"))
