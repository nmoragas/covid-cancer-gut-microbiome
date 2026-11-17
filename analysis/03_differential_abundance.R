# =============================================================================
# 03 - Differential abundance (ANCOM-BC + LinDA)
# -----------------------------------------------------------------------------
# Generic script applied to every comparison in config/comparisons.R.
# Species detected as differentially abundant by BOTH methods are considered
# robust. For comparisons with more than two groups, every level is tested
# against the reference (first) level.
#
#   Rscript analysis/03_differential_abundance.R case_control     # Fig 2b, Supp Table 1
#   Rscript analysis/03_differential_abundance.R timepoint        # Fig 2d, Supp Table 3
#   Rscript analysis/03_differential_abundance.R severity         # Fig 2f, Supp Table 4
#   Rscript analysis/03_differential_abundance.R mortality        # Fig 6b, Supp Table 11
#   Rscript analysis/03_differential_abundance.R mortality_covid  # COVID-19 stratum
#   Rscript analysis/03_differential_abundance.R mortality_noncovid
#
# Input  : species filtered (0.1% relative abundance in >= 5% samples) and
#          zero-replaced (zCompositions SQ) abundance table.
# Output : results/03_differential_abundance/<comparison>/
#            <comparison>_<level>_DA_combined.xlsx   (robust species)
#            <comparison>_<level>_barplot_logFC.pdf
#            ancombc.rds / linda.rds                 (full model objects)
# =============================================================================

source("config/config.R")
source("config/comparisons.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(MicrobiomeStat)
  library(dplyr)
  library(openxlsx)
  library(ggvenn)
})

comparison <- get_arg("case_control")
cmp     <- get_comparison(comparison)
out_dir <- out_dir_for("03_differential_abundance", comparison)
ALPHA   <- if (!is.null(cmp$da_alpha)) cmp$da_alpha else 0.05

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
ps  <- load_marato()
dat <- subset_comparison(ps, cmp)
meta <- dat$meta; group <- dat$group

abund_filt <- filter_prevalence(abundance_table(dat$ps))
zeros      <- replace_zeros(abund_filt)

meta_da <- data.frame(group = group, row.names = rownames(meta))
for (v in cmp$da_covariates) meta_da[[v]] <- as.factor(meta[[v]])
cc <- complete.cases(meta_da)
if (any(!cc)) message("Dropping ", sum(!cc), " samples with missing covariates")
meta_da <- meta_da[cc, , drop = FALSE]
zeros   <- zeros[rownames(meta_da), ]

rhs <- paste(c("group", cmp$da_covariates), collapse = " + ")
message("Model: ~ ", rhs)

# -----------------------------------------------------------------------------
# ANCOM-BC
# -----------------------------------------------------------------------------
ps_da <- phyloseq(otu_table(t(as.matrix(zeros)), taxa_are_rows = TRUE),
                  sample_data(meta_da))

set.seed(SEED)
ancom <- ancombc(data = ps_da, formula = rhs, p_adj_method = "fdr",
                 lib_cut = 0, group = "group", struc_zero = TRUE, neg_lb = TRUE,
                 tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = ALPHA,
                 global = nlevels(group) > 2)
saveRDS(ancom, file.path(out_dir, "ancombc.rds"))

# -----------------------------------------------------------------------------
# LinDA
# -----------------------------------------------------------------------------
linda_out <- linda(feature.dat = as.data.frame(ceiling(t(as.matrix(zeros)))),
                   meta.dat = meta_da, formula = paste("~", rhs),
                   alpha = ALPHA, prev.filter = 0, mean.abund.filter = 0)
saveRDS(linda_out, file.path(out_dir, "linda.rds"))

# -----------------------------------------------------------------------------
# Combine both methods, one contrast at a time
# -----------------------------------------------------------------------------
ancom_column <- function(tab, coef) as.numeric(as.matrix(tab[[coef]]))

for (lv in levels(group)[-1]) {
  coef <- paste0("group", lv)
  contrast <- paste0(lv, "_vs_", levels(group)[1])

  res <- ancom$res
  ancom_df <- data.frame(
    taxon = res$lfc$taxon,
    lfc   = ancom_column(res$lfc, coef),
    se    = ancom_column(res$se, coef),
    p_val = ancom_column(res$p_val, coef),
    q_val = ancom_column(res$q_val, coef),
    diff  = as.logical(ancom_column(res$diff_abn, coef))
  )
  linda_df <- linda_out$output[[coef]]
  linda_df$taxon <- rownames(linda_df)

  sig_ancom <- ancom_df$taxon[ancom_df$diff %in% TRUE]
  sig_linda <- linda_df$taxon[linda_df$reject %in% TRUE]
  robust    <- intersect(sig_ancom, sig_linda)
  message(sprintf("%s: ANCOM-BC %d, LinDA %d, both %d",
                  contrast, length(sig_ancom), length(sig_linda), length(robust)))

  # Venn diagram of both methods
  pv <- ggvenn(list(LinDA = sig_linda, `ANCOM-BC` = sig_ancom),
               fill_color = c("salmon3", "salmon4"), stroke_size = 0.5,
               set_name_size = 5, text_size = 6)
  save_plot(pv, file.path(out_dir, paste0(comparison, "_", lv, "_venn.pdf")), 12, 12)

  if (length(robust) == 0) next

  # NOTE: ANCOM-BC 'lfc' is on the natural-log scale and LinDA 'log2FoldChange'
  # on log2. They are averaged here exactly as in the published analysis.
  combined <- merge(ancom_df[ancom_df$taxon %in% robust, ],
                    linda_df[linda_df$taxon %in% robust,
                             c("taxon", "log2FoldChange", "lfcSE", "pvalue", "padj")],
                    by = "taxon") %>%
    mutate(mean_logFC = (log2FoldChange + lfc) / 2,
           mean_SE    = (lfcSE + se) / 2,
           group      = ifelse(mean_logFC > 0, "Enriched", "Depleted")) %>%
    arrange(desc(abs(mean_logFC)))

  write.xlsx(combined, file.path(out_dir, paste0(comparison, "_", lv, "_DA_combined.xlsx")))

  p <- ggplot(combined, aes(x = reorder(taxon, mean_logFC), y = mean_logFC, fill = group)) +
    geom_col(colour = "black", alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_logFC - mean_SE, ymax = mean_logFC + mean_SE), width = 0.2) +
    scale_fill_manual(values = c(Depleted = "lightblue4", Enriched = "salmon4"), name = NULL) +
    coord_flip() +
    labs(x = "Species", y = "Log fold change (mean ANCOM-BC / LinDA)",
         title = paste("Differential abundance:", contrast)) +
    theme_minimal(base_size = BASE_SIZE) +
    theme(legend.position = "bottom", axis.line = element_line(colour = "black"))
  save_plot(p, file.path(out_dir, paste0(comparison, "_", lv, "_barplot_logFC.pdf")),
            20, max(10, 0.6 * nrow(combined) + 4))
}

message("Results written to ", out_dir)
