# =============================================================================
# 07 - Faecal cytokines
# -----------------------------------------------------------------------------
#   Rscript analysis/07_cytokines.R
#
# Outputs
#   - Pre-processed cytokine table (cytokines_preprocessed.rds)
#   - Fig 4d   : cytokines by COVID-19 severity (ambulatory vs hospitalised)
#   - Fig 4e   : Spearman correlation between severity-associated species
#                (DA, n = 90) and cytokines in hospitalised patients
#   - Supp Fig 3c : cytokines by COVID-19 status x mortality (4 groups)
#
# Cytokines: IL-1b, IL-1a, IL-17A, CCL7, CCL2, IFN-g, TNF-a, IL-6
#
# Pre-processing
#   1. Values below the limit of quantification ("<LOQ") -> LOQ / 2
#   2. Main and pilot batches merged (no batch effect, see PCA check)
#   3. Technical replicates averaged per sample
#   4. Samples with >= 50% missing cytokines removed
#   5. Remaining missing values imputed with kNN (VIM)
#   6. log10(x + 1) transformation and z-score scaling
# =============================================================================

source("config/config.R")
source("config/comparisons.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(VIM)
  library(ggpubr)
  library(Hmisc)
  library(pheatmap)
})

out_dir <- out_dir_for("07_cytokines")

CYTOKINES <- c("IL_1b", "IL_1a", "IL_17A", "CCL7", "CCL2",
               "IFNg_3rd_gen", "TNF_a_2nd_gen", "IL_6_2nd_gen")

# Use nominal p-values (as in the published Fig 4e) or BH-adjusted ones
HEATMAP_SIGNIFICANCE <- "p"     # "p" or "fdr"

# =============================================================================
# 1. Pre-processing
# =============================================================================
main  <- read.xlsx(PATHS$cytokines)
pilot <- read.xlsx(PATHS$cytokines_pilot)
colnames(main)  <- c("ID", CYTOKINES)
colnames(pilot) <- c("ID", CYTOKINES)
main$batch  <- "main"
pilot$batch <- "pilot"
raw <- rbind(main, pilot)

loq_to_numeric <- function(x) {
  x <- gsub(",", ".", as.character(x))
  below <- grepl("<", x)
  v <- suppressWarnings(as.numeric(gsub("<", "", x)))
  v[below] <- v[below] / 2
  v
}
raw[CYTOKINES] <- lapply(raw[CYTOKINES], loq_to_numeric)

# Mean of replicates
cyto <- raw %>%
  group_by(ID) %>%
  summarise(across(all_of(CYTOKINES), ~ mean(.x, na.rm = TRUE)),
            batch = first(batch), .groups = "drop") %>%
  mutate(across(all_of(CYTOKINES), ~ ifelse(is.nan(.x), NA, .x)))

# Remove samples with >= 50% missing values and impute the rest (kNN)
cyto <- as.data.frame(cyto[rowMeans(is.na(cyto[CYTOKINES])) < 0.5, ])
imputed <- kNN(cyto[, c("ID", CYTOKINES)], variable = CYTOKINES, imp_var = FALSE)
cyto[CYTOKINES] <- imputed[CYTOKINES]

# log10 + z-score
cyto[CYTOKINES] <- scale(log10(cyto[CYTOKINES] + 1))

# Batch check (PCA)
pc <- prcomp(cyto[CYTOKINES])
p_batch <- ggplot(data.frame(pc$x[, 1:2], batch = cyto$batch), aes(PC1, PC2, colour = batch)) +
  geom_point() + theme_classic() + ggtitle("Cytokines PCA - batch")
save_plot(p_batch, file.path(out_dir, "cytokines_batch_PCA.pdf"), 14, 12)

# Merge with sample metadata
ps   <- load_marato()
meta <- get_metadata(ps)
meta$ID <- rownames(meta)
cyto_all <- merge(cyto, meta, by = "ID")
cyto_all <- cyto_all[!duplicated(cyto_all$ID), ]
rownames(cyto_all) <- cyto_all$ID
saveRDS(cyto_all, file.path(out_dir, "cytokines_preprocessed.rds"))
message(nrow(cyto_all), " samples with cytokine and microbiome data")

# Wilcoxon test per cytokine
test_cytokines <- function(df, group_var) {
  bind_rows(lapply(CYTOKINES, function(cy) {
    r <- compare_means(as.formula(paste(cy, "~", group_var)), data = df, method = "wilcox.test")
    r$cytokine <- cy
    r
  })) %>% mutate(p_BH_all = p.adjust(p, method = "BH"))
}

boxplot_cytokines <- function(df, group_var, colours, labels) {
  long <- pivot_longer(df, all_of(CYTOKINES), names_to = "Cytokine", values_to = "Level")
  long$Cytokine <- factor(long$Cytokine, levels = CYTOKINES)
  ggboxplot(long, x = "Cytokine", y = "Level", color = group_var, palette = colours) +
    scale_colour_manual(values = colours, labels = labels, name = NULL) +
    labs(x = NULL, y = "Scaled log10 level") +
    theme(legend.position = "top") +
    rotate_x_text(angle = 45)
}

# =============================================================================
# 2. Severity (Fig 4d)
# =============================================================================
sev_cmp <- get_comparison("severity")
cyto_sev <- cyto_all[cyto_all$tipus_kit %in% "id_kit_basal" & cyto_all$COHORT %in% "CASE" &
                     !is.na(cyto_all$covid_severity_modif), ]
cyto_sev$severity <- sev_cmp$group(cyto_sev)
print(table(cyto_sev$severity))

write.xlsx(test_cytokines(cyto_sev, "severity"),
           file.path(out_dir, "cytokines_wilcoxon_severity.xlsx"))
p4d <- boxplot_cytokines(cyto_sev, "severity", sev_cmp$colours, sev_cmp$labels)
save_plot(p4d, file.path(out_dir, "Fig4d_cytokines_severity.pdf"), 15, 12)

# =============================================================================
# 3. Correlation severity-associated species vs cytokines (Fig 4e)
#    Hospitalised patients; 90 DA species from the severity comparison.
# =============================================================================
da_sev <- read.xlsx(file.path(PATHS$results, "03_differential_abundance", "severity",
                              "severity_hospitalized_DA_combined.xlsx"))
abund  <- abundance_table(ps)
hosp   <- rownames(cyto_sev)[cyto_sev$severity == "hospitalized"]
sp     <- intersect(da_sev$taxon, colnames(abund))
message(length(sp), " of ", nrow(da_sev), " severity species found")

cor_res <- rcorr(as.matrix(abund[hosp, sp]), as.matrix(cyto_sev[hosp, CYTOKINES]),
                 type = "spearman")
r_mat <- cor_res$r[sp, CYTOKINES]
p_mat <- cor_res$P[sp, CYTOKINES]
q_mat <- matrix(p.adjust(p_mat, method = "BH"), nrow = nrow(p_mat), dimnames = dimnames(p_mat))

write.xlsx(list(rho = as.data.frame(r_mat), p = as.data.frame(p_mat), fdr = as.data.frame(q_mat)),
           file.path(out_dir, "spearman_severity_species_cytokines.xlsx"), rowNames = TRUE)

sig_mat <- if (HEATMAP_SIGNIFICANCE == "fdr") q_mat else p_mat
rows_sig <- rownames(sig_mat)[apply(sig_mat < 0.05, 1, any, na.rm = TRUE)]
stars <- ifelse(sig_mat < 0.001, "***", ifelse(sig_mat < 0.01, "**",
                ifelse(sig_mat < 0.05, "*", "")))

if (length(rows_sig) > 1) {
  pdf(file.path(out_dir, "Fig4e_heatmap_species_cytokines.pdf"), width = 5, height = 0.18 * length(rows_sig) + 3)
  pheatmap(r_mat[rows_sig, ], display_numbers = stars[rows_sig, ],
           cluster_rows = FALSE, cluster_cols = TRUE, treeheight_col = 0,
           scale = "column", cellwidth = 11, cellheight = 10, fontsize = 9,
           fontsize_number = 7, number_color = "black", border_color = "grey93",
           color = colorRampPalette(c("lightblue4", "#FFFFFF", "salmon4"))(10),
           angle_col = 45,
           main = sprintf("Cytokines - Spearman correlation\n(* %s < 0.05)",
                          ifelse(HEATMAP_SIGNIFICANCE == "fdr", "FDR", "p")))
  dev.off()
}

# =============================================================================
# 4. COVID-19 status x mortality (Supp Fig 3c)
# =============================================================================
mort_cmp <- get_comparison("mortality_4groups")
cyto_mort <- cyto_all[cyto_all$tipus_kit %in% "id_kit_basal" & !is.na(cyto_all$outcome), ]
cyto_mort$outcome_group <- mort_cmp$group(cyto_mort)
print(table(cyto_mort$outcome_group))

write.xlsx(test_cytokines(cyto_mort, "outcome_group"),
           file.path(out_dir, "cytokines_wilcoxon_mortality_4groups.xlsx"))
pS3c <- boxplot_cytokines(cyto_mort, "outcome_group", mort_cmp$colours, mort_cmp$labels)
save_plot(pS3c, file.path(out_dir, "SuppFig3c_cytokines_mortality.pdf"), 19, 12)
