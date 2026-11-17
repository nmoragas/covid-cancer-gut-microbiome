# =============================================================================
# 02 - Alpha and beta diversity
# -----------------------------------------------------------------------------
# Generic script applied to every clinical comparison defined in
# config/comparisons.R.
#
#   Rscript analysis/02_diversity.R case_control      # Fig 2a, Supp Fig 1a
#   Rscript analysis/02_diversity.R timepoint         # Fig 2c, Supp Fig 1b, Supp Table 2
#   Rscript analysis/02_diversity.R severity          # Fig 2e, Supp Fig 1e
#   Rscript analysis/02_diversity.R mortality         # Fig 6a, Supp Table 12
#   Rscript analysis/02_diversity.R mortality_4groups # Supp Fig 3a
#
# Alpha diversity : Shannon (rarefied counts) and Chao1 (unrarefied, rounded)
#                   Wilcoxon (2 groups) / Kruskal-Wallis + pairwise Wilcoxon
# Beta diversity  : Aitchison distance (Euclidean on CLR-transformed data);
#                   PCA of CLR data (= Aitchison PCoA); PERMANOVA (adonis2),
#                   pairwise PERMANOVA (FDR) and dispersion (betadisper).
# =============================================================================

source("config/config.R")
source("config/comparisons.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(vegan)
  library(fossil)
  library(ggpubr)
  library(RVAideMemoire)
  library(openxlsx)
})

comparison <- get_arg("case_control")
cmp     <- get_comparison(comparison)
out_dir <- out_dir_for("02_diversity", comparison)

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
ps  <- load_marato()
dat <- subset_comparison(ps, cmp)
ps <- dat$ps; meta <- dat$meta; group <- dat$group
names(cmp$colours) <- levels(group)
names(cmp$labels)  <- levels(group)

abund <- abundance_table(ps)

# -----------------------------------------------------------------------------
# Alpha diversity
# -----------------------------------------------------------------------------
ps_rare <- rarefy_even_depth(ps, rngseed = SEED, sample.size = min(sample_sums(ps)),
                             replace = TRUE, verbose = FALSE)
rare_tab <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) rare_tab <- t(rare_tab)

alpha <- data.frame(
  sample  = rownames(meta),
  group   = group,
  Shannon = vegan::diversity(rare_tab[rownames(meta), ], index = "shannon"),
  Chao1   = apply(round(abund[rownames(meta), ]), 1, fossil::chao1)
)
write.xlsx(alpha, file.path(out_dir, "alpha_diversity_values.xlsx"))

alpha_tests <- list()
for (idx in c("Shannon", "Chao1")) {
  f <- as.formula(paste(idx, "~ group"))
  global <- if (nlevels(group) == 2) wilcox.test(f, data = alpha) else kruskal.test(f, data = alpha)
  pairwise <- compare_means(f, data = alpha, method = "wilcox.test", p.adjust.method = "BH")
  alpha_tests[[paste0(idx, "_global")]]   <- data.frame(index = idx, method = global$method,
                                                        p = global$p.value)
  alpha_tests[[paste0(idx, "_pairwise")]] <- as.data.frame(pairwise)

  p <- ggboxplot(alpha, x = "group", y = idx, color = "group",
                 palette = unname(cmp$colours), title = paste(idx, "diversity index")) +
    scale_x_discrete(labels = cmp$labels) +
    labs(x = NULL, y = paste(idx, "index")) +
    theme_article()
  if (nlevels(group) == 2) {
    p <- p + stat_compare_means(method = "wilcox.test")
  } else {
    p <- p + stat_compare_means(comparisons = pairwise_levels(group), method = "wilcox.test") +
             stat_compare_means(method = "kruskal.test", label.y.npc = "top")
  }
  save_plot(p, file.path(out_dir, paste0("alpha_", idx, "_", comparison, ".pdf")), 15, 12)
}
write.xlsx(alpha_tests, file.path(out_dir, "alpha_diversity_tests.xlsx"))

# -----------------------------------------------------------------------------
# Beta diversity (Aitchison)
# -----------------------------------------------------------------------------
abund_filt <- filter_prevalence(abund)
zeros      <- replace_zeros(abund_filt)
clr_abund  <- clr_transform(zeros)
clr_dist   <- dist(clr_abund, method = "euclidean")

pca <- prcomp(clr_abund, center = TRUE, scale. = FALSE)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group)

set.seed(SEED)
permanova <- adonis2(clr_dist ~ group, data = data.frame(group = group),
                     permutations = PERMUTATIONS)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = group)) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.68) +
  scale_colour_manual(values = cmp$colours, labels = cmp$labels, name = NULL) +
  labs(title = "clr-abundance PCA (species)",
       subtitle = sprintf("PERMANOVA p = %.4f", permanova$`Pr(>F)`[1]),
       x = paste0("PC1 (", var_expl[1], "%)"), y = paste0("PC2 (", var_expl[2], "%)")) +
  theme_article() + theme(legend.position = "right")
save_plot(p_pca, file.path(out_dir, paste0("beta_PCA_", comparison, ".pdf")), 18, 12)

beta_tests <- list(permanova = as.data.frame(permanova))

# Pairwise PERMANOVA (more than two groups)
if (nlevels(group) > 2) {
  set.seed(SEED)
  pw <- pairwise.perm.manova(clr_dist, group, nperm = PERMUTATIONS, p.method = "fdr")
  beta_tests$pairwise_permanova_fdr <- as.data.frame(pw$p.value)
}

# Homogeneity of multivariate dispersions
disper <- betadisper(clr_dist, group = group)
set.seed(SEED)
beta_tests$betadisper <- as.data.frame(permutest(disper, permutations = PERMUTATIONS)$tab)

# Mortality: additive and interaction models with COVID-19 status (Supp Table 12)
if (isTRUE(cmp$permanova_covid_interaction)) {
  df <- data.frame(outcome = group,
                   COHORT  = factor(meta$COHORT, levels = c("CONTROL", "CASE")))
  df$group4 <- interaction(df$COHORT, df$outcome, sep = "_")
  set.seed(SEED)
  beta_tests$permanova_additive <- as.data.frame(
    adonis2(clr_dist ~ COHORT + outcome, data = df, permutations = PERMUTATIONS, by = "terms"))
  set.seed(SEED)
  beta_tests$permanova_interaction <- as.data.frame(
    adonis2(clr_dist ~ COHORT * outcome, data = df, permutations = PERMUTATIONS, by = "terms"))
  set.seed(SEED)
  beta_tests$betadisper_4groups <- as.data.frame(
    permutest(betadisper(clr_dist, df$group4), permutations = PERMUTATIONS)$tab)
}

write.xlsx(beta_tests, file.path(out_dir, "beta_diversity_tests.xlsx"), rowNames = TRUE)
message("Results written to ", out_dir)
