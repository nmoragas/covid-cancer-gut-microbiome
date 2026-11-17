# =============================================================================
# Shared helper functions
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(zCompositions)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Command-line argument with default
# -----------------------------------------------------------------------------
get_arg <- function(default) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) args[1] else default
}

out_dir_for <- function(step, comparison = NULL) {
  d <- if (is.null(comparison)) file.path(PATHS$results, step)
       else file.path(PATHS$results, step, comparison)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------

#' Load the MARATÓ phyloseq object, drop excluded samples and recode
#' string "NA" values to real NA.
load_marato <- function(path = PATHS$phyloseq, exclude = EXCLUDED_SAMPLES) {
  ps <- readRDS(path)
  ps <- prune_samples(!(sample_names(ps) %in% exclude), ps)
  sd <- data.frame(sample_data(ps), check.names = FALSE)
  sd[] <- lapply(sd, function(x) { if (is.character(x)) x[x == "NA"] <- NA; x })
  sample_data(ps) <- sample_data(sd)
  ps
}

get_metadata <- function(ps) data.frame(sample_data(ps), check.names = FALSE)

#' Subset a phyloseq object to the samples and grouping of a comparison.
#' Returns list(ps, meta, group).
subset_comparison <- function(ps, cmp) {
  meta  <- get_metadata(ps)
  keep  <- cmp$keep(meta)
  keep[is.na(keep)] <- FALSE
  meta  <- meta[keep, , drop = FALSE]
  group <- cmp$group(meta)
  ok    <- !is.na(group)
  meta  <- meta[ok, , drop = FALSE]
  group <- droplevels(group[ok])
  names(group) <- rownames(meta)
  ps    <- prune_samples(rownames(meta), ps)
  ps    <- prune_taxa(taxa_sums(ps) > 0, ps)
  message(sprintf("[%s] %d samples: %s", cmp$name, length(group),
                  paste(names(table(group)), table(group), sep = "=", collapse = ", ")))
  list(ps = ps, meta = meta, group = group)
}

# -----------------------------------------------------------------------------
# Abundance tables
# -----------------------------------------------------------------------------

#' Clean UHGG species names: "s__1XD42-69_sp014287635" -> "1XD42.69_sp014287635"
#' (same result as the data.frame(check.names = TRUE) conversion used originally)
clean_taxa_names <- function(x) gsub(".*__", "", gsub("-", ".", x))

#' Samples x species abundance table with cleaned names.
abundance_table <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  colnames(m) <- clean_taxa_names(colnames(m))
  as.data.frame(m, check.names = FALSE)
}

#' Prevalence/abundance filter (species in columns).
#' Keeps species with relative abundance > relabundance_filter in at least
#' nsamples_filter samples.
filter_matrix <- function(count_matrix, relabundance_filter = 0.001, nsamples_filter = 1) {
  relab <- t(apply(count_matrix, 1, function(x) x / sum(x)))
  keep  <- apply(relab, 2, function(x) length(which(x > relabundance_filter)))
  keep  <- names(keep[keep >= nsamples_filter])
  count_matrix[, keep, drop = FALSE]
}

filter_prevalence <- function(abund,
                              relab = FILTER$relabundance,
                              prevalence = FILTER$prevalence) {
  filt <- filter_matrix(abund, relabundance_filter = relab,
                        nsamples_filter = ceiling(nrow(abund) * prevalence))
  message(sprintf("Species before / after filter: %d / %d", ncol(abund), ncol(filt)))
  filt
}

#' Bayesian-multiplicative zero replacement (zCompositions, SQ prior).
replace_zeros <- function(abund, ...) {
  as.data.frame(cmultRepl(abund, output = "p-counts", method = "SQ", ...),
                check.names = FALSE)
}

#' Centred log-ratio transform (samples in rows). Equivalent to compositions::clr().
clr_transform <- function(x) {
  lx <- log(as.matrix(x))
  sweep(lx, 1, rowMeans(lx), "-")
}

# -----------------------------------------------------------------------------
# Covariate encoding for LASSO models
# -----------------------------------------------------------------------------
#' Numeric, centred encoding of clinical covariates.
#'   age          -> (age / 10) centred at 0
#'   binary vars  -> -0.5 / +0.5
#'   ordinal vars -> -1 / 0 / +1
encode_covariates <- function(meta, vars) {
  enc <- list(
    age                  = function(x) { x <- as.numeric(x) / 10; x - mean(x, na.rm = TRUE) },
    vaccination          = function(x) as.numeric(factor(x, levels = c("no", "yes"))) - 1.5,
    any_comorbidity_data = function(x) as.numeric(factor(x, levels = c("No", "Yes"))) - 1.5,
    smoke                = function(x) as.numeric(factor(x, levels = c("non smoker", "ex smoker", "active smoker"))) - 2,
    tumor_stage          = function(x) as.numeric(factor(x, levels = c("localized", "loco-regional", "metastatic"))) - 2,
    tumor_status         = function(x) as.numeric(factor(x, levels = c("Resmision", "Active"))) - 1.5,
    abx                  = function(x) as.numeric(factor(x, levels = c("FALSE", "TRUE"))) - 1.5
  )
  out <- sapply(vars, function(v) {
    if (!v %in% names(enc)) stop("No encoding defined for covariate '", v, "'")
    enc[[v]](meta[[v]])
  })
  out <- as.data.frame(out)
  rownames(out) <- rownames(meta)
  out
}

#' Any antibiotic (oral or injected) in the previous 3 months.
any_antibiotic <- function(meta) {
  factor(meta$oral_antibiotic != "No" | meta$antibiotic_injection != "No")
}

# -----------------------------------------------------------------------------
# Plot helpers
# -----------------------------------------------------------------------------
theme_article <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(axis.text  = element_text(size = base_size),
          axis.title = element_text(size = base_size),
          plot.title = element_text(hjust = 0),
          legend.position = "none",
          axis.line  = element_line(colour = "black", linewidth = 0.5))
}

save_plot <- function(p, file, width = 12, height = 12) {
  ggsave(file, plot = p, width = width, height = height, units = "cm")
  invisible(file)
}

#' All pairwise combinations of factor levels (for ggpubr comparisons).
pairwise_levels <- function(f) combn(levels(f), 2, simplify = FALSE)
