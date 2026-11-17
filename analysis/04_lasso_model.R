# =============================================================================
# 04 - LASSO prediction model (microbiome + clinical covariates)
# -----------------------------------------------------------------------------
#   Rscript analysis/04_lasso_model.R severity        # Fig 8a
#
# !! RECONSTRUCTED SCRIPT – TO BE CHECKED AGAINST THE ORIGINAL
#    (4_prediction_model/c_severity/5_lasso_severity_ambul_hospi_c11.Rmd).
#    The procedure follows the consensus-LASSO approach used across the
#    project and produces the coefficient table read by 05_external_validation.R.
#
# Procedure
#   1. Candidate predictors: robust DA species of the comparison (CLR values)
#      + clinical covariates (config: model_covariates), centred/scaled.
#   2. LASSO logistic regression (glmnet, alpha = 1), 10-fold CV optimising AUC,
#      repeated with 100 different seeds; coefficients taken at lambda.1se.
#   3. Consensus model: variables selected in >= 90/100 repetitions;
#      coefficient = mean over repetitions. The intercept is dropped.
#   4. Linear predictor = X %*% beta, min-max normalised to [0, 1];
#      ROC/AUC (DeLong 95% CI), sensitivity, specificity and accuracy at the
#      Youden-optimal threshold.
# =============================================================================

source("config/config.R")
source("config/comparisons.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(openxlsx)
})

comparison <- get_arg("severity")
cmp     <- get_comparison(comparison)
out_dir <- out_dir_for("04_prediction_model", comparison)

N_REPEATS     <- 100
MIN_SELECTED  <- 90

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
ps  <- load_marato()
dat <- subset_comparison(ps, cmp)
meta <- dat$meta; group <- dat$group
positive_lv <- levels(group)[2]

clr_abund <- clr_transform(replace_zeros(filter_prevalence(abundance_table(dat$ps))))

da <- read.xlsx(file.path(PATHS$results, "03_differential_abundance", comparison,
                          paste0(comparison, "_", positive_lv, "_DA_combined.xlsx")))
species <- intersect(da$taxon, colnames(clr_abund))
message(length(species), " candidate species")

X <- clr_abund[, species, drop = FALSE]
if (length(cmp$model_covariates)) {
  X <- cbind(as.matrix(encode_covariates(meta, cmp$model_covariates)), X)
}
X <- X[complete.cases(X), , drop = FALSE]
y <- as.numeric(group[rownames(X)] == positive_lv)

# -----------------------------------------------------------------------------
# Repeated cross-validated LASSO
# -----------------------------------------------------------------------------
fits <- lapply(seq_len(N_REPEATS), function(s) {
  set.seed(s)
  cv.glmnet(x = X, y = y, family = "binomial", alpha = 1,
            type.measure = "auc", nfolds = 10)
})
saveRDS(fits, file.path(out_dir, "lasso_cv_fits.rds"))

coefs <- sapply(fits, function(f) as.numeric(coef(f, s = "lambda.1se")))
rownames(coefs) <- rownames(coef(fits[[1]]))
coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]

n_selected <- rowSums(coefs != 0)
keep <- names(n_selected)[n_selected >= MIN_SELECTED]
model_coef <- data.frame(coef = rowMeans(coefs[keep, , drop = FALSE]),
                         n_selected = n_selected[keep], row.names = keep)
print(model_coef)

saveRDS(model_coef[, "coef", drop = FALSE],
        file.path(out_dir, paste0("model_", comparison, "_coefficients.rds")))
write.xlsx(model_coef, file.path(out_dir, paste0("model_", comparison, "_coefficients.xlsx")),
           rowNames = TRUE)

# -----------------------------------------------------------------------------
# Performance in the training cohort
# -----------------------------------------------------------------------------
pred   <- as.numeric(X[, keep, drop = FALSE] %*% model_coef$coef)
pred_n <- (pred - min(pred)) / (max(pred) - min(pred))

roc_obj <- roc(y, pred_n, quiet = TRUE)
ci_auc  <- ci.auc(roc_obj)
best    <- coords(roc_obj, "best", best.method = "youden", ret = "threshold")$threshold[1]

pred_bin <- as.numeric(pred_n >= best)
cm <- table(predicted = factor(pred_bin, 0:1), true = factor(y, 0:1))
perf <- data.frame(
  AUC = as.numeric(auc(roc_obj)), AUC_low = ci_auc[1], AUC_high = ci_auc[3],
  threshold   = best,
  sensitivity = cm[2, 2] / sum(cm[, 2]),
  specificity = cm[1, 1] / sum(cm[, 1]),
  accuracy    = sum(diag(cm)) / sum(cm)
)
print(perf)
write.xlsx(perf, file.path(out_dir, paste0("model_", comparison, "_performance.xlsx")))
saveRDS(best, file.path(out_dir, paste0("model_", comparison, "_threshold.rds")))

p <- ggroc(roc_obj, colour = "salmon4", linewidth = 1) +
  geom_abline(slope = 1, intercept = 1, linewidth = 0.2) +
  annotate("text", x = 0.35, y = 0.15, size = 5,
           label = sprintf("AUC %.2f (95%% CI %.2f-%.2f)", perf$AUC, perf$AUC_low, perf$AUC_high)) +
  labs(title = paste("Training -", cmp$description),
       x = "Specificity", y = "Sensitivity") +
  theme_article()
save_plot(p, file.path(out_dir, paste0("ROC_", comparison, "_training.pdf")), 12, 12)
