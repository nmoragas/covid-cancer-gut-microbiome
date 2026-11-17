# =============================================================================
# 04 - LASSO prediction model of COVID-19 severity
# -----------------------------------------------------------------------------
#   Rscript analysis/04_lasso_model.R          # Fig 8a
#
# Ambulatory (asymptomatic/mild) vs hospitalised (moderate + severe) COVID-19,
# baseline samples. Requires 03_differential_abundance.R severity.
#
# Candidate predictors
#   - 90 species differentially abundant between ambulatory and hospitalised
#     patients (CLR values)
#   - clinical covariates: age, vaccination, any comorbidity, smoking,
#     tumour stage, tumour status, antibiotics (last 3 months)
#
# Procedure
#   1. Prevalence filter, zero replacement and CLR computed on all baseline
#      samples (COVID-19 and non-COVID-19); COVID-19 patients are then selected.
#   2. Missing vaccination / smoking status imputed with missForest.
#   3. LASSO logistic regression (glmnet, alpha = 1, no intercept), 5-fold CV
#      optimising AUC, repeated with seeds 1..100; coefficients at lambda.1se.
#   4. Consensus model: variables selected in >= 90/100 repetitions,
#      coefficient = mean across repetitions.
#   5. Linear predictor = X %*% beta, min-max normalised; threshold 0.45.
#      AUC with DeLong 95% CI.
#
# Published model: age + 9 species (AUC 0.89, 95% CI 0.83-0.96).
# =============================================================================

source("config/config.R")
source("config/comparisons.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(missForest)
  library(pROC)
  library(plotROC)
  library(verification)
  library(openxlsx)
})

out_dir <- out_dir_for("04_prediction_model", "severity")

COVARIATES   <- c("age", "vaccination", "any_comorbidity_data", "smoke",
                  "tumor_stage", "tumor_status", "abx")
N_REPEATS    <- 100
NFOLDS       <- 5        # < 10 observations per fold with 10 folds
MIN_SELECTED <- 90
THRESHOLD    <- 0.45     # chosen from the distribution of normalised predictions

# -----------------------------------------------------------------------------
# 1. Abundance data: all baseline samples
# -----------------------------------------------------------------------------
ps   <- load_marato()
ps   <- prune_samples(get_metadata(ps)$tipus_kit %in% "id_kit_basal", ps)
meta <- get_metadata(ps)
clr_abund <- clr_transform(replace_zeros(filter_prevalence(abundance_table(ps))))

# -----------------------------------------------------------------------------
# 2. Covariates (imputation + encoding on all baseline samples)
# -----------------------------------------------------------------------------
meta$vaccination[meta$vaccination %in% "unknown"] <- NA
meta$smoke[meta$smoke %in% "unknown"] <- NA
meta$abx <- any_antibiotic(meta)

cov_df <- meta[, COVARIATES]
cov_df[] <- lapply(cov_df, function(x) if (is.character(x) || is.logical(x)) factor(x) else x)
cov_df$age <- as.numeric(cov_df$age)
set.seed(SEED)
cov_imp <- missForest(cov_df)$ximp
meta$vaccination <- as.character(cov_imp$vaccination)
meta$smoke       <- as.character(cov_imp$smoke)

X_cov <- as.matrix(encode_covariates(meta, COVARIATES))

# -----------------------------------------------------------------------------
# 3. Severity dataset
# -----------------------------------------------------------------------------
cmp   <- get_comparison("severity")
keep  <- meta$COHORT %in% "CASE" & !is.na(meta$covid_severity_modif) &
         meta$covid_severity_modif != "uninfected"
samples <- rownames(meta)[keep]
y <- cmp$group(meta[samples, ])
print(table(y))

da <- read.xlsx(file.path(PATHS$results, "03_differential_abundance", "severity",
                          "severity_hospitalized_DA_combined.xlsx"))
species <- intersect(da$taxon, colnames(clr_abund))
message(length(species), " of ", nrow(da), " severity species available")

X <- cbind(X_cov[samples, ], clr_abund[samples, species])

# -----------------------------------------------------------------------------
# 4. Repeated cross-validated LASSO
# -----------------------------------------------------------------------------
fits <- lapply(seq_len(N_REPEATS), function(s) {
  set.seed(s)
  cv.glmnet(X, y, family = "binomial", alpha = 1, intercept = FALSE,
            type.measure = "auc", nfolds = NFOLDS)
})
saveRDS(fits, file.path(out_dir, "lasso_cv_fits.rds"))

coefs <- sapply(fits, function(f) as.numeric(coef(f, s = f$lambda.1se)))
rownames(coefs) <- rownames(coef(fits[[1]]))
coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]

n_selected <- rowSums(coefs != 0)
keep_vars  <- names(n_selected)[n_selected >= MIN_SELECTED]
model_coef <- data.frame(ccm = rowMeans(coefs[keep_vars, , drop = FALSE]),
                         n_selected = n_selected[keep_vars], row.names = keep_vars)
print(model_coef)

saveRDS(as.matrix(model_coef[, "ccm", drop = FALSE]),
        file.path(out_dir, "model_severity_coefficients.rds"))
saveRDS(THRESHOLD, file.path(out_dir, "model_severity_threshold.rds"))
write.xlsx(model_coef, file.path(out_dir, "model_severity_coefficients.xlsx"), rowNames = TRUE)

# -----------------------------------------------------------------------------
# 5. Performance (training cohort)
# -----------------------------------------------------------------------------
pred   <- as.numeric(X[, keep_vars, drop = FALSE] %*% model_coef$ccm)
pred_n <- (pred - min(pred)) / (max(pred) - min(pred))
y_bin  <- as.numeric(y == "hospitalized")

roc_obj <- roc(y_bin, pred, quiet = TRUE)
ci_auc  <- ci.auc(roc_obj)
p_auc   <- verification::roc.area(y_bin, pred_n)$p.value

cm <- table(predicted = factor(as.numeric(pred_n > THRESHOLD), 0:1), true = factor(y_bin, 0:1))
perf <- data.frame(
  AUC = as.numeric(auc(roc_obj)), AUC_low = ci_auc[1], AUC_high = ci_auc[3], p_value = p_auc,
  threshold   = THRESHOLD,
  sensitivity = cm[2, 2] / sum(cm[, 2]),
  specificity = cm[1, 1] / sum(cm[, 1]),
  accuracy    = sum(diag(cm)) / sum(cm)
)
print(cm); print(perf)
write.xlsx(perf, file.path(out_dir, "model_severity_performance.xlsx"))

pred_df <- data.frame(sample = samples, group = y, pred = pred, pred_n = pred_n)
write.xlsx(pred_df, file.path(out_dir, "model_severity_predictions.xlsx"))

# Density of normalised predictions
p_dens <- ggplot(pred_df, aes(pred_n, colour = group, fill = group)) +
  geom_density(linewidth = 0.8, alpha = 0.2) +
  scale_fill_manual(values = cmp$colours, labels = cmp$labels, name = NULL) +
  scale_colour_manual(values = cmp$colours, labels = cmp$labels, name = NULL) +
  labs(x = "Model prediction", y = "Density",
       caption = paste0("t-test p = ", signif(t.test(pred ~ group, data = pred_df)$p.value, 3))) +
  theme_article() + theme(legend.position = "bottom")
save_plot(p_dens, file.path(out_dir, "Fig8a_density_severity.pdf"), 12, 12)

# ROC curve
roc_df <- data.frame(prediction = pred_n, labels = y_bin, id = "hospitalized")
p_roc <- ggplot(roc_df, aes(m = prediction, d = labels, colour = id)) +
  geom_roc(labels = FALSE, pointsize = 0, n.cuts = 0) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.1) +
  scale_colour_manual(values = c(hospitalized = "salmon4")) +
  annotate("text", x = 0.65, y = 0.20, size = 5,
           label = sprintf("AUC: %.2f 95%% CI (%.2f - %.2f)", perf$AUC, perf$AUC_low, perf$AUC_high)) +
  annotate("text", x = 0.54, y = 0.10, size = 5,
           label = paste("p-value:", format(p_auc, scientific = TRUE, digits = 3))) +
  labs(x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)") +
  theme_article()
save_plot(p_roc, file.path(out_dir, "Fig8a_ROC_severity.pdf"), 12, 12)
