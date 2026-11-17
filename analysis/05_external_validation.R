# =============================================================================
# 05 - External validation of the severity model (cohort CRA003945)
# -----------------------------------------------------------------------------
#   Rscript analysis/05_external_validation.R            # Fig 8b
#
# The CRA003945 raw reads (NGDC Genome Sequence Archive) were processed with
# preprocessing/shotgun_pipeline.sh, exactly as the MARATÓ cohort.
#
# Groups mapped to the MARATÓ definition:
#   "Taken from asyptomatic infected people"     -> ambulatory   (n = 13)
#   "Taken from people being cured of COVID-19"  -> hospitalized (n = 24)
#
# Notes
# - The prevalence filter is NOT applied here because it would remove some of
#   the model species; zero replacement + CLR are applied to all species.
# - Age is not available for this cohort: the (centred) age term is set to 0,
#   i.e. every subject is given the mean age, so prediction relies on the
#   microbiome terms only.
# - The classification threshold is the one selected in the training cohort.
# =============================================================================

source("config/config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(pROC)
  library(plotROC)
  library(verification)
  library(openxlsx)
})

out_dir <- out_dir_for("05_external_validation")

model_dir <- file.path(PATHS$results, "04_prediction_model", "severity")
ccm   <- readRDS(file.path(model_dir, "model_severity_coefficients.rds"))
beta  <- setNames(ccm[, 1], rownames(ccm))
# Threshold selected in the training cohort (0.45 in the published model)
TH <- if (file.exists(file.path(model_dir, "model_severity_threshold.rds")))
        readRDS(file.path(model_dir, "model_severity_threshold.rds")) else 0.45

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
ps    <- readRDS(PATHS$phyloseq_valid)      # genome-length normalised abundances
meta  <- get_metadata(ps)

groups <- c("Taken from asyptomatic infected people"    = "ambulatory",
            "Taken from people being cured of COVID-19" = "hospitalized")
meta  <- meta[meta$Public.description %in% names(groups), , drop = FALSE]
meta$severity <- factor(groups[meta$Public.description], levels = c("ambulatory", "hospitalized"))
print(table(meta$severity))

abund <- abundance_table(ps)
clr_abund <- clr_transform(replace_zeros(abund))     # cmultRepl may drop samples
samples   <- intersect(rownames(meta), rownames(clr_abund))
meta      <- meta[samples, , drop = FALSE]

missing <- setdiff(names(beta), c("age", colnames(clr_abund)))
if (length(missing)) stop("Model species not found in validation data: ",
                          paste(missing, collapse = ", "))

X <- cbind(age = 0, clr_abund[samples, , drop = FALSE])[, names(beta), drop = FALSE]

# -----------------------------------------------------------------------------
# Prediction
# -----------------------------------------------------------------------------
pred   <- as.numeric(X %*% beta)
pred_n <- (pred - min(pred)) / (max(pred) - min(pred))
y      <- as.numeric(meta$severity == "hospitalized")

res <- data.frame(sample = samples, severity = meta$severity, pred = pred, pred_n = pred_n)
write.xlsx(res, file.path(out_dir, "validation_predictions.xlsx"))

roc_obj <- roc(y, pred_n, quiet = TRUE)
ci_auc  <- ci.auc(roc_obj)
p_auc   <- verification::roc.area(y, pred_n)$p.value

pred_bin <- as.numeric(pred_n >= TH)
cm <- table(predicted = factor(pred_bin, 0:1), true = factor(y, 0:1))
perf <- data.frame(
  n_ambulatory = sum(y == 0), n_hospitalized = sum(y == 1),
  AUC = as.numeric(auc(roc_obj)), AUC_low = ci_auc[1], AUC_high = ci_auc[3], p_value = p_auc,
  threshold   = TH,
  sensitivity = cm[2, 2] / sum(cm[, 2]),
  specificity = cm[1, 1] / sum(cm[, 1]),
  accuracy    = sum(diag(cm)) / sum(cm)
)
print(cm); print(perf)
write.xlsx(perf, file.path(out_dir, "validation_performance.xlsx"))

# -----------------------------------------------------------------------------
# ROC curve
# -----------------------------------------------------------------------------
roc_df <- data.frame(prediction = pred_n, labels = y, id = "hospitalized")
p <- ggplot(roc_df, aes(m = prediction, d = labels, colour = id)) +
  geom_roc(labels = FALSE, pointsize = 0, n.cuts = 0) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.1) +
  scale_colour_manual(values = c(hospitalized = "salmon4")) +
  annotate("text", x = 0.65, y = 0.20, size = 5,
           label = sprintf("AUC: %.2f 95%% CI (%.2f - %.2f)", perf$AUC, perf$AUC_low, perf$AUC_high)) +
  annotate("text", x = 0.54, y = 0.10, size = 5,
           label = paste("p-value:", format(p_auc, scientific = TRUE, digits = 3))) +
  labs(title = "Cohort CRA003945 - ambulatory vs hospitalized",
       x = "False Positive Rate (1-Specificity)", y = "True Positive Rate (Sensitivity)") +
  theme_article()
save_plot(p, file.path(out_dir, "ROC_severity_CRA003945.pdf"), 12, 12)
