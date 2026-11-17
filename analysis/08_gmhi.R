# =============================================================================
# 08 - Gut Microbiome Health Index (GMHI)
# -----------------------------------------------------------------------------
#   Rscript analysis/08_gmhi.R
#
#   Fig 7a          : GMHI by mortality (alive vs deceased) and by
#                     COVID-19 status x mortality (4 groups)
#   Fig 7b          : GMHI in healthy GCAT individuals vs non-COVID-19 vs
#                     COVID-19 cancer patients (baseline samples)
#   Supp Table 14   : age and sex, GCAT healthy vs cancer cohort
#
# Healthy GCAT individuals: no cancer, no adenoma/polyp and none of the listed
# comorbidities (n = 122).
# =============================================================================

source("config/config.R")
source("config/comparisons.R")
source("R/utils.R")
source("R/gmhi.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggpubr)
  library(rstatix)
  library(gtsummary)
  library(openxlsx)
})

out_dir <- out_dir_for("08_gmhi")

gmhi_boxplot <- function(df, group_var, colours, labels, title) {
  f  <- as.formula(paste("GMHI ~", group_var))
  pw <- df %>% wilcox_test(f, p.adjust.method = "BH") %>%
    add_xy_position(x = group_var, step.increase = 0.08)
  p <- ggboxplot(df, x = group_var, y = "GMHI", color = group_var, palette = unname(colours),
                 title = title) +
    stat_pvalue_manual(pw, label = "p.adj", tip.length = 0.01) +
    scale_x_discrete(labels = labels) +
    labs(x = NULL, y = "Gut Microbiome Health Index (GMHI)") +
    theme_article() + rotate_x_text(angle = 45)
  if (nlevels(df[[group_var]]) > 2) p <- p + stat_compare_means(method = "kruskal.test", label.y.npc = "bottom")
  list(plot = p, stats = pw)
}

# =============================================================================
# 1. MARATÓ cohort - mortality (Fig 7a)
# =============================================================================
ps   <- load_marato()
meta <- get_metadata(ps)
gmhi <- compute_gmhi(species_profile(ps))
meta$GMHI <- gmhi[rownames(meta), "GMHI"]
write.xlsx(meta[, c("COHORT", "tipus_kit", "basal_line", "outcome", "GMHI")],
           file.path(out_dir, "GMHI_marato.xlsx"), rowNames = TRUE)

stats <- list()
for (name in c("mortality", "mortality_4groups")) {
  cmp <- get_comparison(name)
  df  <- meta[cmp$keep(meta), ]
  df$group <- cmp$group(df)
  res <- gmhi_boxplot(df, "group", cmp$colours, cmp$labels, "GMHI - Mortality")
  save_plot(res$plot, file.path(out_dir, paste0("Fig7a_GMHI_", name, ".pdf")), 12, 13)
  stats[[name]] <- as.data.frame(res$stats[, setdiff(names(res$stats), "groups")])
}

# =============================================================================
# 2. GCAT healthy reference vs cancer cohort (Fig 7b)
# =============================================================================
ps_gcat   <- readRDS(PATHS$phyloseq_gcat)
meta_gcat <- get_metadata(ps_gcat)

comorbidities <- c("Obesity", "Diabetes", "HTA", "Hypercholesterolemia",
                   "Cardiobascular_disease", "Digestive_disease", "Breast_disease",
                   "Psoriasis", "Arthritis", "Autoimmune_diseases", "Respiratory_diseases",
                   "Rinitis", "EPOC", "Asma")
unhealthy <- meta_gcat$cancer %in% "yes" | meta_gcat$adenoma_polyp %in% "yes" |
             rowSums(sapply(meta_gcat[comorbidities], function(x) x %in% "Yes")) > 0
meta_gcat <- meta_gcat[!unhealthy, ]
message(nrow(meta_gcat), " healthy GCAT individuals")

# GCAT taxa names are full lineages separated by ";" -> keep "s__<species>"
sp_gcat <- species_profile(ps_gcat)[, rownames(meta_gcat)]
rownames(sp_gcat) <- paste0("s__", sub(".*;", "", rownames(sp_gcat)))

# MARATÓ baseline samples
meta_base <- meta[meta$tipus_kit %in% "id_kit_basal", ]
sp_mar    <- species_profile(ps)[, rownames(meta_base)]

# Union of species; absent species = 0
all_sp <- union(rownames(sp_mar), rownames(sp_gcat))
fill <- function(m) {
  out <- matrix(0, length(all_sp), ncol(m), dimnames = list(all_sp, colnames(m)))
  out[rownames(m), ] <- m
  out
}
sp_all  <- cbind(fill(sp_mar), fill(sp_gcat))
gmhi_all <- compute_gmhi(sp_all)

combined <- bind_rows(
  data.frame(sample = rownames(meta_base), id = as.character(meta_base$nhc_sap),
             age = as.numeric(meta_base$age), sex = as.character(meta_base$sex),
             cohort = "Cancer", group = as.character(meta_base$COHORT)),
  data.frame(sample = rownames(meta_gcat), id = as.character(meta_gcat$ID),
             age = as.numeric(meta_gcat$Age), sex = as.character(meta_gcat$Sex),
             cohort = "Healthy", group = "Healthy")
) %>%
  mutate(sex    = recode(sex, M = "Male", F = "Female"),
         cohort = factor(cohort, levels = c("Healthy", "Cancer")),
         group  = factor(group, levels = c("Healthy", "CONTROL", "CASE")),
         GMHI   = gmhi_all[sample, "GMHI"])
print(table(combined$group))
write.xlsx(combined, file.path(out_dir, "GMHI_GCAT_marato.xlsx"))

res <- gmhi_boxplot(combined, "group",
                    c("lightblue4", "darkslategrey", "salmon4"),
                    c("Healthy (GCAT)", "Non-COVID-19", "COVID-19"),
                    "GMHI - Healthy vs cancer")
save_plot(res$plot, file.path(out_dir, "Fig7b_GMHI_healthy_cancer.pdf"), 12, 13)
stats$gcat <- as.data.frame(res$stats[, setdiff(names(res$stats), "groups")])
write.xlsx(stats, file.path(out_dir, "GMHI_tests.xlsx"))

# =============================================================================
# 3. Supp Table 14 - age and sex, healthy vs cancer
# =============================================================================
tbl <- combined %>%
  select(cohort, age, sex) %>%
  tbl_summary(by = cohort, statistic = list(age ~ "{mean} ± {sd}"),
              label = list(age ~ "Age (years)", sex ~ "Sex")) %>%
  add_p(test = list(age ~ "t.test", sex ~ "chisq.test")) %>%
  modify_header(label ~ "**Characteristic**") %>%
  bold_labels()
write.xlsx(as_tibble(tbl), file.path(out_dir, "SuppTable14_age_sex_GCAT_vs_cancer.xlsx"))
