# =============================================================================
# 06 - Longitudinal follow-up: active infection vs recovery
# -----------------------------------------------------------------------------
#   Rscript analysis/06_follow_up.R                 # Supp Fig 1d
#
# Patients with a baseline sample during active infection ("Activo") and a
# follow-up sample after recovery ("Resuelto"). Patients who died during
# follow-up are excluded (final n = 15 patients, 30 samples).
#
# Taxa analysed (candidate permissive / protective taxa, see 03b_da_overlap.R):
#   - 4 species enriched in all three contrasts (COVID-19, active, severity)
#   - 8 species enriched in both the COVID-19 and active-infection contrasts
#   - Romboutsia timonensis (only species depleted in active infection)
#
# Test: paired Wilcoxon signed-rank test (patient as pairing unit), BH-adjusted.
# =============================================================================

source("config/config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(openxlsx)
})

out_dir <- out_dir_for("06_follow_up")

TAXA <- c(
  # enriched in COVID-19, active infection and severity
  "Enterococcus_A_avium", "Enterococcus_B_faecium", "Enterococcus_faecalis",
  "Hungatella_sp005845265",
  # enriched in COVID-19 and active infection
  "Alistipes_A_sp900240235", "Bacteroides_sp900765785", "CAG.279_sp000437795",
  "Enterocloster_bolteae", "Paraprevotella_clara", "Paraprevotella_xylaniphila",
  "Phascolarctobacterium_faecium", "Prevotella_sp900765465",
  # depleted in active infection
  "Romboutsia_timonensis"
)

# -----------------------------------------------------------------------------
# Patient selection
# -----------------------------------------------------------------------------
ps   <- load_marato()
meta <- get_metadata(ps)
meta$sample <- rownames(meta)
meta$Fecha_muestra <- as.Date(meta$Fecha_muestra)

paired <- meta %>%
  filter((tipus_kit == "id_kit_basal"     & Status_.PN. == "Activo") |
         (tipus_kit == "id_kit_follow_up" & Status_.PN. == "Resuelto")) %>%
  arrange(nhc_sap, tipus_kit, Fecha_muestra) %>%
  distinct(nhc_sap, tipus_kit, .keep_all = TRUE) %>%     # first sample per kit
  group_by(nhc_sap) %>% filter(n() == 2) %>% ungroup() %>%
  filter(outcome != "Mort" | is.na(outcome))             # exclude deaths

intervals <- paired %>%
  select(nhc_sap, tipus_kit, Fecha_muestra) %>%
  pivot_wider(names_from = tipus_kit, values_from = Fecha_muestra) %>%
  mutate(days = as.numeric(id_kit_follow_up - id_kit_basal))
message(sprintf("Patients: %d | interval median %.1f days [IQR %.1f-%.1f]",
                nrow(intervals), median(intervals$days),
                quantile(intervals$days, 0.25), quantile(intervals$days, 0.75)))
write.xlsx(intervals, file.path(out_dir, "follow_up_patients.xlsx"))

# -----------------------------------------------------------------------------
# CLR abundances (zero replacement on all species of the selected samples)
# -----------------------------------------------------------------------------
ps_sub <- prune_samples(paired$sample, ps)
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
clr_abund <- clr_transform(replace_zeros(abundance_table(ps_sub)))

missing <- setdiff(TAXA, colnames(clr_abund))
if (length(missing)) warning("Taxa not found: ", paste(missing, collapse = ", "))

dat_long <- as.data.frame(clr_abund[, intersect(TAXA, colnames(clr_abund))]) %>%
  mutate(sample = rownames(.)) %>%
  left_join(paired %>% select(sample, nhc_sap, Status_.PN.), by = "sample") %>%
  pivot_longer(-c(sample, nhc_sap, Status_.PN.), names_to = "taxon", values_to = "CLR") %>%
  mutate(status = factor(Status_.PN., levels = c("Resuelto", "Activo"),
                         labels = c("Recovered", "Active")),
         taxon  = factor(taxon, levels = TAXA))

# -----------------------------------------------------------------------------
# Paired Wilcoxon test per taxon
# -----------------------------------------------------------------------------
stats <- dat_long %>%
  select(nhc_sap, taxon, status, CLR) %>%
  pivot_wider(names_from = status, values_from = CLR) %>%
  group_by(taxon) %>%
  summarise(n_pairs = sum(complete.cases(Recovered, Active)),
            median_active    = median(Active, na.rm = TRUE),
            median_recovered = median(Recovered, na.rm = TRUE),
            p = wilcox.test(Active, Recovered, paired = TRUE, exact = FALSE)$p.value,
            .groups = "drop") %>%
  mutate(p_BH = p.adjust(p, method = "BH")) %>%
  arrange(p)
print(stats)
write.xlsx(stats, file.path(out_dir, "follow_up_paired_wilcoxon.xlsx"))

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
labels_p <- stats %>% mutate(label = paste0("p = ", signif(p, 2)))
p <- ggplot(dat_long, aes(status, CLR, colour = status)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_line(aes(group = nhc_sap), colour = "grey80", linewidth = 0.3) +
  geom_point(size = 1, colour = "grey30") +
  geom_text(data = labels_p, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 3.5) +
  facet_wrap(~ taxon, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = c(Recovered = "lightblue4", Active = "salmon4")) +
  labs(x = NULL, y = "CLR-transformed abundance") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(face = "italic"))
save_plot(p, file.path(out_dir, "follow_up_boxplots.pdf"), 30, 30)
