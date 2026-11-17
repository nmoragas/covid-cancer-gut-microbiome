# =============================================================================
# Gut Microbiome Health Index (GMHI)
# -----------------------------------------------------------------------------
# Adapted from Gupta VK et al., Nat Commun 2020;11:4635
# (https://github.com/jaeyunsung/GMHI_2020). Original author: Vinod K. Gupta.
#
# Input: species x samples abundance table with "s__" species names.
# The GMHI health-prevalent (MH, n = 7) and health-scarce (MN, n = 43) species
# are named after the MetaPhlAn2/NCBI taxonomy; UHGG/GTDB synonyms are mapped
# below. Species without a UHGG equivalent are not included.
# =============================================================================

# UHGG (GTDB) name -> GMHI (NCBI) name
GMHI_SYNONYMS <- c(
  "s__Lancefieldella_parvula"          = "s__Atopobium_parvulum",
  "s__Enterocloster_asparagiformis"    = "s__Clostridium_asparagiforme",
  "s__Enterocloster_bolteae"           = "s__Clostridium_bolteae",
  "s__Enterocloster_citroniae"         = "s__Clostridium_citroniae",
  "s__Enterocloster_clostridioformis"  = "s__Clostridium_clostridioforme",
  "s__Hungatella_A_hathewayi_A"        = "s__Clostridium_hathewayi",
  "s__Erysipelatoclostridium_ramosum"  = "s__Clostridium_ramosum",
  "s__Ligilactobacillus_salivarius"    = "s__Lactobacillus_salivarius",
  "s__Ruminococcus_B_gnavus"           = "s__Ruminococcus_gnavus",
  "s__Bulleidia_moorei"                = "s__Solobacterium_moorei",
  "s__Streptococcus_infantis_I"        = "s__Streptococcus_infantis",
  "s__Streptococcus_sanguinis_H"       = "s__Streptococcus_sanguinis"
)

compute_gmhi <- function(species_profile,
                         MH_species = scan(PATHS$gmhi_MH, what = "", sep = "\n", quiet = TRUE),
                         MN_species = scan(PATHS$gmhi_MN, what = "", sep = "\n", quiet = TRUE),
                         MH_prime = 7, MN_prime = 31) {
  sp <- as.matrix(species_profile)

  # remove unclassified and viral entries
  sp <- sp[!grepl("unclassified|virus", rownames(sp)), , drop = FALSE]

  # map UHGG names to the GMHI (NCBI) names
  hit <- rownames(sp) %in% names(GMHI_SYNONYMS)
  rownames(sp)[hit] <- GMHI_SYNONYMS[rownames(sp)[hit]]

  # relative abundances; values < 1e-5 set to 0
  sp <- sweep(sp, 2, colSums(sp), "/")
  sp[sp < 1e-5] <- 0

  message(sprintf("GMHI species found: MH %d/%d, MN %d/%d",
                  sum(MH_species %in% rownames(sp)), length(MH_species),
                  sum(MN_species %in% rownames(sp)), length(MN_species)))

  MH <- sp[rownames(sp) %in% MH_species, , drop = FALSE]
  MN <- sp[rownames(sp) %in% MN_species, , drop = FALSE]

  shannon <- function(x) -sum(log(x[x > 0]) * x[x > 0])
  psi_MH <- (colSums(MH > 0) / MH_prime) * apply(MH, 2, shannon)
  psi_MN <- (colSums(MN > 0) / MN_prime) * apply(MN, 2, shannon)

  data.frame(sample = colnames(sp),
             GMHI = log10((psi_MH + 1e-5) / (psi_MN + 1e-5)),
             row.names = colnames(sp))
}

#' Species x samples table with "s__" names (no further cleaning)
species_profile <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) m <- t(m)
  m
}
