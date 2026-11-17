# =============================================================================
# Post-Kraken processing: abundance table -> phyloseq object (faecal samples)
# -----------------------------------------------------------------------------
#   Rscript preprocessing/post_kraken.R
#
#   1. Sample clean-up
#      - remove low-quality sequencing replicates (re-sequenced by Novogene)
#      - fix sample naming errors
#      - simplify column names to the sample ID (e.g. "F_1372")
#   2. Normalisation by genome length: Bracken reads / genome length (Mb)
#      of the UHGG v2.0.2 species representative genome
#   3. phyloseq object: OTU table (species x samples), taxonomy
#      (Phylum ... Species) and sample metadata
#
# Input : out/species/bracken_abundance_species_mpa.txt (shotgun_pipeline.sh)
#         UHGG genome metadata with species names matching Kraken output
#         clinical metadata (one row per kit, id_kit without the "F_" prefix)
# Output: phyloseq_obj_F.rds
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(readxl)
  library(dplyr)
  library(readr)
})

# --- Configuration -----------------------------------------------------------
IN_MPA        <- "out/species/bracken_abundance_species_mpa.txt"
GENOME_LENGTH <- "genomes_all_metadata_equiparacio_kraken.tsv"   # UHGG v2.0.2 metadata
METADATA      <- "Metadata_COVID_ICO_v5.xlsx"
OUT_RDS       <- "phyloseq_obj_F.rds"
OUT_DIR       <- "post_kraken"
dir.create(OUT_DIR, showWarnings = FALSE)

# Low-quality replicates (lane L1) re-sequenced in a later batch
REMOVE_REGEX <- c("^F_1726_.*_L1", "^F_1619_.*_L1", "^F_1076_.*_L1", "^F_1235_.*_L1",
                  "^F_1375_.*_L1", "^F_1633_.*_L1", "^F_225_.*_L1",  "^F_906_.*_L1",
                  "^F_910_.*_L1",  "^F_1308_.*_L1")
# Duplicated runs: the run with fewer reads is removed
REMOVE_COLUMNS <- c(
  "F_1372_EKDN230003858.1A_HNWJTDSX5_L4_nohuman_Q20.report_bracken.txt",
  "F_190_EKDN230006935.1A_HVFH2DSX5_L1_nohuman_Q20.report_bracken.txt",
  "F_190_EKDN230006935.1A_HW7GHDSX5_L3_nohuman_Q20.report_bracken.txt",
  "F_594_EKDN230006896.1A_HVGLKDSX5_L2_nohuman_Q20.report_bracken.txt"
)
# Naming errors from the sequencing provider
RENAME_COLUMNS <- c(
  "F_652_EKDN230020851.1A_HFHJNDSX7_L3_nohuman_Q20.report_bracken.txt" =
    "F_651_EKDN230020851.1A_HFHJNDSX7_L3_nohuman_Q20.report_bracken.txt",
  "CTRL_NEG_F2_ZKDN230031154.1A_HVVLCDSX7_L1_nohuman_Q20.report_bracken.txt" =
    "CTRL_F2_ZKDN230031154.1A_HVVLCDSX7_L1_nohuman_Q20.report_bracken.txt"
)

# =============================================================================
# 1. Sample clean-up
# =============================================================================
abund <- read.delim(IN_MPA, check.names = FALSE, comment.char = "")
rownames(abund) <- abund[[1]]          # first column: "#Classification" (lineage)
abund <- abund[, -1]

drop <- grepl(paste(REMOVE_REGEX, collapse = "|"), colnames(abund)) |
        colnames(abund) %in% REMOVE_COLUMNS
abund <- abund[, !drop]
hit <- colnames(abund) %in% names(RENAME_COLUMNS)
colnames(abund)[hit] <- RENAME_COLUMNS[colnames(abund)[hit]]

# Faecal samples and negative controls only (F_*, CTRL_F*)
abund <- abund[, grepl("^F_|^CTRL_F", colnames(abund))]
colnames(abund) <- sub("^([A-Z0-9]+_[A-Z0-9]+).*", "\\1", colnames(abund))
stopifnot(!any(duplicated(colnames(abund))))
abund <- abund[rowSums(abund) > 0, ]

write.table(data.frame(sample = colnames(abund), bracken_reads = colSums(abund)),
            file.path(OUT_DIR, "bracken_reads_per_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
# 2. Normalisation by genome length
# =============================================================================
species <- sub(".*s__([^|]+).*", "s__\\1", rownames(abund))
lineage <- setNames(rownames(abund), species)
rownames(abund) <- species

genome_length <- read_tsv(GENOME_LENGTH, show_col_types = FALSE) %>%
  distinct(Lineage, .keep_all = TRUE) %>%
  mutate(length_mb = as.double(Length) / 1e6)
len <- setNames(genome_length$length_mb, genome_length$species)

missing <- setdiff(rownames(abund), names(len))
if (length(missing)) stop(length(missing), " species without genome length, e.g. ", missing[1])

abund_norm <- sweep(as.matrix(abund), 1, len[rownames(abund)], "/")

# =============================================================================
# 3. phyloseq object
# =============================================================================
# Taxonomy from the lineage string "p__...|c__...|...|s__..."
ranks <- c(Phylum = "p", Class = "c", Order = "o", Family = "f", Genus = "g", Species = "s")
taxmat <- sapply(ranks, function(r) {
  x <- lineage[rownames(abund_norm)]
  out <- sub(paste0(".*", r, "__([^|;]+).*"), paste0(r, "__\\1"), x)
  out[out == x] <- NA
  out
})
rownames(taxmat) <- rownames(abund_norm)

metadata <- as.data.frame(read_excel(METADATA))
metadata$Fecha_dx_COVID <- as.character(metadata$Fecha_dx_COVID, format = "%Y-%m-%d")
metadata[metadata == "<NA>"] <- NA
metadata <- metadata[metadata$Feces_sequenced %in% "Yes", ]
metadata$id_kit <- ifelse(startsWith(metadata$id_kit, "CTRL_"), metadata$id_kit,
                          paste0("F_", metadata$id_kit))
rownames(metadata) <- metadata$id_kit

not_in_meta <- setdiff(colnames(abund_norm), rownames(metadata))
if (length(not_in_meta)) message("Samples without metadata (dropped): ",
                                 paste(not_in_meta, collapse = ", "))
samples <- intersect(colnames(abund_norm), rownames(metadata))

ps <- phyloseq(otu_table(abund_norm[, samples], taxa_are_rows = TRUE),
               tax_table(taxmat),
               sample_data(metadata[samples, ]))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
print(ps)
saveRDS(ps, OUT_RDS)
