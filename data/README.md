# Input data

Patient-level data are not included in this repository. The analysis scripts
expect the following files in this folder (paths can be changed in
`config/config.R`):

| File | Content | Produced by / source |
|---|---|---|
| `genomes_all_metadata_equiparacio_kraken.tsv` | UHGG v2.0.2 genome metadata (`Lineage`, `Length`) with a `species` column matching the Kraken2 species names | UHGG v2.0.2 `genomes-all_metadata.tsv` (used by `preprocessing/post_kraken.R`) |
| `Metadata_COVID_ICO_v5.xlsx` | Clinical metadata, one row per kit | REDCap export (used by `preprocessing/post_kraken.R`) |
| `phyloseq_obj_F_uhgg_sel.rds` | phyloseq object: species abundances (Bracken, UHGG v2.0.2, normalised by genome length), UHGG taxonomy and sample metadata | `preprocessing/shotgun_pipeline.sh` + `preprocessing/post_kraken.R` |
| `phyloseq_obj_CRA003945.rds` | Same for the external validation cohort | Raw reads from the NGDC Genome Sequence Archive, accession CRA003945 |
| `GCAT_phyloseqobject_BEC.rds` | GCAT Microbiome cohort (healthy reference) | GCAT \| Genomes for Life (access upon request) |
| `cytokine_viral_load_covid.xlsx` | Faecal cytokines (sheet 1) and faecal SARS-CoV-2 RNA (sheet 2) | Laboratory measurements |
| `Cytokines data_pilot study.xlsx` | Faecal cytokines, pilot batch | Laboratory measurements |
| `gmhi/MH_species.txt`, `gmhi/MN_species.txt` | GMHI health-prevalent / health-scarce species | Gupta et al. 2020, https://github.com/jaeyunsung/GMHI_2020 |

## Main metadata variables

| Variable | Values |
|---|---|
| `COHORT` | `CASE` (COVID-19) / `CONTROL` (non-COVID-19) |
| `tipus_kit` | `id_kit_basal` / `id_kit_follow_up` / `id_kit_adicional` |
| `basal_line` | `Yes` for baseline samples |
| `Status_.PN.` | `No covid` / `Activo` (0–7 d) / `En resolución` (8–89 d) / `Resuelto` (≥90 d) |
| `covid_severity_modif` | `asymptomatic/mild disease` / `moderate_disease` / `severe_disease` |
| `outcome` | `Viu` (alive) / `Mort` (deceased) |
| `nhc_sap` | Patient identifier (pseudonymised) |
