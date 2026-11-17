# Gut microbiota dysbiosis is associated with COVID-19 severity and mortality in cancer patients

Code for the microbiome analyses of:

> Santiago A\*, Moragas N\*, Ruiz-Pace F, *et al.* **Gut microbiota dysbiosis is associated with COVID-19 severity and mortality in cancer patients.** *(in preparation, 2026)*
> \* These authors contributed equally.

Prospective multicentre study (Vall d'Hebron University Hospital and Catalan Institute of Oncology, 2020–2023) of patients with solid tumours with and without COVID-19. Stool samples were analysed by shotgun metagenomic sequencing, together with faecal cytokines and faecal SARS-CoV-2 RNA. Funded by La Marató de TV3.

---

## Repository structure

```
├── config/
│   ├── config.R              # paths, filters, seeds
│   └── comparisons.R         # definition of every clinical comparison
├── R/
│   ├── utils.R               # data loading, filtering, zero replacement, CLR, plot helpers
│   └── gmhi.R                # Gut Microbiome Health Index
├── preprocessing/
│   ├── shotgun_pipeline.sh   # raw reads -> species abundance table
│   └── config.example.sh
├── analysis/
│   ├── 02_diversity.R             # alpha / beta diversity        (generic)
│   ├── 03_differential_abundance.R# ANCOM-BC + LinDA               (generic)
│   ├── 03b_da_overlap.R           # overlap of DA species across contrasts
│   ├── 04_lasso_model.R           # LASSO prediction model
│   ├── 05_external_validation.R   # validation in CRA003945
│   ├── 06_follow_up.R             # longitudinal active vs recovered
│   ├── 07_cytokines.R             # faecal cytokines
│   └── 08_gmhi.R                  # GMHI, GCAT healthy reference
├── data/README.md            # required input files
└── install_packages.R
```

The scripts `02_diversity.R` and `03_differential_abundance.R` are generic: the
comparison is passed as an argument and defined in `config/comparisons.R`
(samples included, grouping variable, reference level, colours, covariates).

| Comparison | Description |
|---|---|
| `case_control` | COVID-19 vs non-COVID-19 (baseline samples) |
| `timepoint` | non-COVID-19 vs recovered vs recovering vs active infection (all samples) |
| `severity` | ambulatory (asymptomatic/mild) vs hospitalised COVID-19 |
| `severity_3levels` | asymptomatic/mild vs moderate vs severe |
| `mortality` | alive vs deceased (whole cohort) |
| `mortality_4groups` | COVID-19 status × mortality |
| `mortality_covid`, `mortality_noncovid` | alive vs deceased within each stratum |

## Workflow

### 1. Pre-processing (HPC, SGE)

```bash
cd preprocessing
cp config.example.sh config.sh     # edit paths
qsub shotgun_pipeline.sh
```

| Step | Tool | Settings |
|---|---|---|
| Human read removal | Bowtie2 2.3.4 + samtools 1.8 | GRCh38, `--very-sensitive-local`, keep pairs with both mates unmapped |
| QC | FastQC + MultiQC 1.10.1 | before and after trimming |
| Deduplication | BBTools clumpify 38.26 | `dedupe subs=1 k=11 passes=3` |
| Trimming / contaminants | BBTools bbduk 38.26 | adapters, artefacts, PhiX; `qtrim=rl trimq=20 minlength=75` |
| Taxonomic profiling | Kraken2 2.1.2 (UHGG v2.0.2) | `--confidence 0.1` |
| Abundance re-estimation | Bracken 2.2 | species level, `-r 150 -t 10` |
| Abundance table | KrakenTools | `kreport2mpa.py`, `combine_mpa.py` |

The external validation cohort (CRA003945, NGDC Genome Sequence Archive) was
downloaded and processed with the same pipeline.

### 2. Statistical analysis (R 4.3.3)

```bash
Rscript install_packages.R

for c in case_control timepoint severity mortality mortality_4groups; do
  Rscript analysis/02_diversity.R $c
done
for c in case_control timepoint severity mortality mortality_covid mortality_noncovid; do
  Rscript analysis/03_differential_abundance.R $c
done
Rscript analysis/03b_da_overlap.R
Rscript analysis/04_lasso_model.R severity
Rscript analysis/05_external_validation.R
Rscript analysis/06_follow_up.R
Rscript analysis/07_cytokines.R
Rscript analysis/08_gmhi.R
```

All scripts are run from the repository root; results are written to `results/`.

Common processing of the abundance tables:
species with relative abundance > 0.1% in ≥ 5% of samples are retained;
zeros are replaced with a Bayesian-multiplicative method (`zCompositions::cmultRepl`, SQ);
data are CLR-transformed. Alpha diversity is computed on rarefied (Shannon) and
unrarefied (Chao1) counts.

## Figures and tables

| Figure / table | Script |
|---|---|
| Fig 2a, Supp Fig 1a | `02_diversity.R case_control` |
| Fig 2b, Supp Table 1 | `03_differential_abundance.R case_control` |
| Fig 2c, Supp Fig 1b, Supp Table 2 | `02_diversity.R timepoint` |
| Fig 2d, Supp Table 3 | `03_differential_abundance.R timepoint` |
| Fig 2e, Supp Fig 1e | `02_diversity.R severity` |
| Fig 2f, Supp Table 4 | `03_differential_abundance.R severity` |
| Fig 2g | `03b_da_overlap.R` |
| Supp Fig 1d | `06_follow_up.R` |
| Fig 4d, Fig 4e | `07_cytokines.R` |
| Fig 6a, Supp Table 12 | `02_diversity.R mortality` |
| Fig 6b, Supp Table 11 | `03_differential_abundance.R mortality` |
| Supp Fig 3a | `02_diversity.R mortality_4groups` |
| Supp Fig 3c | `07_cytokines.R` |
| Fig 7a, Fig 7b, Supp Table 14 | `08_gmhi.R` |
| Fig 8a | `04_lasso_model.R severity` |
| Fig 8b | `05_external_validation.R` |

### Analyses not included in this repository

<!-- TODO: confirm the author of each analysis -->
The following analyses were performed by **Fiorella Ruiz-Pace** (Paolo Nuciforo's
group, Vall d'Hebron Institute of Oncology) and are not part of this repository:

- Table 1 and Table 2 (clinical characteristics, systemic inflammatory markers)
- Systemic inflammatory markers and their correlation with species (Fig 3, Fig 6c, Supp Fig 3b)
- Faecal SARS-CoV-2 RNA (Fig 4a–c)
- Multivariable and Firth logistic regression models of mortality (Fig 5, Supp Fig 2, Supp Tables 8–10)
- Faecal cytokines in the CRC subgroup (Supp Fig 3d)
- External colorectal cancer cohorts (Fig 7c–d)
- Mortality prediction models (Fig 8c–e)
- Trend analysis across infection stages (Supp Fig 1c)

## Data availability

Raw sequencing data: *[accession to be added]*. Clinical data are available
from the corresponding authors upon reasonable request. See `data/README.md`
for the files expected by the scripts.

## Acknowledgements

Analysis functions adapted from the Unit of Biomarkers and Susceptibility
(Catalan Institute of Oncology): prevalence filter (`filter_matrix`, J. Mas).
GMHI implementation adapted from Gupta *et al.*, *Nat Commun* 2020.

## Contact

Núria Moragas — *[email]*

## License

MIT (see `LICENSE`).
