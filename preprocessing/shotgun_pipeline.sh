#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -pe smp 8
#$ -l mf=28G
#$ -N shotgun_pipeline
#$ -e logs
#$ -o logs
# =============================================================================
# Shotgun metagenomics pre-processing pipeline (species level)
# -----------------------------------------------------------------------------
#   1. Human read removal ............ Bowtie2 (--very-sensitive-local) + samtools
#   2. Quality control (raw) ......... FastQC + MultiQC
#   3. Deduplication + trimming ...... BBTools clumpify + bbduk
#   4. Quality control (clean) ....... FastQC + MultiQC
#   5. Taxonomic profiling ........... Kraken2 (UHGG v2.0.2) + Bracken (species)
#   6. Abundance table ............... KrakenTools kreport2mpa + combine_mpa
#
# Usage (SGE):   qsub shotgun_pipeline.sh            (reads ./config.sh)
#       (local): bash shotgun_pipeline.sh config.sh
#
# Output: ${OUTDIR}/species/bracken_abundance_species_mpa.txt
#         ${OUTDIR}/read_counts_clean.tsv
# Next step: Rscript post_kraken.R (genome-length normalisation + phyloseq)
# The same pipeline was used for the MARATÓ cohort and for the external
# validation cohort (CRA003945).
# =============================================================================

set -euo pipefail

CONFIG="${1:-config.sh}"
# shellcheck source=/dev/null
source "${CONFIG}"

THREADS="${NSLOTS:-8}"
MEM="28g"

module load apps/java-9.0.4
module load apps/bowtie2-2.3.4
module load apps/samtools-1.8
module load apps/bbmap-38.26
module load apps/multiqc-1.10.1
module load apps/kraken2-2.1.2
module load apps/bracken-2.2

D1="${WORKDIR}/1_human_remove/nohuman"
D2="${WORKDIR}/2_QC_before_trim"
D3="${WORKDIR}/3_dedup_trim"
D3_TMP="${D3}/intermediate"
D3_OUT="${D3}/clean"
D4="${WORKDIR}/4_QC_after_trim"
D5_K2="${WORKDIR}/5_kraken/kraken2"
D5_BR="${WORKDIR}/5_kraken/bracken"
D5_MPA="${WORKDIR}/5_kraken/mpa"

mkdir -p logs "${D1}" "${D2}" "${D3_TMP}" "${D3_OUT}" "${D4}" \
         "${D5_K2}/outputs" "${D5_K2}/reports" "${D5_BR}" "${D5_MPA}" "${OUTDIR}/species"

# Check that the expected number of files exists, otherwise stop
check_files () {   # $1 = folder, $2 = glob pattern, $3 = expected number
  local n
  n=$(find "$1" -maxdepth 1 -type f -name "$2" | wc -l)
  if [[ "$n" -ne "$3" ]]; then
    echo "ERROR: expected $3 files '$2' in $1, found $n" >&2
    exit 1
  fi
  echo "Checkpoint OK: $n files '$2' in $1"
}

# -----------------------------------------------------------------------------
# Sample list (one line per sample: <id> <R1> <R2>)
# -----------------------------------------------------------------------------
SAMPLES="${WORKDIR}/samples.txt"
: > "${SAMPLES}"
for r1 in "${RAW}"/*/*_1.*; do
  r2="${r1/_1./_2.}"
  id=$(basename "${r1}"); id="${id%%_1.*}"
  [[ -f "${r2}" ]] || { echo "Missing mate for ${r1}" >&2; exit 1; }
  echo "${id} $(realpath "${r1}") $(realpath "${r2}")" >> "${SAMPLES}"
done
N=$(wc -l < "${SAMPLES}")
echo "Samples found: ${N}"

# -----------------------------------------------------------------------------
# 1. Human read removal: keep pairs where both mates are unmapped (-f 12)
# -----------------------------------------------------------------------------
while read -r id r1 r2; do
  [[ -s "${D1}/${id}_nohuman.m1.fq" ]] && continue
  bowtie2 -x "${HUMAN_INDEX}" -p "${THREADS}" -1 "${r1}" -2 "${r2}" \
          --very-sensitive-local -k 1 \
    | samtools fastq -f 12 \
        -1 "${D1}/${id}_nohuman.m1.fq" \
        -2 "${D1}/${id}_nohuman.m2.fq" -
done < "${SAMPLES}"
check_files "${D1}" "*_nohuman.m?.fq" $((N * 2))

# -----------------------------------------------------------------------------
# 2. QC before trimming
# -----------------------------------------------------------------------------
fastqc -t "${THREADS}" -o "${D2}" "${D1}"/*.fq
multiqc "${D2}" --outdir "${D2}"

# -----------------------------------------------------------------------------
# 3. Deduplication (clumpify) + adapter/artefact/PhiX removal and
#    quality trimming (bbduk)
# -----------------------------------------------------------------------------
while read -r id _ _; do
  [[ -s "${D3_OUT}/${id}_Q${QUALITY}.m1.fq" ]] && continue

  clumpify.sh -Xmx${MEM} \
    in1="${D1}/${id}_nohuman.m1.fq" in2="${D1}/${id}_nohuman.m2.fq" \
    out1="${D3_TMP}/${id}_deduped.m1.fq" out2="${D3_TMP}/${id}_deduped.m2.fq" \
    dedupe subs=1 k=11 passes=3 \
    &> "${D3_TMP}/${id}_deduplicationstats.txt"

  bbduk.sh -Xmx${MEM} \
    in1="${D3_TMP}/${id}_deduped.m1.fq" in2="${D3_TMP}/${id}_deduped.m2.fq" \
    out1="${D3_OUT}/${id}_Q${QUALITY}.m1.fq" out2="${D3_OUT}/${id}_Q${QUALITY}.m2.fq" \
    ref=adapters,artifacts,phix \
    k=25 mink=6 hdist=1 hdist2=0 ktrim=r \
    qtrim=rl trimq="${QUALITY}" minlength=75 \
    stats="${D3_TMP}/${id}_adapters.txt" dump="${D3_TMP}/${id}_kmers.txt" \
    &> "${D3_TMP}/${id}_cleaning.txt"
done < "${SAMPLES}"
check_files "${D3_OUT}" "*_Q${QUALITY}.m?.fq" $((N * 2))

# Number of clean reads per file (input to Kraken2)
for fq in "${D3_OUT}"/*.fq; do
  echo -e "$(basename "${fq}")\t$(( $(wc -l < "${fq}") / 4 ))"
done > "${OUTDIR}/read_counts_clean.tsv"

# -----------------------------------------------------------------------------
# 4. QC after trimming
# -----------------------------------------------------------------------------
fastqc -t "${THREADS}" -o "${D4}" "${D3_OUT}"/*.fq
multiqc "${D4}" --outdir "${D4}"

# -----------------------------------------------------------------------------
# 5. Kraken2 classification + Bracken species re-estimation
# -----------------------------------------------------------------------------
while read -r id _ _; do
  base="${id}_Q${QUALITY}"
  kraken2 --db "${KRAKEN_DB}" \
    --confidence 0.1 --threads "${THREADS}" --use-names \
    --output "${D5_K2}/outputs/${base}.output.txt" \
    --report "${D5_K2}/reports/${base}.report.txt" \
    --paired "${D3_OUT}/${base}.m1.fq" "${D3_OUT}/${base}.m2.fq"

  bracken -d "${KRAKEN_DB}" -i "${D5_K2}/reports/${base}.report.txt" \
    -r "${READ_LENGTH}" -t "${BRACKEN_THRESHOLD}" -l S \
    -o "${D5_BR}/${base}_report_species.txt"
done < "${SAMPLES}"

# Bracken writes the re-estimated kraken-style report next to the input report
mv "${D5_K2}"/reports/*_bracken*.txt "${D5_BR}/"

# -----------------------------------------------------------------------------
# 6. Convert to MetaPhlAn-style (mpa) format and merge all samples
# -----------------------------------------------------------------------------
for rep in "${D5_BR}"/*.report_bracken*.txt; do
  base=$(basename "${rep}"); base="${base%%.report_bracken*}"
  python "${KRAKENTOOLS}/kreport2mpa.py" -r "${rep}" \
    -o "${D5_MPA}/${base}_mpa.txt" --display-header
done

python "${KRAKENTOOLS}/combine_mpa.py" -i "${D5_MPA}"/*_mpa.txt \
  -o "${OUTDIR}/species/combined_species_mpa.txt"
grep -E "(s__)|(#Classification)" "${OUTDIR}/species/combined_species_mpa.txt" \
  > "${OUTDIR}/species/bracken_abundance_species_mpa.txt"

echo "Done: ${OUTDIR}/species/bracken_abundance_species_mpa.txt"
