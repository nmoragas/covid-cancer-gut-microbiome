# =============================================================================
# Configuration for shotgun_pipeline.sh
# Copy to config.sh and edit:  cp config.example.sh config.sh
# =============================================================================

# Folder with raw paired-end FASTQ files, one sub-folder per sample:
#   ${RAW}/<sample>/<sample>_1.fq.gz  and  ${RAW}/<sample>/<sample>_2.fq.gz
RAW="/path/to/raw_fastq"

# Bowtie2 index of the human reference genome (GRCh38)
HUMAN_INDEX="/path/to/references/human/hg38/bowtie2/genome"

# Kraken2 database built from UHGG v2.0.2 (must include Bracken kmer distributions)
KRAKEN_DB="/path/to/uhgg/2.0.2/kraken2_uhgg_v2_0_2"

# Path to KrakenTools scripts (https://github.com/jenniferlu717/KrakenTools)
KRAKENTOOLS="/path/to/KrakenTools"

# Quality-trimming threshold (bbduk trimq)
QUALITY=20

# Bracken parameters
READ_LENGTH=150
BRACKEN_THRESHOLD=10

# Working / output folders
WORKDIR="./temp"
OUTDIR="./out"
