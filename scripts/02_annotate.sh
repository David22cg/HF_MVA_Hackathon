#!/usr/bin/env bash
#   It attaches to each variant (in INFO/CSQ) evidence: 
#     - gene (SYMBOL) and consequence  -> which gene and how severe
#     - gnomAD (AF)                    -> rarity: a rare-disease cause is rare
#     - ClinVar (CLNSIG)               -> already-curated clinical classification
#     - AlphaMissense (am_patho...)    -> missense pathogenicity (VUS)
#   all with local databases
#
#   Requires: VEP GRCh38 cache + FASTA + ClinVar + AlphaMissense
set -euo pipefail

#paths
WORK="$HOME/mva/work"
CACHE="$HOME/mva/cache"
REF="$HOME/mva/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DATA="$HOME/mva/data"
PLUGINS="$HOME/.vep/Plugins"
THREADS=4

IN="$WORK/all.norm.vcf.gz"
OUT="$WORK/annotated.vcf.gz"
[[ -f "$IN" ]] || { echo "ERROR: $IN does not exist (run 01_prep.sh first)"; exit 1; }

VEP_ARGS=(
  --offline --cache --dir_cache "$CACHE" --fasta "$REF"   # local
  --species homo_sapiens --assembly GRCh38                # challenge build
  --input_file "$IN" --format vcf
  --vcf --compress_output bgzip                           # annotated+bgzip VCF out (keeps GT/AD/DP)
  --output_file "$OUT" --force_overwrite --stats_text
  --fork "$THREADS" --buffer_size 20000                   # speed on >500k variants
  --symbol --canonical                                    # gene SYMBOL + canonical transcript
  --af_gnomadg --af_gnomade                               # gnomAD frequencies (genome and exome)
  --sift b --polyphen b                                   # classic predictors (support)
  --pick_allele_gene                                      # 1 consequence per allele and gene (manageable CSQ)
  --dir_plugins "$PLUGINS"
  --custom "$DATA/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN"
  --plugin AlphaMissense,file="$DATA/AlphaMissense_hg38.ens.tsv.gz"
)

# SpliceAI: predicts creation/loss of splice sites (the typical mechanism of the BUB1B 2nd allele). 
#It needs 2 precomputed VCFs (~30 GB total) downloaded from Illumina/BaseSpace with registration.
USE_SPLICEAI=1
SPLICEAI_SNV="$DATA/spliceai_scores.raw.snv.hg38.vcf.gz"
SPLICEAI_INDEL="$DATA/spliceai_scores.raw.indel.hg38.vcf.gz"
#
# LOFTEE: confirms a LoF is really loss-of-function. Lighter than SpliceAI but needs the LOFTEE repo + human_ancestor.fa.gz (~1 GB).
USE_LOFTEE=1
LOFTEE_PATH="$HOME/.vep/Plugins/loftee"
HUMAN_ANCESTOR="$DATA/human_ancestor.fa.gz"

if [[ "${USE_SPLICEAI:-0}" == "1" ]]; then
  [[ -f "$SPLICEAI_SNV" && -f "$SPLICEAI_INDEL" ]] || { echo "ERROR: SpliceAI VCFs missing"; exit 1; }
  VEP_ARGS+=(--plugin "SpliceAI,snv=$SPLICEAI_SNV,indel=$SPLICEAI_INDEL")
  echo "[vep] SpliceAI ENABLED (03's rescue track will use spliceai_min)"
fi
if [[ "${USE_LOFTEE:-0}" == "1" ]]; then
  VEP_ARGS+=(--plugin "LoF,loftee_path:$LOFTEE_PATH,human_ancestor_fa:$HUMAN_ANCESTOR")
  echo "[vep] LOFTEE ENABLED"
fi


vep "${VEP_ARGS[@]}"
tabix -p vcf "$OUT"
echo "[vep] done -> $OUT"

