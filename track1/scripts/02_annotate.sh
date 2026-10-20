#!/usr/bin/env bash
#   It attaches to each variant (in INFO/CSQ) evidence: 
#     - gene (SYMBOL) and consequence  -> which gene and how severe
#     - gnomAD (AF)                    -> rarity: a rare-disease cause is rare
#     - ClinVar (CLNSIG)               -> already-curated clinical classification
#     - AlphaMissense (am_patho...)    -> missense pathogenicity (VUS)
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
  --symbol --canonical                                    # gene symbol + canonical transcript
  --af_gnomadg --af_gnomade                               # gnomAD frequencies (genome and exome)
  --sift b --polyphen b                                   # classic predictors (support)
  --pick_allele_gene                                      # 1 consequence per allele and gene (manageable CSQ)
  --dir_plugins "$PLUGINS"
  --custom "$DATA/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN"
  --plugin AlphaMissense,file="$DATA/AlphaMissense_hg38.ens.tsv.gz"
)



#predicts creation/loss of splice sites
USE_SPLICEAI=1
SPLICEAI_SNV="$DATA/spliceai.hg38.vcf.gz"
SPLICEAI_INDEL="$DATA/spliceai.hg38.vcf.gz"
#
# LOFTEE: confirms a LoF is really loss-of-function 
 USE_LOFTEE=1
 LOFTEE_PATH="$HOME/.vep/Plugins/loftee"
 HUMAN_ANCESTOR="$DATA/human_ancestor.fa.gz"
 GERP_BW="$DATA/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
 LOFTEE_SQL="$DATA/loftee.sql"
 
if [[ "${USE_SPLICEAI:-0}" == "1" ]]; then
  [[ -f "$SPLICEAI_SNV" && -f "$SPLICEAI_INDEL" ]] || { echo "error: SpliceAI VCFs missing"; exit 1; }
  VEP_ARGS+=(--plugin "SpliceAI,snv=$SPLICEAI_SNV,indel=$SPLICEAI_INDEL")
  echo "[vep] SpliceAI enabled "
fi
if [[ "${USE_LOFTEE:-0}" == "1" ]]; then
  
  for f in "$HUMAN_ANCESTOR" "$GERP_BW" "$LOFTEE_SQL"; do
    [[ -f "$f" ]] || { echo "error: LOFTEE file missing: $f"; exit 1; }
  done
  command -v samtools >/dev/null || { echo "error: samtools not on PATH (conda install -n vep -c bioconda samtools)"; exit 1; }
  export PERL5LIB="${PERL5LIB:-}:$LOFTEE_PATH"      # so VEP finds loftee's modules
  VEP_ARGS+=(--plugin "LoF,loftee_path:$LOFTEE_PATH,human_ancestor_fa:$HUMAN_ANCESTOR,gerp_bigwig:$GERP_BW,conservation_file:$LOFTEE_SQL")
  echo "[vep] LOFTEE enabled "
fi
 
vep "${VEP_ARGS[@]}"
tabix -p vcf "$OUT"
echo "[vep] done -> $OUT"
 
 