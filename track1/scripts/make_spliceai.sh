set -euo pipefail

WORK="$HOME/mva/work"
DATA="$HOME/mva/data"
FASTA="$HOME/mva/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
WINDOW=50000          # +-bp around each candidate position to score

CAND="$WORK/candidates.tsv"
NORM="$WORK/all.norm.vcf.gz"
BED="$DATA/cand.bed"
SUB="$DATA/cand.vcf"
OUT="$DATA/spliceai.hg38.vcf"

command -v spliceai >/dev/null || { echo "error: spliceai not installed -> pip install spliceai tensorflow"; exit 1; }
command -v bcftools >/dev/null || { echo "error: bcftools not found (conda activate mva)"; exit 1; }
[[ -f "$CAND"  ]] || { echo "error: $CAND not found -> run 01->02(no plugin)->03 first"; exit 1; }
[[ -f "$NORM"  ]] || { echo "error: $NORM not found -> run 01_prep.sh first"; exit 1; }
[[ -f "$FASTA" ]] || { echo "error: FASTA not found: $FASTA"; exit 1; }
[[ -f "${FASTA}.fai" ]] || { echo "[spliceai] indexing FASTA..."; samtools faidx "$FASTA"; }
mkdir -p "$DATA"

awk -F'\t' -v w="$WINDOW" 'NR>1 && $8!="" {s=$9-w; if(s<0)s=0; print $8"\t"s"\t"$9+w}' \
    "$CAND" | sort -k1,1 -k2,2n \
  | awk 'BEGIN{OFS="\t"} {if(NR==1){c=$1;s=$2;e=$3} else if($1==c&&$2<=e){if($3>e)e=$3} else {print c,s,e;c=$1;s=$2;e=$3}} END{if(NR)print c,s,e}' \
    > "$BED"
echo "[spliceai] regions to score: $(wc -l < "$BED") merged windows -> $BED"

bcftools view -R "$BED" "$NORM" -Ov -o "$SUB"
echo "[spliceai] variants to score: $(grep -vc '^#' "$SUB")"

# run SpliceAI 
echo "[spliceai] running SpliceAI (this is the slow step; a few thousand variants = minutes)..."
spliceai -I "$SUB" -O "$OUT" -R "$FASTA" -A grch38 -D "$((WINDOW/1000))"

bcftools sort "$OUT" -Oz -o "${OUT}.gz"
tabix -p vcf "${OUT}.gz"
rm -f "$OUT" "$SUB"
echo "[spliceai] done -> ${OUT}.gz"

