#!/usr/bin/env bash
# 01_prep.sh  —  Get the RAW VCF ready to annotate 
#   (1) subset of sample + main chromosomes (1-22, X, Y).
#       Scaffolds/alt contigs and the mitochondrion do not contain MVA panel genes.
#   (2) reconciles the naming
#   (3) normalizes with 'bcftools norm':
#         -m-any   splits multiallelics (ALT "A,T") into one row per allele
#                  -> each allele can be matched separately (essential to phase
#                  a compound-het).
#         -f FASTA left-aligns indels -> the SAME variant is always written the
#                  same way as in ClinVar -> matching works.
# =====================================================================
set -euo pipefail

# edit paths
VCF_IN=""
FASTA="$HOME/mva/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTDIR="$HOME/mva/work"

mkdir -p "$OUTDIR"; cd "$OUTDIR"

#check parametres
[[ -n "$VCF_IN" ]] || { echo "ERROR: set VCF_IN"; exit 1; }
[[ -f "$FASTA"  ]] || { echo "ERROR: FASTA not found: $FASTA"; exit 1; }
#index
if [[ ! -f "${VCF_IN}.tbi" && ! -f "${VCF_IN}.csi" ]]; then
  echo "[prep] no index -> creating ${VCF_IN}.tbi"
  bcftools index -t -f "$VCF_IN"
fi

# sample
SAMPLE=$(bcftools query -l "$VCF_IN" | head -n1)
echo "[prep] sample = $SAMPLE"

#renaming if necessary 
FIRST_CHR=$(bcftools view -H "$VCF_IN" 2>/dev/null | head -n1 | cut -f1)
echo "[prep] first chromosome = $FIRST_CHR"

if [[ "$FIRST_CHR" == chr* ]]; then
  echo "[prep] 'chr' detected -> restricting to chr1-22,X,Y and stripping it for annotation"
  REGIONS=$(for i in $(seq 1 22) X Y; do printf 'chr%s,' "$i"; done | sed 's/,$//')
  : > rename.txt
  for i in $(seq 1 22) X Y; do printf 'chr%s\t%s\n' "$i" "$i" >> rename.txt; done
  RENAME=(bcftools annotate --rename-chrs rename.txt)
  echo 1 > .stripped_chr          # informational
else
  echo "[prep] no 'chr': restricting to 1-22,X,Y"
  REGIONS=$(for i in $(seq 1 22) X Y; do printf '%s,' "$i"; done | sed 's/,$//')
  RENAME=(cat)                    # passthrough: already matches the FASTA
  echo 0 > .stripped_chr
fi

#control
N_RAW=$(bcftools view -H "$VCF_IN" | wc -l)
N_MAIN=$(bcftools view -H -r "$REGIONS" "$VCF_IN" | wc -l)
echo "[prep] raw variants              : $N_RAW"
echo "[prep] on main chromosomes       : $N_MAIN   (dropped $((N_RAW - N_MAIN)) from scaffolds/mito, no panel genes)"

#Subset sample+regions -> rename -> normalize
bcftools view -s "$SAMPLE" -r "$REGIONS" "$VCF_IN" \
  | "${RENAME[@]}" \
  | bcftools norm -m-any -f "$FASTA" -Oz -o all.norm.vcf.gz 2> norm.log
bcftools index -t -f all.norm.vcf.gz

N_OUT=$(bcftools view -H all.norm.vcf.gz | wc -l)
echo "[prep] variants after normalize  : $N_OUT   (rises vs $N_MAIN: 'norm -m-any' splits multiallelics)"

# loss alert
if grep -qi "mismatch\|skipp" norm.log; then
  echo "[prep] norm reported REF-mismatch/skipped. Summary:"
  grep -i "mismatch\|skipp\|total" norm.log | tail -n 5
  echo "[prep]    -> usually means the wrong-build FASTA or a naming mismatch."
  echo "[prep]    -> check the FASTA is GRCh38 Ensembl. Full log: $OUTDIR/norm.log"
else
  echo "[prep] OK: norm reported no REF-mismatch"
fi

echo "[prep] done -> $OUTDIR/all.norm.vcf.gz"