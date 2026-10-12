#!/usr/bin/env bash
# =====================================================================
# 01_prep.sh  —  
# ---------------------------------------------------------------------
#  (1) RECONCILIA el naming de cromosomas. Si el VCF usa 'chr1' (estilo
#      Broad/GATK) lo pasa a '1' (estilo Ensembl), que es lo que esperan
#      el cache de VEP, ClinVar y AlphaMissense. Si no coincide, el cruce
#      falla SIN avisar y perderias la variante causal. Al final del
#      pipeline se re-anade 'chr' porque la plantilla pide 'chr15'.
#  (2) NORMALIZA con 'bcftools norm':
#        -m-any   parte los sitios multialelicos (ej. ALT "A,T") en una
#                 fila por alelo -> cada alelo se puede cruzar por separado.
#        -f FASTA alinea los indels a la izquierda (left-align): la misma
#                 variante queda siempre escrita igual -> el matching cuadra.
#
# NO filtra por PASS: conservamos tambien las no-PASS, porque los mosaicos
# (VAF baja) suelen caer ahi. 
#
# Requiere: bcftools (env 'mva') y el FASTA Ensembl GRCh38 (sin 'chr').
# =====================================================================
set -euo pipefail

# ============ EDITA ESTAS RUTAS ============
VCF_IN=""
FASTA="$HOME/mva/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTDIR="$HOME/mva/work"
# ===========================================

mkdir -p "$OUTDIR"; cd "$OUTDIR"

# Muestra: la unica del VCF (autodetectada).
SAMPLE=$(bcftools query -l "$VCF_IN" | head -n1)
echo "[prep] muestra = $SAMPLE"

# ¿El VCF usa prefijo 'chr'?
FIRST_CHR=$(bcftools view -H "$VCF_IN" | head -n1 | cut -f1)
echo "[prep] primer cromosoma = $FIRST_CHR"

if [[ "$FIRST_CHR" == chr* ]]; then
  echo "[prep] 'chr' detectado -> lo quito para anotar (se re-anade al final)"
  # Mapa chr1->1 ... chrX->X, chrY->Y, chrM->MT (Ensembl llama 'MT' a la mitocondria).
  : > chr2nochr.txt
  for i in $(seq 1 22) X Y; do printf 'chr%s\t%s\n' "$i" "$i" >> chr2nochr.txt; done
  printf 'chrM\tMT\n' >> chr2nochr.txt
  RENAME=(bcftools annotate --rename-chrs chr2nochr.txt)
  echo 1 > .stripped_chr          # marca para que 05_make_csv.py vuelva a poner 'chr'
else
  echo "[prep] sin 'chr': ya es estilo Ensembl"
  RENAME=(cat)                    # passthrough: no renombra nada
  echo 0 > .stripped_chr
fi

# Subset de la muestra  ->  renombrar contigs  ->  normalizar  ->  comprimir+indexar.
bcftools view -s "$SAMPLE" "$VCF_IN" \
  | "${RENAME[@]}" \
  | bcftools norm -m-any -f "$FASTA" -Oz -o all.norm.vcf.gz
bcftools index -t all.norm.vcf.gz

echo "[prep]  -> $OUTDIR/all.norm.vcf.gz"
echo "[prep] variantes tras normalizar: $(bcftools view -H all.norm.vcf.gz | wc -l)"
