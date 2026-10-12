#!/usr/bin/env bash
# =====================================================================
# 02_annotate.sh  —  Anota cada variante con VEP
# ---------------------------------------------------------------------
# x variante -> gen y consecuencia, frecuencia gnomAD,
# clasificacion ClinVar y score AlphaMissense.
# ningun dato del paciente sale de la maquina. // mirar
#
# NO ejecutar hasta tener: cache VEP + FASTA + ClinVar +
# AlphaMissense
#
# flags:
#   --offline --cache     modo local, sin ninguna llamada de red
#   --dir_cache / --fasta  donde estan el cache y la referencia
#   --tab --output_file    salida en tabla (facil de leer con pandas)
#   --fork 4               usa 4 nucleos (mas rapido en >500k variantes)
#   --symbol               nombre del gen (SYMBOL) -> por el filtramos el panel
#   --canonical            marca el transcrito canonico
#   --af_gnomadg --af_gnomade  frecuencias gnomAD (genomas/exomas) del cache
#   --sift --polyphen      predictores clasicos (de apoyo a AlphaMissense)
#   --pick_allele_gene     1 consecuencia por alelo y gen (salida manejable)
#   --custom ClinVar       cruza con ClinVar (CLNSIG = clasificacion clinica)
#   --plugin AlphaMissense anade el score de patogenicidad missense
# =====================================================================
set -euo pipefail

WORK="$HOME/mva/work"
CACHE="$HOME/mva/cache"
REF="$HOME/mva/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DATA="$HOME/mva/data"

vep \
  --offline --cache --dir_cache "$CACHE" --fasta "$REF" \
  --species homo_sapiens --assembly GRCh38 \
  --input_file "$WORK/all.norm.vcf.gz" --format vcf \
  --tab --output_file "$WORK/annotated.tsv" --force_overwrite \
  --fork 4 --buffer_size 20000 \
  --symbol --canonical \
  --af_gnomadg --af_gnomade \
  --sift b --polyphen b \
  --pick_allele_gene \
  --custom "$DATA/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN" \
  --plugin AlphaMissense,file="$DATA/AlphaMissense_hg38.tsv.gz"

# --- FASE 2 (opcional, pesan mucho; anadir dentro del comando de arriba) --- -> esto interesante
#   --plugin SpliceAI,snv=...,indel=...   # efecto en splicing (2o alelo tipico de BUB1B)
#   --plugin LoF,...                      # LOFTEE: confirma perdida de funcion

echo "[vep] listo -> $WORK/annotated.tsv"
