#!/usr/bin/env python3
# =====================================================================
# 00_inspect.py 
# ---------------------------------------------------------------------
#   Antes de anotar/filtrar necesitamos conocer el VCF. Este script lee
#   solo la ESTRUCTURA (cabecera + columnas) e imprime un resumen. Ninguna
#   linea de salida contiene genotipos ni coordenadas: son propiedades de
#   formato (naming, build, campos, conteos), que SI se pueden compartir.
#
#   1. Nombre de la muestra  -> es el proband_id del CSV de entrega.
#   2. Prefijo 'chr' (chr1 vs 1) -> hay que reconciliarlo con VEP/ClinVar
#      (estilo Ensembl, sin 'chr') o el cruce falla EN SILENCIO.
#   3. Build (longitud de chr1)  -> el reto exige GRCh38.
#   4. Campos GT/AD/DP  -> sin AD/DP no hay VAF, ni mosaico, ni mapa BAF.
#   5. Distribucion de FILTER  -> si es '.', NO se puede filtrar por PASS.
#   6. % indels y multialelicos -> justifican 'bcftools norm'.
#
# USO:   python3 00_inspect.py /ruta/al/xxx.vcf.gz
# =====================================================================
import sys, gzip, collections

if len(sys.argv) != 2:
    sys.exit("Uso: python3 00_inspect.py <vcf.gz>")
vcf = sys.argv[1]

# Longitud de chr1 en cada build, para distinguir GRCh38 de GRCh37.
CHR1_LEN = {248956422: "GRCh38", 249250621: "GRCh37/hg19"}

samples, caller, info_ids, format_fields = [], set(), set(), set()
filters = collections.Counter()
chr1_len = first_chrom = None
n_sites = n_indel_alleles = n_multiallelic = 0

# gzip.open(..., "rt") lee el .vcf.gz directamente (texto), sin descomprimir a disco.
with gzip.open(vcf, "rt") as fh:
    for line in fh:

        if line.startswith("#"): # cabecera (metadatos + nombres de columnas)
            low = line.lower()
            # Las lineas de contig tienen esta forma:  ##contig=<ID=chr1,length=248956422>
            # De ahi sacamos la longitud de chr1 -> build.
            if line.startswith("##contig") and ("id=chr1," in low or "id=1," in low):
                try:
                    chr1_len = int(line.split("length=")[1].split(">")[0].split(",")[0])
                except Exception:
                    pass
            elif line.startswith("##INFO=<ID="):          # IDs de campos INFO disponibles
                info_ids.add(line.split("ID=")[1].split(",")[0])
            for t in ("gatk", "haplotypecaller", "dragen", "deepvariant", "strelka", "mutect"):
                if t in low:
                    caller.add(t)                          # programa que llamo las variantes
            if line.startswith("#CHROM"):                  # ultima linea de cabecera:
                samples = line.rstrip("\n").split("\t")[9:] # columnas 10+ = muestras
            continue

        # una variante por linea
        f = line.split("\t", 9)                            # no hace falta partir mas alla de FORMAT
        chrom, ref, alt, flt, fmt = f[0], f[3], f[4], f[6], f[8]
        if first_chrom is None:
            first_chrom = chrom
        filters[flt] += 1
        format_fields.update(fmt.split(":"))
        alts = alt.split(",")
        if len(alts) > 1:
            n_multiallelic += 1                            # sitio con >1 alelo alterado
        for a in alts:
            if len(ref) != len(a):
                n_indel_alleles += 1                       # insercion o delecion
        n_sites += 1

# informe 
chr_prefix = str(first_chrom).startswith("chr")
ad, dp = "AD" in format_fields, "DP" in format_fields
print(f"Muestra(s)           : {samples}   -> usar como proband_id")
print(f"Caller               : {', '.join(sorted(caller)) or 'no evidente'}")
print(f"Primer cromosoma     : {first_chrom}")
print(f"Prefijo 'chr'        : {'SI -> pasar a estilo Ensembl para anotar' if chr_prefix else 'NO -> ya estilo Ensembl'}")
build = CHR1_LEN.get(chr1_len, f"indeterminado (chr1={chr1_len})")
print(f"Build                : {build}   {'[OK, el reto pide GRCh38]' if chr1_len == 248956422 else '[!! el reto pide GRCh38]'}")
print(f"GT / AD / DP         : GT={'GT' in format_fields}  AD={ad}  DP={dp}   "
      f"{'(OK: hay VAF -> mosaico + mapa BAF)' if ad and dp else '(FALTA AD/DP: sin VAF ni mapa BAF)'}")
print(f"FILTER               : {dict(filters.most_common())}")
if set(filters) <= {"."}:
    print("   -> VCF sin filtrar (FILTER='.'): NO uses '-f PASS' o borras todo.")
print(f"Variantes totales    : {n_sites}")
print(f"Sitios multialelicos : {n_multiallelic}   -> por eso 'bcftools norm -m-any'")
print(f"Alelos indel         : {n_indel_alleles}   -> por eso left-align con -f FASTA")
print(f"INFO tags            : {sorted(info_ids)}")
print("=" * 64)
