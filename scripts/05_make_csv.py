#!/usr/bin/env python3
# INPUT : work/candidates.tsv  
# OUTPUT: out/<hf_username>_<approach_name>.csv  
#   1. Takes the best candidates
#   2. Assigns epcr .
#   4. Marks finding_type: primary (panel gene) / secondary (incidental).
#
# USAGE:  python3 05_make_csv.py            # uses ./config.yaml
import sys, os, csv, yaml

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config.yaml"))
WORK = os.path.expanduser(cfg["paths"]["work"])
OUT  = os.path.expanduser(cfg["paths"]["out"]); os.makedirs(OUT, exist_ok=True)
CAND = os.path.join(WORK, "candidates.tsv")

# Optional keys with defaults (so it does not fail if they are missing in config.yaml).
PID      = cfg.get("proband_id", "PROBAND")
HF       = cfg.get("hf_username", "user")
APPROACH = cfg.get("approach_name", "offline-3track")
PREFIX   = "chr" if cfg.get("add_chr_prefix", True) else ""
E = cfg.get("epcr", {"top_clinvar_or_biallelic_lof": 0.95, "strong_contender": 0.55,
                     "weak_contender": 0.20, "secondary_floor": 0.10})
fname = f"{HF}_{APPROACH}.csv"
dest = os.path.join(OUT, fname)

# read candidates
rows = list(csv.DictReader(open(CAND), delimiter="\t"))

#2nd-copy findings they flag a gene as probably biallelic
loh_genes, loh_dens = set(), {}
loh_path = os.path.join(OUT, "loh_scan.tsv")
if os.path.exists(loh_path):
    for lr in csv.DictReader(open(loh_path), delimiter="\t"):
        if "SUSPECTED_BIALLELIC" in (lr.get("flag", "")):
            loh_genes.add(lr["gene"]); loh_dens[lr["gene"]] = lr.get("rel_density", "")

def strong_evidence(r):
    clnsig = (r["clnsig_1"] or "").lower()
    is_clinvar = "pathogenic" in clnsig and "conflicting" not in clnsig
    biallelic_core = r["kind"] in ("homozygous", "compound_het") and r["tier"] == "core"
    loh_biallelic = r["symbol"] in loh_genes      # 2nd copy detected by LOH => biallelic
    return is_clinvar or biallelic_core or loh_biallelic


PRIMARY_TIERS = ("core", "ext", "phenocopy", "pheno_hpo", "core_rescue")
TIER_RANK = {"core": 5, "ext": 4, "phenocopy": 3, "pheno_hpo": 2, "core_rescue": 1}

def is_path(r):
    s = (r["clnsig_1"] or "").lower()
    return "pathogenic" in s and "conflicting" not in s

#   1) tier (core > ext > phenocopy > HPO > rescue)  
#   2) 2nd copy (LOH) = biallelic  -> only rises within the same tier
#   3) ClinVar pathogenic
#   4) score
def prim_key(r):
    return (TIER_RANK.get(r["tier"], 0), r["symbol"] in loh_genes, is_path(r), float(r["score"]))
primaries   = sorted([r for r in rows if r["tier"] in PRIMARY_TIERS], key=prim_key, reverse=True)
secondaries = sorted([r for r in rows if r["tier"] not in PRIMARY_TIERS],
                     key=lambda r: (r["tier"] == "acmg_sf", is_path(r), float(r["score"])), reverse=True)
MAX_ROWS, RESERVE_SEC = 10, 3
n_prim = min(len(primaries), MAX_ROWS - min(len(secondaries), RESERVE_SEC))
n_sec  = min(len(secondaries), MAX_ROWS - n_prim)
selected = primaries[:n_prim] + secondaries[:n_sec]

out_rows = []
prev_epcr = 1.0          # ordered by confidence
for i, r in enumerate(selected):
    primary = r["tier"] in PRIMARY_TIERS
    finding = "primary" if primary else "secondary"

    if i == 0:
        epcr = E["top_clinvar_or_biallelic_lof"] if strong_evidence(r) else E["strong_contender"]
    elif i <= 2:
        epcr = E["strong_contender"] if strong_evidence(r) else E["weak_contender"]
    else:
        epcr = E["weak_contender"]
    if r["tier"] in ("phenocopy", "pheno_hpo"):      
        cap = E["strong_contender"] if r["symbol"] in loh_genes else E["weak_contender"]
        epcr = min(epcr, cap)
    if r["tier"] == "core_rescue":                   
        epcr = min(epcr, 0.15)
    if finding == "secondary":
        epcr = min(epcr, E["secondary_floor"])       # <=0.1
    # the ranking already goes from more to less confidence -> epcr never rises going down
    epcr = min(epcr, prev_epcr)
    prev_epcr = epcr

    #justification note 
    if finding == "secondary":
        zyg = {"homozygous": "homozygous (biallelic)", "compound_het": "compound-het (biallelic)"}.get(r["kind"], "heterozygous (carrier if the gene is recessive)")
        act = "ACMG-actionable" if r["tier"] == "acmg_sf" else "ClinVar pathogenic / non-actionable gene"
        note = (f"secondary finding ({act}); zygosity={zyg}; " f"gene={r['symbol']}; ClinVar={r['clnsig_1']}; cons={r['consequence_1']}")
    else:
        bits = [r["kind"], f"gene={r['symbol']}"]
        if r["clnsig_1"]: bits.append(f"ClinVar={r['clnsig_1']}")
        if r["am_1"]:     bits.append(f"AM={r['am_1']}")
        if r["mosaic"] == "True": bits.append("low-VAF(mosaic)")
        if r["symbol"] in loh_genes:
            bits.append(f"probable 2nd copy (LOH, rel_density={loh_dens.get(r['symbol'], '?')})")
        if r.get("review"):
            bits.append(r["review"])
        bits.append(f"cons={r['consequence_1']}")
        note = "; ".join(bits)

    #second variant 
    c2 = p2 = rf2 = al2 = ""
    if r["chrom_2"]:
        c2, p2, rf2, al2 = PREFIX + r["chrom_2"], r["pos_2"], r["ref_2"], r["alt_2"]

    out_rows.append(dict(
        proband_id=PID,
        chrom_1=PREFIX + r["chrom_1"], pos_1=r["pos_1"], ref_1=r["ref_1"], alt_1=r["alt_1"],
        chrom_2=c2, pos_2=p2, ref_2=rf2, alt_2=al2,
        epcr=round(float(epcr), 3), finding_type=finding, notes=note,
    ))

cols = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes"]
with open(dest, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(out_rows)

print(f"[csv] {len(out_rows)} rows -> {dest}")
for r in out_rows:
    print(f"   {r['chrom_1']}:{r['pos_1']} {r['ref_1']}>{r['alt_1']}"
          f"{'  + '+r['chrom_2']+':'+str(r['pos_2']) if r['chrom_2'] else ''}"
          f"  epcr={r['epcr']}  {r['finding_type']}")