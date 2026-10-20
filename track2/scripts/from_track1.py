#!/usr/bin/env python3
# Reads Track 1's ranked candidates.tsv and writes the primary gene symbols to data/seed_genes.txt
#
# Usage:
#   python3 from_track1.py /path/to/track1/work/candidates.tsv
#   python3 from_track1.py            # uses paths.track1_candidates from config
import sys, os, csv, yaml

cfg = {}
if os.path.exists("config_2.yaml"):
    cfg = yaml.safe_load(open("config_2.yaml")) or {}
P = cfg.get("paths", {})

src = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(P.get("track1_candidates", ""))
if not src or not os.path.exists(src):
    sys.exit("usage: python3 from_track1.py <track1 work/candidates.tsv>\n"
             "   (or set paths.track1_candidates in config_2.yaml)")

OUT = os.path.expanduser(P.get("seed_genes", "./data/seed_genes.txt"))
# tiers that are disease-relevant (not secondary/incidental findings)
PRIMARY = {"core", "ext", "phenocopy", "pheno_hpo", "core_rescue"}

seen, seeds = set(), []
for r in csv.DictReader(open(src), delimiter="\t"):
    sym = (r.get("symbol") or "").strip().upper()
    tier = (r.get("tier") or "").strip()
    try:
        score = float(r.get("score", "0"))
    except ValueError:
        score = 0.0
    if not sym or sym in seen:
        continue
    if tier in PRIMARY:
        seeds.append((sym, tier, score))
        seen.add(sym)

seeds.sort(key=lambda x: x[2], reverse=True)   # best-scoring gene first

os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("# disease seed gene(s) from Track 1, best score first (generated; re-run to refresh)\n")
    for sym, tier, score in seeds:
        fh.write(f"{sym}\t# {tier}, track1_score={score}\n")

print(f"[bridge] read {src}")
print(f"[bridge] {len(seeds)} primary gene(s) -> {OUT}")
for sym, tier, score in seeds:
    print(f"   {sym:<10} ({tier}, score={score})")
if not seeds:
    print("[bridge] warning: no primary genes found; Track 2 falls back to the static checkpoint list")
