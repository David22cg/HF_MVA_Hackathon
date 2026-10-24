#!/usr/bin/env python3
# Global comparison of RankVar's ranking and this pipeline's ranking-> Benchmarking
# Inputs:
#   rank_var.tsv     RankVar output
#   candidates.tsv   output of our 03_filter_rank.py
#   [panel.txt]      optional, e.g. ../genes/mva_g.txt, to tag MVA genes
#
# Usage:
#   python3 compare_rankvar.py rank_var.tsv candidates.tsv ../genes/mva_g.txt
import sys, csv

RANKVAR = sys.argv[1] if len(sys.argv) > 1 else "rank_var.tsv"
OURS    = sys.argv[2] if len(sys.argv) > 2 else "candidates.tsv"
PANEL   = sys.argv[3] if len(sys.argv) > 3 else ""


def norm_chr(c):
    return str(c).replace("chr", "").replace("Chr", "").strip()


def split_genes(field):
    out = set()
    for part in str(field).replace(";", ",").replace("|", ",").split(","):
        g = part.strip().upper()
        if g and g not in (".", "NA", "NONE"):
            out.add(g)
    return out


def load_panel(path):
    d = {}
    if not path:
        return d
    try:
        for line in open(path):
            line = line.split("#")[0].strip()
            if line:
                p = line.split()
                d[p[0].upper()] = p[1] if len(p) > 1 else "ext"
    except OSError:
        pass
    return d


def pearson(xs, ys):
    # Pearson on the rank values equals Spearman of the two rankings
    n = len(xs)
    if n < 2:
        return None
    mx = sum(xs) / n; my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs); vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return None
    return cov / (vx * vy) ** 0.5


panel = load_panel(PANEL)
tag = lambda g: (panel.get(g, "off-panel") if panel else "?")


# RankVar: gene -> best rank, ordered genes, and a locus index for variant match
rv_rank, rv_rows, rv_locus = {}, [], {}
with open(RANKVAR) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    cols = {c.lower(): c for c in (r.fieldnames or [])}
    c_rank = cols.get("rank")
    
    c_gene = (cols.get("gene.refgene")
              or next((cols[k] for k in cols if k.startswith("gene") and "func" not in k), None))
    c_chr  = cols.get("chr") or cols.get("chromosome")
    c_pos  = cols.get("start") or cols.get("pos") or cols.get("position")
    if not c_rank:
        sys.exit("[compare] no 'Rank' column in the RankVar file. Columns: "
                 + ", ".join(r.fieldnames or []))
    for row in r:
        try:
            rank = float(row[c_rank])
        except (KeyError, ValueError, TypeError):
            continue
        genes = split_genes(row.get(c_gene, "")) if c_gene else set()
        rv_rows.append((rank, genes))
        for g in genes:
            if g not in rv_rank or rank < rv_rank[g]:
                rv_rank[g] = rank
        if c_chr and c_pos:
            key = (norm_chr(row.get(c_chr, "")), str(row.get(c_pos, "")).strip())
            if key not in rv_locus or rank < rv_locus[key]:
                rv_locus[key] = rank
rv_rows.sort(key=lambda x: x[0])
rv_top_genes = []
for _, genes in rv_rows:
    for g in genes:
        if g not in rv_top_genes:
            rv_top_genes.append(g)


# gene-> rank, tier , locus
our_top_genes, our_rank, our_tier, our_locus = [], {}, {}, {}
for i, row in enumerate(csv.DictReader(open(OURS), delimiter="\t"), 1):
    g = (row.get("symbol") or "").upper()
    if not g or g in our_rank:
        continue
    our_rank[g] = i
    our_tier[g] = row.get("tier", "")
    our_locus[g] = (norm_chr(row.get("chrom_1", "")), str(row.get("pos_1", "")).strip())
    our_top_genes.append(g)


shared = [g for g in our_rank if g in rv_rank]
rho = pearson([our_rank[g] for g in shared], [rv_rank[g] for g in shared])


def overlap(k):
    a = set(our_top_genes[:k]); b = set(rv_top_genes[:k])
    return len(a & b), sorted(a & b)

print("RankVar  comparison")
print(f"genes ranked: ours={len(our_rank)}  RankVar={len(rv_rank)}  shared={len(shared)}")
print("rank correlation on shared genes , spearman: "
      + (f"{rho:.2f}" if rho is not None else "n/a (too few shared genes)"))
for k in (5, 10, 20):
    n, g = overlap(k)
    print(f"top-{k:<2} overlap: {n}/{k}" + (f" -> {g}" if g else ""))

#  global distribution
buckets = {"top10": 0, "11-50": 0, "51+": 0, "absent": 0}
for g in our_top_genes:
    rr = rv_rank.get(g)
    if rr is None:
        buckets["absent"] += 1
    elif rr <= 10:
        buckets["top10"] += 1
    elif rr <= 50:
        buckets["11-50"] += 1
    else:
        buckets["51+"] += 1

print(f"where our {len(our_top_genes)} genes land in RankVar: "
      f"top10={buckets['top10']}  11-50={buckets['11-50']}  51+={buckets['51+']}  absent={buckets['absent']}")

# panel composition
if panel:
    top20 = rv_top_genes[:20]
    on = sum(1 for g in top20 if g in panel)
    print(f"RankVar top-20 composition: MVA-panel={on}/{len(top20)}  off-panel={len(top20) - on}")

# gene-level summary
print("ranking with RankVar position and variant check")
print(f"  {'gene':<12}{'ours':>6}{'RankVar':>9}  {'tier':<10} variant-in-RankVar")
for g in our_top_genes[:15]:
    rr = rv_rank.get(g)
    loc = our_locus.get(g)
    vfound = ""
    if rv_locus and loc and loc[1]:
        vr = rv_locus.get(loc)
        vfound = (f"yes (#{int(vr)})" if vr is not None else "ABSENT")
    print(f"  {g:<12}{('#'+str(our_rank[g])):>6}"
          f"{(('#'+str(int(rr))) if rr is not None else 'absent'):>9}  {our_tier.get(g,''):<10} {vfound}")

# variant-level summary over candidates
if rv_locus:
    present = absent = 0
    for g in our_top_genes[:15]:
        loc = our_locus.get(g)
        if not loc or not loc[1]:
            continue
        if rv_locus.get(loc) is not None:
            present += 1
        else:
            absent += 1
    print(f"variant-level: of our candidate variants, {present} are in RankVar, "
          f"{absent} are absent from its table")