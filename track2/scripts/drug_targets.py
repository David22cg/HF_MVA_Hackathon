#!/usr/bin/env python3
# Build a real drug->target table from the DGIdb interactions export
# Output: drug <tab> targets(; sep) <tab> status <tab> pediatric <tab> note
# Get the data (public): the interactions TSV from https://www.dgidb.org/downloads
# Usage:
#   python3 build_drug_targets.py interactions.tsv           
#   python3 build_drug_targets.py interactions.tsv out.tsv
import sys, csv, os
from collections import defaultdict

if len(sys.argv) < 2:
    sys.exit("usage: python3 build_drug_targets.py <DGIdb interactions.tsv> [out.tsv]")
SRC = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join("data", "drug_targets_real.tsv")

with open(SRC) as fh:
    all_lines = fh.readlines()
hdr_idx = next((i for i, ln in enumerate(all_lines)
                if ln.strip() and not ln.lstrip().startswith("#")), None)
if hdr_idx is None:
    sys.exit(f"[targets] error: no header row in {SRC} (only comments/blank lines?)")
head = all_lines[hdr_idx]
data_lines = all_lines[hdr_idx + 1:]
delim = "\t" if head.count("\t") >= head.count(",") else ","
cols = [c.strip() for c in head.rstrip("\n").split(delim)]
low = [c.lower() for c in cols]
if hdr_idx > 0:
    print(f"[targets] skipped {hdr_idx} preamble line(s) before the header")


def pick(cands):
    # first exact column-name match, then substring match
    for want in cands:
        for i, c in enumerate(low):
            if c == want:
                return cols[i]
    for want in cands:
        for i, c in enumerate(low):
            if want in c:
                return cols[i]
    return None


drug_col = pick(["drug_name", "drug_claim_name", "drug"])
gene_col = pick(["gene_name", "gene_claim_name", "gene"])
appr_col = pick(["approved", "drug_approved", "approval"])
print(f"[targets] columns -> drug: {drug_col!r}  gene: {gene_col!r}  approval: {appr_col!r}")
if not drug_col or not gene_col:
    sys.exit(f"[targets] error: could not find drug/gene columns in: {cols}")

targets = defaultdict(set)
approved = {}
n = 0
for r in csv.DictReader(data_lines, fieldnames=cols, delimiter=delim):
    drug = (r.get(drug_col) or "").strip()
    gene = (r.get(gene_col) or "").strip().upper()
    if not drug or not gene:
        continue
    targets[drug].add(gene)
    if appr_col:
        val = (r.get(appr_col) or "").strip().lower()
        approved[drug] = approved.get(drug, False) or val in ("true", "1", "yes", "approved")
    n += 1

os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("drug\ttargets\tstatus\tpediatric\tnote\n")
    for drug in sorted(targets):
        tg = ";".join(sorted(targets[drug]))
        status = "approved" if approved.get(drug) else "investigational"
        # pediatric safety is not in DGIdb; leave blank and curate the top hits by hand
        fh.write(f"{drug}\t{tg}\t{status}\t\tDGIdb\n")

print(f"[targets] {n} interactions -> {len(targets)} drugs -> {OUT}")
print(f"[targets] approved flagged: {sum(1 for v in approved.values() if v)} "
      "(pediatric left blank; curate the top hits)")
