# Build drug signatures
# Output (long format): drug <tab> gene <tab> zscore  (+1 up-gene, -1 down-gene).
#
# Get the data (public; wget works despite robots.txt):
#   base='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName='
#   wget -O lincs_up.gmt   "${base}LINCS_L1000_Chem_Pert_up"
#   wget -O lincs_down.gmt "${base}LINCS_L1000_Chem_Pert_down"
# GMT line: term <tab> description <tab> gene1 <tab> gene2 ...  The drug name is
# taken from the term; the script prints a few samples so you can check parsing.
#
# Usage:
#   python3 build_drug_signatures.py lincs_up.gmt lincs_down.gmt
#   python3 build_drug_signatures.py up.gmt down.gmt --out data/drug_signatures_real.tsv
#   python3 build_drug_signatures.py up.gmt down.gmt --name-regex '^([a-zA-Z0-9\-]+)'
import sys, os, re, argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("up_gmt")
ap.add_argument("down_gmt")
ap.add_argument("--out", default=os.path.join("data", "drug_signatures_real.tsv"))
ap.add_argument("--name-regex", default=r"^\s*([A-Za-z0-9][A-Za-z0-9+\-]*)",
                help="group(1) is the drug label inside the GMT term (default: first token)")
args = ap.parse_args()
rx = re.compile(args.name_regex)


def drug_label(term):
    # pull the drug name out of the perturbation term
    m = rx.search(term)
    return (m.group(1) if m else term).strip().lower()


# sig[(drug, gene)] accumulates +1/-1 across perturbations; averaged at the end
sig = defaultdict(float)
cnt = defaultdict(int)
samples = []


def ingest(path, value):
    # add every (drug, gene) from one GMT with the given direction (+1 up, -1 down)
    terms = 0
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        term = parts[0]
        drug = drug_label(term)
        if not drug:
            continue
        if len(samples) < 8:
            samples.append((term[:50], drug))
        genes = [g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()]
        for g in genes:
            sig[(drug, g)] += value
            cnt[(drug, g)] += 1
        terms += 1
    return terms


tu = ingest(args.up_gmt, +1.0)
td = ingest(args.down_gmt, -1.0)
print(f"[sig] terms read: up={tu}, down={td}")
print("[sig] sample term -> drug label (check these look like drug names):")
for term, drug in samples:
    print(f"   {term:<52} -> {drug}")

drugs = sorted({d for (d, _g) in sig})
os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
rows = 0
with open(args.out, "w") as fh:
    fh.write("drug\tgene\tzscore\n")
    for (drug, g), s in sig.items():
        z = s / cnt[(drug, g)]           # mean direction across this drug's perturbations
        if z == 0:                       # up and down cancelled: ambiguous, drop
            continue
        fh.write(f"{drug}\t{g}\t{z:.3f}\n")
        rows += 1

print(f"[sig] wrote {rows} rows for {len(drugs)} drugs -> {args.out}")
