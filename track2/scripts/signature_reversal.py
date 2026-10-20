# transcriptomic reversion (CMap / LINCS L1000 connectivity).
# for each drug, compare its expression z-scores over the disease-signature genes
# against the disease direction (+1 up / -1 down) with a cosine:
#     connectivity = cos(drug, disease)   (+1 mimics the disease, -1 reverses it)
#     reversal = -connectivity            (higher is better: a drug that reverses it)
#
# Input  : config paths disease_signature + drug_signatures
# Output : out/signature_reversal.tsv
# Usage  : python3 signature_reversal.py [config_2.yaml]
import sys, os, math, yaml
from collections import defaultdict

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config_2.yaml"))
P = cfg["paths"]
SIG = cfg.get("signature", {})
MIN_OVERLAP = int(SIG.get("min_overlap", 5))
SIG_CONN    = float(SIG.get("significant_conn", 0.20))
OUT = os.path.expanduser(P.get("out", "./out")); os.makedirs(OUT, exist_ok=True)

# refuse the synthetic demo signatures unless allowed, so real L1000 is used
_sigpath = os.path.expanduser(P["drug_signatures"])
if not cfg.get("allow_demo", False):
    try:
        if "SYNTHETIC DEMO" in "".join(open(_sigpath).readlines()[:12]):
            sys.exit("[sig] refusing the demo signatures. Build real ones:\n"
                     "  python3 drug_signatures.py lincs_up.gmt lincs_down.gmt\n"
                     "  then set paths.drug_signatures: ./data/drug_signatures_real.tsv "
                     "(or allow_demo: true).")
    except OSError:
        pass


def rows_of(path):
    # yield tab-split rows, skipping comments and the header
    for ln in open(os.path.expanduser(path)):
        ln = ln.rstrip("\n")
        if not ln or ln.lstrip().startswith("#"):
            continue
        parts = ln.split("\t")
        if parts and parts[0].strip().lower() in ("gene", "drug"):
            continue
        yield [p.strip() for p in parts]


# disease signature: gene -> +1 (up) / -1 (down)
disease = {}
for r in rows_of(P["disease_signature"]):
    if len(r) < 2:
        continue
    disease[r[0].upper()] = 1.0 if r[1].lower() in ("1", "+1", "up", "+") else -1.0
n_up = sum(1 for v in disease.values() if v > 0)
print(f"[sig] disease signature: {len(disease)} genes ({n_up} up, {len(disease) - n_up} down)")

# drug signatures: drug -> {gene: z}
drug_z = defaultdict(dict)
for r in rows_of(P["drug_signatures"]):
    if len(r) < 3:
        continue
    try:
        drug_z[r[0]][r[1].upper()] = float(r[2])
    except ValueError:
        continue
print(f"[sig] drug signatures: {len(drug_z)} drugs")


def connectivity(zmap):
    # reversal = -cosine(drug vector, disease direction) over overlapping genes
    genes = [g for g in disease if g in zmap]
    if len(genes) < MIN_OVERLAP:
        return None
    d = [disease[g] for g in genes]
    x = [zmap[g] for g in genes]
    dot = sum(di * xi for di, xi in zip(d, x))
    nd = math.sqrt(sum(di * di for di in d))
    nx = math.sqrt(sum(xi * xi for xi in x))
    if nd == 0 or nx == 0:
        return None
    cos = dot / (nd * nx)
    up_rev = sum(1 for g in genes if disease[g] > 0 and zmap[g] < 0)   # up-genes pushed down
    dn_rev = sum(1 for g in genes if disease[g] < 0 and zmap[g] > 0)   # down-genes pushed up
    return (-cos, len(genes), up_rev, dn_rev)


rows = []
for drug, zmap in drug_z.items():
    res = connectivity(zmap)
    if res is None:
        rows.append(dict(drug=drug, reversal=float("nan"), n_overlap=len(set(zmap) & set(disease)),
                         up_reversed=0, down_reversed=0,
                         call="UNSCORABLE", reason=f"<{MIN_OVERLAP} signature genes measured"))
        continue
    reversal, n_ov, up_rev, dn_rev = res
    if reversal >= SIG_CONN:
        call = "REVERSES (candidate)"
        reason = f"reverses the signature (reversal={reversal:+.2f}); {up_rev} up down, {dn_rev} down up of {n_ov}"
    elif reversal <= -SIG_CONN:
        call = "MIMICS (avoid)"
        reason = f"mimics the signature (reversal={reversal:+.2f}); would deepen the stress state"
    else:
        call = "NEUTRAL"
        reason = f"no clear connectivity (reversal={reversal:+.2f})"
    rows.append(dict(drug=drug, reversal=round(reversal, 3), n_overlap=n_ov,
                     up_reversed=up_rev, down_reversed=dn_rev, call=call, reason=reason))

# strongest reversers first 
rows.sort(key=lambda r: (r["reversal"] == r["reversal"],
                         r["reversal"] if r["reversal"] == r["reversal"] else -9), reverse=True)

cols = ["rank", "drug", "call", "reversal", "n_overlap", "up_reversed", "down_reversed", "reason"]
dest = os.path.join(OUT, "signature_reversal.tsv")
with open(dest, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for i, r in enumerate(rows, 1):
        fh.write("\t".join([str(i)] + [str(r[c]) for c in cols[1:]]) + "\n")

print(f"[sig] ranked {len(rows)} drugs -> {dest}")
for i, r in enumerate(rows[:15], 1):      # show only the top 15
    print(f"  #{i:<2} {r['drug']:<16} {r['call']:<22} reversal={r['reversal']}")
