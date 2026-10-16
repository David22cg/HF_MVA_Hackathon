#!/usr/bin/env python3
# USAGE:  python3 00b_panel_audit.py            
# uses ./config.yaml

import sys, os, gzip, yaml
from collections import defaultdict

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config.yaml"))
WORK = os.path.expanduser(cfg["paths"]["work"])
ANN  = os.path.join(WORK, "annotated.vcf.gz")
MAX_AF_REC = cfg["freq"]["max_af_recessive"]
AM_MIN     = cfg["deleterious"]["alphamissense_min"]
HIGH_IMP   = set(cfg["deleterious"]["high_impact"])


def load_panel(path, default_tier):
    #Returns {gene: tier} 
    d = {}
    p = os.path.expanduser(path or "")
    if p and os.path.exists(p):
        for line in open(p):
            line = line.split("#")[0].strip()
            if line:
                parts = line.split()
                d[parts[0]] = parts[1] if len(parts) > 1 else default_tier
    return d


#MVA panel (core/ext) + phenocopies
panel = load_panel(cfg["paths"]["panel"], "ext")
for g, t in load_panel(cfg["paths"].get("phenocopy"), "phenocopy").items():
    panel.setdefault(g, "phenocopy")

#per-gene each key is a funnel stage
stage = defaultdict(lambda: dict(total=0, in_sample=0, rare=0, strict=0, moderate=0))


def af_float(x):
    try: return float(x)
    except (ValueError, TypeError): return 0.0


def in_sample(fmt, sample):
    #Does the sample carry the allele 
    d = dict(zip(fmt.split(":"), sample.split(":")))
    gt = d.get("GT", "./.").replace("|", "/")
    return "1" in gt.split("/")


csq_fields = []
with gzip.open(ANN, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = line.split("Format:")[1].strip().strip('">').split("|")
            continue
        if line.startswith("#CHROM"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 10:
            continue
        info, fmt, sample = c[7], c[8], c[9]
        csq = next((f[4:] for f in info.split(";") if f.startswith("CSQ=")), "")
        if not csq:
            continue
        carried = in_sample(fmt, sample) # does the sample carry this variant
        seen_here = set()
        for e in (dict(zip(csq_fields, part.split("|"))) for part in csq.split(",")):
            sym = e.get("SYMBOL", "")
            if sym not in panel or sym in seen_here:
                continue
            seen_here.add(sym)
            st = stage[sym]
            st["total"] += 1                        #annotated to this gene whether or not the sample carries it
            if not carried:
                continue
            st["in_sample"] += 1                    #the sample carries this variant

            af = max(af_float(e.get("gnomADg_AF")), af_float(e.get("gnomADe_AF")))
            if af <= MAX_AF_REC:
                st["rare"] += 1                     #rare in population (recessive requirement)

            cons = e.get("Consequence", "")
            impact = e.get("IMPACT", "")
            am = e.get("am_pathogenicity", "")
            am = float(am) if am not in ("", ".") else None
            clnsig = (e.get("ClinVar_CLNSIG", "") or "").lower()
            is_lof = impact == "HIGH" or any(x in cons for x in HIGH_IMP)
            is_mis_path = am is not None and am >= AM_MIN
            is_clinvar_path = "pathogenic" in clnsig and "conflicting" not in clnsig
            if af <= MAX_AF_REC and (is_lof or is_mis_path or is_clinvar_path):
                st["strict"] += 1                   # rare + strongly damaging 
            elif af <= MAX_AF_REC and am is not None and am >= 0.34:
                st["moderate"] += 1                 # rare + VUS (possible 2nd copy)

#report
order = sorted(panel, key=lambda g: ({"core": 0, "ext": 1}.get(panel[g], 2), g))
print(f"panel audit over {ANN}")
print("Funnel per gene:  total annotated -> in sample -> rare -> {strong | VUS}")
print(f"{'gene':<10}{'tier':<11}{'total':>7}{'sample':>9}{'rare':>7}{'strong':>9}{'VUS':>6}")
for g in order:
    st = stage.get(g, dict(total=0, in_sample=0, rare=0, strict=0, moderate=0))
    flag = ""
    if st["total"] == 0:
        flag = "  <- 0 annotated: check prep/names/SYMBOL"
    print(f"{g:<10}{panel[g]:<11}{st['total']:>7}{st['in_sample']:>9}"
          f"{st['rare']:>7}{st['strict']:>9}{st['moderate']:>6}{flag}")

# per-tier summary
print("-" * 78)
for tier in ("core", "ext", "phenocopy"):
    gs = [g for g in panel if panel[g] == tier]
    seen = sum(1 for g in gs if stage.get(g, {}).get("total", 0) > 0)
    strong = sum(stage.get(g, {}).get("strict", 0) for g in gs)
    print(f"{tier:<11}: {len(gs)} genes | with >=1 annotated variant: {seen} | strong candidates: {strong}")

print("note: if the ext genes show 'total'>0 but 'strong'=0, nothing was lost in prep.")
print("the patient prob has no rare+damaging variant in those genes ")