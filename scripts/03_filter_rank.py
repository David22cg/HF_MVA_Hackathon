#!/usr/bin/env python3
# INPUT : work/annotated.vcf.gz  (from 02 = VCF with VEP annotation in
#         INFO/CSQ + the sample's GT/AD/DP genotype)
# OUTPUT: work/candidates.tsv    (ranked table read by 05_make_csv.py)
#   1. Reads each variant: gene, consequence, gnomAD, ClinVar, AlphaMissense,
#      SpliceAI , and from the genotype derives zygosity (het/hom) and VAF.
#   2. Separates PRIMARY (MVA panel genes) from SECONDARY (actionable ACMG SF
#      genes, only if ClinVar-pathogenic or LoF) and drops the rest.
#   3. Splits into tracks (recessive / dominant / mosaic) and pairs compound-het within the same gene.
#   4. Scores with a additive rubric: panel gene + phenotype match (HPO) + recessive completeness + molecular damage + rarity + mosaic. 
#
# USAGE:  python3 03_filter_rank.py            # uses ./config.yaml
import sys, gzip, os, yaml
from collections import defaultdict

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config.yaml"))
WORK = os.path.expanduser(cfg["paths"]["work"])
ANN  = os.path.join(WORK, "annotated.vcf.gz")
OUT  = os.path.join(WORK, "candidates.tsv")

MAX_AF_REC = cfg["freq"]["max_af_recessive"]
AM_MIN     = cfg["deleterious"]["alphamissense_min"]
HIGH_IMP   = set(cfg["deleterious"]["high_impact"])
SPLICEAI_MIN = cfg["deleterious"].get("spliceai_min", 0.5)
VAF_LO, VAF_HI = cfg["mosaic"]["vaf_min"], cfg["mosaic"]["vaf_max"]
MIN_DP     = cfg["mosaic"]["min_dp"]
SUBTYPE_HINT = str(cfg.get("subtype_hint", "none")).upper()
PHE_PER   = cfg.get("phenotype", {}).get("points_per_match", 0.5)
PHE_MAX   = cfg.get("phenotype", {}).get("max_points", 3)


RESCUE_CONS = ("splice_region_variant", "splice_donor_region_variant", "splice_donor_5th_base_variant", "splice_polypyrimidine_tract_variant",
               "5_prime_UTR_variant", "3_prime_UTR_variant", "upstream_gene_variant","synonymous_variant", "intron_variant")

def rescue_priority(consequence, af):
    c = consequence
    if any(s in c for s in ("splice_region", "splice_donor", "splice_acceptor",
                            "splice_polypyrimidine")):
        return 3                                   # splice edge
    if "synonymous" in c:
        return 2                                   # exonic: may break an exonic splice enhancer
    if ("UTR" in c) or ("upstream" in c):
        return 2                                   # promoter / regulatory region
    if ("intron" in c) and af == 0.0:
        return 1                                   # deep intronic but novel /absent in gnomAD
    return 0                                       # common intronic -> noise, not rescued

def load_gene_list(path):
    d = {}
    if not path or not os.path.exists(os.path.expanduser(path)):
        return d
    for line in open(os.path.expanduser(path)):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        d[p[0]] = p[1] if len(p) > 1 else "ext"
    return d

panel   = load_gene_list(cfg["paths"]["panel"])

for g in load_gene_list(cfg["paths"].get("phenocopy")):
    panel.setdefault(g, "phenocopy")
acmg_sf = set(load_gene_list(cfg["paths"].get("acmg_sf")))

# phenotype 
hpo_patient = set()
hpo_path = os.path.expanduser(cfg["paths"].get("hpo_terms", ""))
if hpo_path and os.path.exists(hpo_path):
    for line in open(hpo_path):
        t = line.split("#")[0].strip()
        if t.startswith("HP:"):
            hpo_patient.add(t)

# load annotated VCF and filter variants 
gene_hpo = defaultdict(set)
hg_path = os.path.expanduser(cfg["paths"].get("hpo_genes", ""))
if hg_path and os.path.exists(hg_path):
    for line in open(hg_path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3 and p[2].startswith("HP:"):
            gene_hpo[p[1]].add(p[2])   # p[1]=gene_symbol, p[2]=hpo_id

def phenotype_points(sym):
    if not hpo_patient or sym not in gene_hpo:
        return 0.0
    matches = len(hpo_patient & gene_hpo[sym])
    return min(matches * PHE_PER, PHE_MAX)

def af_float(x):
    try: return float(x)
    except (ValueError, TypeError): return 0.0

def zygosity_and_vaf(fmt, sample):
    d = dict(zip(fmt.split(":"), sample.split(":")))
    gt = d.get("GT", "./.").replace("|", "/")
    alleles = gt.split("/")
    is_hom = alleles.count("1") == 2
    is_het = ("1" in alleles) and not is_hom
    vaf = dp = None
    if d.get("AD", ".") not in (".", ""):
        ad = [int(x) for x in d["AD"].split(",") if x.isdigit()]
        if len(ad) >= 2 and sum(ad) > 0:
            vaf, dp = ad[1] / sum(ad), sum(ad)
    if dp is None and d.get("DP", ".").isdigit():
        dp = int(d["DP"])
    return is_het, is_hom, vaf, dp

# annotated vcf
csq_fields = []
variants = []
with gzip.open(ANN, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = line.split("Format:")[1].strip().strip('">').split("|")
            continue
        if line.startswith("#CHROM"):
            continue
        c = line.rstrip("\n").split("\t")
        chrom, pos, _id, ref, alt, _q, flt, info, fmt, sample = c[:10]
        csq = next((f[4:] for f in info.split(";") if f.startswith("CSQ=")), "")
        if not csq:
            continue
        entries = [dict(zip(csq_fields, e.split("|"))) for e in csq.split(",")]
        ann = next((e for e in entries
                    if e.get("SYMBOL") in panel or e.get("SYMBOL") in acmg_sf), entries[0])
        symbol = ann.get("SYMBOL", "")

        is_het, is_hom, vaf, dp = zygosity_and_vaf(fmt, sample)
        if not (is_het or is_hom):
            continue

        af = max(af_float(ann.get("gnomADg_AF")), af_float(ann.get("gnomADe_AF")))
        am = ann.get("am_pathogenicity", "")
        am = float(am) if am not in ("", ".") else None
        clnsig = ann.get("ClinVar_CLNSIG", "") or ""
        consequence = ann.get("Consequence", "")
        impact = ann.get("IMPACT", "")

        # SpliceAI: max of the DS_* delta scores
        ds = []
        for k, v in ann.items():
            if "SpliceAI" in k and "DS" in k and v not in ("", "."):
                try: ds.append(float(v))
                except ValueError: pass
        is_splice = bool(ds) and max(ds) >= SPLICEAI_MIN

        is_lof = impact == "HIGH" or any(x in consequence for x in HIGH_IMP)
        is_mis_path = am is not None and am >= AM_MIN
        is_clinvar_path = ("pathogenic" in clnsig.lower()) and ("conflicting" not in clnsig.lower())

        # strong : breaks the protein, pathogenic missense, ClinVar or splice.
        strict = is_lof or is_mis_path or is_clinvar_path or is_splice
        # moderate: possible second copy of a compound-het. Panel genes only.
        # Real VUS: missense with AlphaMissense >=0.34 (ambiguous or worse), splice_region, or in-frame indel. We exclude benign (ClinVar benign or AlphaMissense <0.34).
        is_missense = "missense_variant" in consequence
        clinvar_benign = "benign" in clnsig.lower()
        moderate = (not strict) and (not clinvar_benign) and (
            (is_missense and am is not None and am >= 0.34)
            or "splice_region_variant" in consequence
            or "inframe_insertion" in consequence or "inframe_deletion" in consequence
            or "protein_altering_variant" in consequence)

        if af > MAX_AF_REC:
            continue

        tier_sym = panel.get(symbol)             # 'core'/'ext'/'phenocopy' or None
        in_panel = tier_sym is not None
        is_mva   = tier_sym in ("core", "ext")
        in_acmg  = symbol in acmg_sf
        pheno_rel = phenotype_points(symbol) > 0  
        is_core   = tier_sym == "core"
        resc_prio = rescue_priority(consequence, af) if is_core else 0
        is_rescue = (is_core and resc_prio > 0 and not strict and not moderate and not clinvar_benign)
        #   - MVA (core/ext): strong AND moderate (moderate = candidate 2nd copy)
        #   - Phenocopies: strong only
        #   - any phenotypically relevant gene (HPO overlap): strong  <- genome-wide
        #   - any ClinVar pathogenic (whole genome): always reported
        #   - ACMG SF: LoF (actionable) -> secondary
        keep = ((is_mva and (strict or moderate))
                or (tier_sym == "phenocopy" and strict)
                or (pheno_rel and strict)
                or is_clinvar_path
                or (in_acmg and is_lof)
                or is_rescue)                      # <- cryptic 2nd copy in core genes
        if not keep:
            continue

        # FILTER: PASS or . is OK; non-PASS is OK if mosaic-like (VAF 0.05-0.35) or ClinVar-pathogenic
        mosaic_like = (vaf is not None and dp is not None
                       and VAF_LO <= vaf <= VAF_HI and dp >= MIN_DP)
        if flt not in ("PASS", ".") and not (mosaic_like or is_clinvar_path):
            continue

        # tier by priority: MVA/phenocopy panel > HPO-phenotype (genome-wide) > ACMG > ClinVar incidental
        if   in_panel:          tier = tier_sym          # core / ext / phenocopy
        elif pheno_rel:         tier = "pheno_hpo"        # phenotypically relevant gene (primary)
        elif in_acmg:           tier = "acmg_sf"          # actionable (secondary)
        else:                   tier = "clinvar_incidental"  # ClinVar pathogenic, no phenotype link (secondary)
        variants.append(dict(
            chrom=chrom, pos=int(pos), ref=ref, alt=alt, symbol=symbol, tier=tier,
            in_panel=in_panel, is_het=is_het, is_hom=is_hom, vaf=vaf, dp=dp, af=af,
            am=am, clnsig=clnsig, consequence=consequence, impact=impact, is_lof=is_lof,
            is_splice=is_splice, is_clinvar_path=is_clinvar_path, strict=strict, flt=flt,
            is_rescue=is_rescue, resc_prio=resc_prio,
        ))

#group by gene 
by_gene = defaultdict(list)
for v in variants:
    by_gene[v["symbol"]].append(v)

def deleterious_points(v):
    if v["is_clinvar_path"]: return 5
    if v["is_lof"] or v["is_splice"]: return 3
    if v["am"] is not None:
        if v["am"] >= AM_MIN: return 2          # pathogenic missense (AlphaMissense)
        if v["am"] >= 0.34:   return 1          # VUS missense (possible 2nd copy)
        return 0.5
    if not v["strict"]: return 1                # non-missense moderate (splice_region, inframe)
    return 0

def rarity_points(v):
    return 2 if v["af"] == 0.0 else 1

def gene_points(sym):
    tier = panel.get(sym, "acmg_sf" if sym in acmg_sf else "offpanel")
    base = {"core": 3, "ext": 1}.get(tier, 0)
    if SUBTYPE_HINT != "NONE" and sym == SUBTYPE_HINT:
        base += 2
    return base + phenotype_points(sym)    # + phenotype match (HPO)

rows = []
BIALLELIC_BONUS, SINGLE_HET_BONUS = 6, 1

def biallelic_bonus(a, b): # bonus for two strong alleles in the same gene (compound-het or homozygous)
    return BIALLELIC_BONUS if (a["strict"] and b["strict"]) else BIALLELIC_BONUS - 3

for sym, vs in by_gene.items():
    het  = [v for v in vs if v["is_het"]]
    het1 = [v for v in het if v["strict"]]                                 # strong
    hetm = [v for v in het if (not v["strict"]) and not v["is_rescue"]]    # moderate (real VUS)
    resc = [v for v in het if v["is_rescue"]]                              # rescue (core only)
    het2 = hetm + resc                                                     # any possible 2nd copy
    hom1 = [v for v in vs if v["is_hom"] and v["strict"]]

    # only in curated genes (MVA/phenocopy). Genome-wide, with no trio and no phasing
    # a compound-het is speculative and must not inflate the ranking
    allow_compound = panel.get(sym) in ("core", "ext", "phenocopy")
    pair = None
    if allow_compound and len(het1) >= 2:
        pair = tuple(sorted(het1, key=deleterious_points, reverse=True)[:2])
    elif allow_compound and len(het1) == 1 and het2:
        # 2nd copy: first a real moderate VUS; if none, a rescue variant
        second = max(hetm, key=deleterious_points) if hetm else max(resc, key=lambda v: v["resc_prio"])
        pair = (het1[0], second)
    elif allow_compound and len(het1) == 0 and len(hetm) >= 2:
        pair = tuple(sorted(hetm, key=deleterious_points, reverse=True)[:2])
    if pair:
        a, b = pair
        if a["is_rescue"] or b["is_rescue"]:
            rev = "2nd_allele_rescue_splicing(review)"
        elif not (a["strict"] and b["strict"]):
            rev = "2nd_allele_VUS(review)"           # unconfirmed 2nd allele
        else:
            rev = ""
        score = (gene_points(sym) + biallelic_bonus(a, b)
                 + deleterious_points(a) + deleterious_points(b)
                 + rarity_points(a) + rarity_points(b))
        rows.append(dict(kind="compound_het", score=score, v1=a, v2=b, symbol=sym, review=rev))

    #strong homozygous
    for v in hom1:
        score = gene_points(sym) + BIALLELIC_BONUS + deleterious_points(v) + rarity_points(v)
        rows.append(dict(kind="homozygous", score=score, v1=v, v2=None, symbol=sym, review=""))

    #lone strong allele 
    for v in het1:
        score = gene_points(sym) + SINGLE_HET_BONUS + deleterious_points(v) + rarity_points(v)
        rows.append(dict(kind="single_het", score=score, v1=v, v2=None, symbol=sym, review=""))

    #best lone moderate VUS (low-confidence backup)
    used = {id(x) for x in (pair or ())}
    spare = [v for v in hetm if id(v) not in used]
    if spare:
        v = max(spare, key=deleterious_points)
        score = gene_points(sym) + 0.5 + deleterious_points(v) + rarity_points(v)
        rows.append(dict(kind="het_vus", score=score, v1=v, v2=None, symbol=sym, review=""))

    #up to 2 intronic/regulatory variants per gene, with their own low ranking tier, so they do not push out real candidates
    spare_resc = sorted([v for v in resc if id(v) not in used],
                        key=lambda v: (v["resc_prio"], rarity_points(v)), reverse=True)[:2]
    for v in spare_resc:
        score = gene_points(sym) + deleterious_points(v) + rarity_points(v)
        rows.append(dict(kind="core_rescue", score=score, v1=v, v2=None, symbol=sym,
                         rtier="core_rescue", review="review_cryptic_splicing"))

for r in rows:
    v = r["v1"]
    r["mosaic"] = (v["vaf"] is not None and v["dp"] is not None
                   and VAF_LO <= v["vaf"] <= VAF_HI and v["dp"] >= MIN_DP)
    if r["mosaic"]:
        r["score"] += 1
    if v["flt"] == "PASS":
        r["score"] += 0.5

# ranking by evidence:
#   1) primary before secondary
#   2) gene prior: known MVA gene > phenocopy > HPO-phenotype > (secondary)
#   3) ClinVar pathogenic evidence
#   4) score (tie-break)
PRIMARY_TIERS = {"core", "ext", "phenocopy", "pheno_hpo", "core_rescue"}
TIER_RANK = {"core": 5, "ext": 4, "phenocopy": 3, "pheno_hpo": 2, "core_rescue": 1, "acmg_sf": 1, "clinvar_incidental": 1}
def etier(r):
    return r.get("rtier") or r["v1"]["tier"]
def rank_key(r):
    v = r["v1"]; t = etier(r)
    return (t in PRIMARY_TIERS, TIER_RANK.get(t, 0),
            bool(v["is_clinvar_path"]), r["score"])
rows.sort(key=rank_key, reverse=True)

def fmt_v(v): return f"{v['chrom']}:{v['pos']}:{v['ref']}>{v['alt']}"

with open(OUT, "w") as out:
    out.write("rank\tkind\tscore\tsymbol\ttier\tin_panel\tmosaic\t"
              "chrom_1\tpos_1\tref_1\talt_1\tchrom_2\tpos_2\tref_2\talt_2\t"
              "vaf_1\taf_1\tam_1\tclnsig_1\tconsequence_1\treview\n")
    for i, r in enumerate(rows, 1):
        a, b = r["v1"], r["v2"]
        c2 = [b["chrom"], b["pos"], b["ref"], b["alt"]] if b else ["", "", "", ""]
        out.write("\t".join(str(x) for x in [ i, r["kind"], round(r["score"], 1), r["symbol"], etier(r), a["in_panel"], r["mosaic"],
            a["chrom"], a["pos"], a["ref"], a["alt"], *c2,
            (round(a["vaf"], 3) if a["vaf"] is not None else ""), a["af"],
            (a["am"] if a["am"] is not None else ""), a["clnsig"], a["consequence"],
            r.get("review", ""),]) + "\n")

LABEL = {"pheno_hpo": " [HPO-phenotype]", "phenocopy": " [phenocopy]", "core_rescue": " [splicing-rescue-review]", "acmg_sf": " [ACMG secondary]", "clinvar_incidental": " [ClinVar incidental]"}
print(f"[rank] candidates: {len(rows)}  ->  {OUT}")
for i, r in enumerate(rows[:10], 1):
    a = r["v1"]
    tags = ("+mosaic" if r["mosaic"] else "") + LABEL.get(etier(r), "")
    print(f"  #{i:<2} {r['symbol']:<9} {r['kind']:<13} score={r['score']:<5} {fmt_v(a)} {tags}")