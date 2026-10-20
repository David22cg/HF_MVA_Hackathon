#   (A) mosaic aneuploidy score:
#       Estimates the mosaicism fraction per chromosome from the B-allele frequency  deviation of heterozygous SNPs
#
#   (B) cryptic 2nd copy detector:
#       In the candidate genes, it looks for loss of heterozygosity, a depth drop or allelic imbalance
import sys, os, csv, subprocess, yaml
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config.yaml"))
WORK = os.path.expanduser(cfg["paths"]["work"])
OUT  = os.path.expanduser(cfg["paths"]["out"]); os.makedirs(OUT, exist_ok=True)
VCF  = os.path.join(WORK, "all.norm.vcf.gz")
CAND = os.path.join(WORK, "candidates.tsv")
MIN_DP = cfg["mosaic"]["min_dp"]

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
CHR_LEN = {"1":248956422,"2":242193529,"3":198295559,"4":190214555,"5":181538259,
    "6":170805979,"7":159345973,"8":145138636,"9":138394717,"10":133797422,
    "11":135086622,"12":133275309,"13":114364328,"14":107043718,"15":101991189,
    "16":90338345,"17":83257441,"18":80373285,"19":58617616,"20":64444167,
    "21":46709983,"22":50818468,"X":156040895,"Y":57227415}

#extract heterozygous SNPs (pos, BAF, depth) per chromosome
def get_het_table():
    cached = os.path.join(WORK, "het_snps.tsv")
    if os.path.exists(cached):
        src = open(cached)
    else:
        p1 = subprocess.Popen(["bcftools", "view", "-v", "snps", "-m2", "-M2", VCF],
                              stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT\t%AD]\n"],
                              stdin=p1.stdout, stdout=subprocess.PIPE, text=True)
        p1.stdout.close()
        src = p2.stdout
    baf = {c: [] for c in CHROMS}; dep = {c: [] for c in CHROMS}; pos = {c: [] for c in CHROMS}
    for line in src:
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        c, p, gt, ad = f[0], f[1], f[2].replace("|", "/"), f[3]
        if c not in baf or gt not in ("0/1", "1/0"):
            continue
        try:
            r, a = (int(x) for x in ad.split(",")[:2])
        except ValueError:
            continue
        tot = r + a
        if tot < MIN_DP:
            continue
        baf[c].append(a / tot); dep[c].append(tot); pos[c].append(int(p))
    return ({c: np.array(baf[c]) for c in CHROMS},
            {c: np.array(dep[c]) for c in CHROMS},
            {c: np.array(pos[c]) for c in CHROMS})

baf, dep, pos = get_het_table()
all_dep = np.concatenate([dep[c] for c in CHROMS if dep[c].size]) if any(dep[c].size for c in CHROMS) else np.array([1.0])
genome_med_depth = float(np.median(all_dep))
total_het = int(sum(baf[c].size for c in CHROMS))
genome_het_density = total_het / (sum(CHR_LEN.values()) / 1e6)   # het per Mb

#a)mean |BAF-0.5| per chromosome; baseline = median of the autosomes (almost all normal)
devs = {c: (float(np.mean(np.abs(baf[c] - 0.5))) if baf[c].size else np.nan) for c in CHROMS}
auto_devs = [devs[c] for c in [str(i) for i in range(1, 23)] if not np.isnan(devs[c])]
baseline = float(np.median(auto_devs)) if auto_devs else 0.0

def mosaic_fraction(d):
    exc = max(0.0, d - baseline)
    if exc <= 0:
        return 0.0
    return float(min(1.0, 4 * exc / max(1e-6, 1 - 2 * exc)))

burden = []
for c in CHROMS:
    if baf[c].size < 50:
        burden.append((c, baf[c].size, np.nan, np.nan, np.nan, "")); continue
    d = devs[c]
    frac = mosaic_fraction(d)
    dr = float(np.median(dep[c]) / genome_med_depth)
    flag = "MOSAIC_ANEUPLOIDY" if (frac >= 0.10 or dr > 1.15 or dr < 0.85) else ""
    burden.append((c, baf[c].size, round(d, 3), round(frac, 3), round(dr, 3), flag))

with open(os.path.join(OUT, "aneuploidy_burden.tsv"), "w") as fh:
    fh.write("chrom\tn_het\tbaf_dev\test_mosaic_fraction\tdepth_ratio\tflag\n")
    fh.write(f"# baseline baf_dev (diploid noise) = {round(baseline,3)}; genome median depth = {round(genome_med_depth,1)}\n")
    for row in burden:
        fh.write("\t".join(str(x) for x in row) + "\n")

flagged = [c for c, n, d, fr, dr, fl in burden if fl]
max_frac = max([fr for c, n, d, fr, dr, fl in burden if isinstance(fr, float) and not np.isnan(fr)] or [0])
print(f"[aneuploidy] chromosomes with mosaic signal: {flagged or 'none'}")
print(f"[aneuploidy] max estimated mosaic fraction: {max_frac:.2f}  "
      f"(detection floor ~0.10-0.15; low/variegated mosaicism may be invisible in bulk)")

# b) 2nd-copy detector via LOH / allelic imbalance
W = 50000                   
MIN_EXPECTED = 15            
chrom_density = {c: (pos[c].size / (CHR_LEN[c] / 1e6)) if pos[c].size else 0.0 for c in CHROMS}  # het/Mb

def fnum(x):                  # format NaN -> "" for the TSV
    return "" if isinstance(x, float) and np.isnan(x) else x

loh_rows = []
if os.path.exists(CAND):
    seen = set()
    for r in csv.DictReader(open(CAND), delimiter="\t"):
        if r.get("tier") not in ("core", "ext", "phenocopy", "pheno_hpo"):   # primary candidates only
            continue
        sym = r["symbol"]; c = str(r["chrom_1"]).replace("chr", "")
        if (sym, c) in seen:
            continue
        seen.add((sym, c))
        try:
            p = int(r["pos_1"])
        except (ValueError, KeyError):
            continue
        if c not in pos:
            loh_rows.append((sym, c, 0, "", "", "", "", "chromosome out of range")); continue
        expected = chrom_density[c] * (2 * W / 1e6)             # expected het in the window
        m = (pos[c] >= p - W) & (pos[c] <= p + W)
        n = int(m.sum())
        densr = (n / expected) if expected > 0 else float("nan")  # <1 = het depletion = LOH
        if n >= 4:
            dr  = float(np.median(dep[c][m]) / genome_med_depth)
            imb = float(abs(np.mean(baf[c][m]) - 0.5))
        else:
            dr = imb = float("nan")
        strong_depletion = (expected >= MIN_EXPECTED and densr < 0.30)
        depth_drop = (n >= 4 and not np.isnan(dr) and dr < 0.75)
        imbalance  = (n >= 4 and not np.isnan(imb) and imb > 0.15)
        # excess het + excess depth = repetitive region / segmental duplication / copy gain 
        repetitive = (n >= 4 and not np.isnan(dr) and dr > 2.0 and not np.isnan(densr) and densr > 3.0)
        if expected < MIN_EXPECTED:
            flag = "too few expected het in the window (inconclusive)"
        elif depth_drop or (strong_depletion and n <= 3):
            flag = "SUSPECTED_BIALLELIC: deleted 2nd allele (LOH + depth drop)"
        elif strong_depletion:
            flag = "SUSPECTED_BIALLELIC: 2nd allele by LOH (deletion or UPD; depth inconclusive)"
        elif repetitive:
            flag = "repetitive/gain region (low confidence; variants here unreliable)"
        elif imbalance:
            flag = "SUSPECTED_BIALLELIC: allelic imbalance (possible 2nd allele)"
        else:
            flag = "no 2nd-allele signal (normal heterozygous region)"
        loh_rows.append((sym, c, n, round(expected, 1),
                         round(dr, 2) if not np.isnan(dr) else float("nan"),
                         round(densr, 2) if not np.isnan(densr) else float("nan"),
                         round(imb, 3) if not np.isnan(imb) else float("nan"), flag))

with open(os.path.join(OUT, "loh_scan.tsv"), "w") as fh:
    fh.write("gene\tchrom\tn_het_window\tn_het_expected\tdepth_ratio\trel_density\tallelic_imbalance\tflag\n")
    for row in loh_rows:
        fh.write("\t".join(str(fnum(x)) for x in row) + "\n")

hits = [r[0] for r in loh_rows if "SUSPECTED_BIALLELIC" in str(r[-1])]
print(f"[2nd-copy] genes with suspected biallelic (2nd allele by LOH/deletion): {hits or 'none'}  -> out/loh_scan.tsv")

# figures
fig, axes = plt.subplots(4, 6, figsize=(18, 10), sharey=True)
for ax, c in zip(axes.flat, CHROMS):
    if baf[c].size:
        x = pos[c]; y = baf[c]
        if x.size > 5000:
            idx = np.random.choice(x.size, 5000, replace=False); x, y = x[idx], y[idx]
        ax.scatter(x / 1e6, y, s=1, alpha=0.2)
    ax.axhline(0.5, color="grey", lw=0.6)
    ax.set_title(f"chr{c}", fontsize=9); ax.set_ylim(0, 1); ax.set_xticks([])
for ax in axes.flat[len(CHROMS):]:
    ax.axis("off")
fig.suptitle("B-allele frequency per chromosome (het SNPs). Band split away from 0.5 = mosaic aneuploidy")
fig.supylabel("BAF = alt/(ref+alt)")
fig.tight_layout(rect=[0.01, 0, 1, 0.97])
fig.savefig(os.path.join(OUT, "baf_map.png"), dpi=120)

fig2, ax2 = plt.subplots(figsize=(12, 4))
ratios = [dr if isinstance(dr, float) and not np.isnan(dr) else 0 for _, _, _, _, dr, _ in burden]
colors = ["#c0392b" if r > 1.15 else "#2980b9" if (0 < r < 0.85) else "#7f8c8d" for r in ratios]
ax2.bar(range(len(CHROMS)), ratios, color=colors)
ax2.axhline(1.0, color="black", lw=0.8)
ax2.set_xticks(range(len(CHROMS))); ax2.set_xticklabels([f"chr{c}" for c in CHROMS], rotation=90, fontsize=8)
ax2.set_ylabel("median depth / genome"); ax2.set_title("Depth karyotype (red=gain, blue=loss)")
fig2.tight_layout(); fig2.savefig(os.path.join(OUT, "depth_karyotype.png"), dpi=120)

print(f"[maps] -> {OUT}/baf_map.png, depth_karyotype.png, aneuploidy_burden.tsv, loh_scan.tsv")