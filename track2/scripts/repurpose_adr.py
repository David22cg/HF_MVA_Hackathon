# aneuploidy-directed, checkpoint-safe network repurposing 
# Classic network proximity (Guney 2016) ranks a drug by how close its targets sit to the disease gene;
# Usage: python3 repurpose_adr.py [config_2.yaml]
# Needs: networkx, numpy, matplotlib, pyyaml, and the STRING files.
import sys, os, gzip, csv, yaml, random
from collections import deque, defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config_2.yaml"))
P = cfg["paths"]
MIN_SCORE = int(cfg.get("string_min_score", 700))
NRAND     = int(cfg.get("n_random", 1000))      # degree-matched null samples per drug
W = cfg.get("weights", {"stress": 1.0, "dependency": 1.5, "kinetic": 1.25, "approved": 1.0, "danger": 1.0})
W.setdefault("kinetic", 1.25)
# mitotic kinases that could slow mitosis at low dose but are vetoed by default;
# flagged only as a dose-dependent hypothesis, never ranked as safe.
MITOTIC_SLOWERS = {"AURKB", "AURKA", "PLK1", "CDK1", "TTK"}
SIG_Z = float(cfg.get("significant_z", -1.64))  # z <= this is "significantly proximal"
OUT = os.path.expanduser(P.get("out", "./out")); os.makedirs(OUT, exist_ok=True)
random.seed(7); np.random.seed(7)

# refuse the demo drug table unless allowed, so the search uses the real database
_dt = os.path.basename(os.path.expanduser(P["drug_targets"])).lower()
if "demo" in _dt and not cfg.get("allow_demo", False):
    sys.exit("[adr] refusing the demo drug table. Build the real one:\n"
             "  python3 drug_targets.py interactions.tsv\n"
             "  then set paths.drug_targets: ./data/drug_targets_real.tsv (or allow_demo: true).")


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def load_set(path):
    # one gene symbol per line, '#' comments ignored
    s = set()
    for line in open(os.path.expanduser(path)):
        line = line.split("#")[0].strip()
        if line:
            s.add(line.split()[0].upper())
    return s


def load_string_graph():
    # build a gene-symbol graph from STRING edges with combined score >= MIN_SCORE
    links = os.path.expanduser(P["string_links"]); info = os.path.expanduser(P["string_info"])
    if not (os.path.exists(links) and os.path.exists(info)):
        sys.exit(f"ERROR: STRING files not found ({links} / {info}); same files as Track 1.")
    import networkx as nx
    prot2sym = {}
    with opener(info) as fh:
        next(fh, None)
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) >= 2:
                prot2sym[a[0]] = a[1].upper()
    G = nx.Graph()
    with opener(links) as fh:
        next(fh, None)
        for line in fh:
            a = line.split()
            if len(a) < 3:
                continue
            try:
                sc = int(a[2])
            except ValueError:
                continue
            if sc < MIN_SCORE:
                continue
            s1, s2 = prot2sym.get(a[0]), prot2sym.get(a[1])
            if s1 and s2 and s1 != s2:
                G.add_edge(s1, s2)
    return G


def dist_to_module(G, module):
    # multi-source BFS: shortest distance from every node to the nearest module gene
    seeds = [g for g in module if g in G]
    dist = {g: 0 for g in seeds}
    dq = deque(seeds)
    while dq:
        u = dq.popleft()
        for v in G.neighbors(u):
            if v not in dist:
                dist[v] = dist[u] + 1
                dq.append(v)
    return dist, seeds


def degree_bins(G, nbins=10):
    # group nodes into degree bins so the null draws degree-matched random targets
    nodes = list(G.nodes()); deg = np.array([G.degree(n) for n in nodes])
    order = np.argsort(deg)
    bin_of = {}; binmembers = defaultdict(list)
    for rank, idx in enumerate(order):
        b = int(rank * nbins / len(nodes))
        bin_of[nodes[idx]] = b; binmembers[b].append(nodes[idx])
    return bin_of, binmembers


def proximity_z(targets_in_G, dist, bin_of, binmembers, big):
    # Guney "closest" distance with a degree-preserving z-score and an empirical
    # one-sided p-value. Returns (d_c, z, p); small z/p means closer than random.
    d_obs = np.mean([dist.get(t, big) for t in targets_in_G])
    null = np.empty(NRAND)
    for i in range(NRAND):
        samp = [random.choice(binmembers[bin_of[t]]) for t in targets_in_G]
        null[i] = np.mean([dist.get(s, big) for s in samp])
    mu, sd = null.mean(), null.std()
    z = (d_obs - mu) / sd if sd > 0 else 0.0
    p = (int(np.sum(null <= d_obs)) + 1) / (NRAND + 1)   # +1 pseudocount so p != 0
    return float(d_obs), float(z), float(p)


def benjamini_hochberg(pvals):
    # BH false-discovery-rate q-values for a list of p-values 
    idx = [i for i, p in enumerate(pvals) if p == p]
    m = len(idx)
    q = [float("nan")] * len(pvals)
    if m == 0:
        return q
    order = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        prev = q[i] = min(prev, pvals[i] * m / k)
    return q


def module_lcc(G, module, bin_of, binmembers, nrand):
    # disease-module test (Menche 2015): is the largest connected component of the
    # module bigger than for degree-matched random gene sets of the same size?
    import networkx as nx
    seeds = [g for g in module if g in G]
    if len(seeds) < 2:
        return None
    size = lambda nodes: max((len(c) for c in nx.connected_components(G.subgraph(nodes))), default=0)
    def rand_set():
        out = set()
        for g in seeds:
            for _ in range(10):
                c = random.choice(binmembers[bin_of[g]])
                if c not in out:
                    out.add(c); break
            else:
                out.add(c)
        return out
    obs = size(seeds)
    null = np.array([size(rand_set()) for _ in range(nrand)])
    sd = null.std()
    z = (obs - null.mean()) / sd if sd > 0 else 0.0
    p = (int(np.sum(null >= obs)) + 1) / (nrand + 1)
    return dict(n=len(seeds), lcc=obs, rand=round(float(null.mean()), 1), z=round(float(z), 2), p=round(float(p), 4))


G = load_string_graph()
BIG = G.number_of_nodes()                      # distance for targets not reachable
stress = load_set(P["stress_module"]); checkpoint = load_set(P["checkpoint_module"])
deps   = load_set(P["dependencies"])
kinetic = load_set(P["kinetic_module"]) if P.get("kinetic_module") else set()

# disease anchor from Track 1-> merge the predicted genes into the veto set
seed_genes = load_set(P["seed_genes"]) if P.get("seed_genes") and os.path.exists(os.path.expanduser(P["seed_genes"])) else set()
if seed_genes:
    checkpoint |= seed_genes
    print(f"[adr] Track-1 seed(s): {sorted(seed_genes)} added to the veto anchor")
else:
    print("[adr] no Track-1 seeds (run from_track1.py); using the static checkpoint list only")

CHK = cfg.get("checkpoint", {})
if CHK.get("from_seed") and seed_genes:
    radius = int(CHK.get("radius", 1))
    derived = set(g for g in seed_genes if g in G)
    frontier = set(derived)
    for _ in range(max(1, radius)):
        nxt = set()
        for u in frontier:
            nxt |= set(G.neighbors(u))
        frontier = nxt - derived
        derived |= nxt
    derived |= seed_genes
    before = len(checkpoint)
    checkpoint = derived if CHK.get("replace_static", True) else (checkpoint | derived)
    print(f"[adr] veto set derived from seed+STRING (radius {radius}): "
          f"{len(checkpoint)} genes (was {before} static)")

dist_stress, seed_s = dist_to_module(G, stress)
dist_check,  seed_c = dist_to_module(G, checkpoint)
bin_of, binmembers = degree_bins(G)
HUB_CUT = float(np.percentile([d for _, d in G.degree()], 99))   # top 1% degree = hub
print(f"[adr] STRING graph: {G.number_of_nodes()} genes, {G.number_of_edges()} edges; "
      f"stress seeds={len(seed_s)}, checkpoint seeds={len(seed_c)}")

# disease-module connectivity check (writes module_validation.txt)
with open(os.path.join(OUT, "module_validation.txt"), "w") as mf:
    mf.write("disease-module connectivity: largest connected component vs degree-matched random\n")
    for name, mod in (("stress", stress), ("checkpoint", checkpoint)):
        v = module_lcc(G, mod, bin_of, binmembers, min(NRAND, 500))
        if v:
            ln = f"{name:<11} genes={v['n']} LCC={v['lcc']} random={v['rand']} z={v['z']} p={v['p']}"
            mf.write(ln + "\n"); print("[adr] module", ln)

rows = []
for r in csv.DictReader(open(os.path.expanduser(P["drug_targets"])), delimiter="\t"):
    drug = r["drug"]; status = r.get("status", ""); ped = r.get("pediatric", "")
    tg = [t.strip().upper() for t in r["targets"].split(";") if t.strip()]
    in_G = [t for t in tg if t in G]
    on_checkpoint = [t for t in tg if t in checkpoint]
    on_deps       = [t for t in tg if t in deps]
    on_kinetic    = [t for t in tg if t in kinetic]
    on_hub        = [t for t in in_G if G.degree(t) >= HUB_CUT]   # high-degree, likely toxic to inhibit
    if not in_G:
        rows.append(dict(drug=drug, status=status, pediatric=ped, targets=";".join(tg),
                         z_stress=float("nan"), z_check=float("nan"), closeness_stress=float("nan"),
                         dependency=bool(on_deps), kinetic=bool(on_kinetic), hub=bool(on_hub),
                         approved=status.startswith("approved"),
                         veto=bool(on_checkpoint), score=float("-inf"),
                         p_stress=float("nan"), q_stress=float("nan"),
                         recommendation="UNSCORABLE (no target in network)",
                         reason="targets absent from STRING")); continue
    _, z_stress, p_stress = proximity_z(in_G, dist_stress, bin_of, binmembers, BIG)
    _, z_check,  _        = proximity_z(in_G, dist_check,  bin_of, binmembers, BIG)
    closeness_stress = -z_stress               # higher = closer to the stress module
    dependency = bool(on_deps)
    kinetic_rescue = bool(on_kinetic)
    approved = status.startswith("approved")
    approved_score = (1.0 if approved else 0.0) + (0.5 if ped == "yes" else 0.0)
    veto = bool(on_checkpoint)                 # drug acts on the mitotic-fidelity machinery

    # a drug scores via any of three routes (stress proximity, dependency, kinetic);
    # distance from the stress module is neutral, not a penalty. Safety is the veto.
    stress_reward = max(0.0, closeness_stress)
    score = (W["stress"] * stress_reward
             + W["dependency"] * (1.0 if dependency else 0.0)
             + W["kinetic"]   * (1.0 if kinetic_rescue else 0.0)
             + W["approved"]  * approved_score)
    if veto:                                   # safety first
        on_slow = [t for t in on_checkpoint if t in MITOTIC_SLOWERS]
        rec = "SAFETY-VETO"
        reason = (f"acts on the mitotic-fidelity module ({','.join(on_checkpoint)}); "
                  "inhibiting it worsens missegregation in a SAC-deficient patient")
        if on_slow:
            reason += ("  [dose-dependent kinetic-rescue hypothesis at sub-lethal dose; "
                       "high-risk, needs validation, not ranked as safe]")
        score = -100 + score                   # sinks below everything, order kept among vetoed
    elif dependency:
        rec = "PRIORITIZE (aneuploidy dependency)"
        reason = (f"targets a synthetic-lethal aneuploidy dependency ({','.join(on_deps)}); "
                  "kills aneuploid cells without breaking the checkpoint")
    elif kinetic_rescue:
        rec = "PRIORITIZE (kinetic rescue: slow the cell cycle)"
        reason = (f"slows the cell cycle via {','.join(on_kinetic)} (G1/S or mTOR) without "
                  "inhibiting the SAC; gives error-correction more time (Hyp.1)")
    elif z_stress <= SIG_Z and approved:
        rec = "PRIORITIZE (approved, stress-proximal)"
        reason = f"significantly proximal to the stress module (z={z_stress:.2f}) and approved"
    elif z_stress <= SIG_Z:
        rec = "CONSIDER (stress-proximal, not approved)"
        reason = f"significantly proximal to the stress module (z={z_stress:.2f}) but investigational"
    elif approved:
        rec = "CONSIDER (approved, weak proximity)"
        reason = f"approved but not significantly stress-proximal (z={z_stress:.2f})"
    else:
        rec = "LOW"
        reason = f"not stress-proximal (z={z_stress:.2f})"
    if on_hub:                                 # centrality is not therapeutic importance
        reason += f"  [hub target ({','.join(on_hub)}), may be toxic to inhibit]"
    rows.append(dict(drug=drug, status=status, pediatric=ped, targets=";".join(tg),
                     z_stress=round(z_stress, 2), z_check=round(z_check, 2),
                     closeness_stress=round(closeness_stress, 2),
                     p_stress=round(p_stress, 4), q_stress=float("nan"),
                     dependency=dependency, kinetic=kinetic_rescue, hub=bool(on_hub),
                     approved=approved, veto=veto,
                     score=round(score, 2), recommendation=rec, reason=reason))

# multiple-testing correction over scorable, non-vetoed drugs (vetoed are safety calls)
elig = [r for r in rows if not r["veto"] and r["p_stress"] == r["p_stress"]]
for r, q in zip(elig, benjamini_hochberg([r["p_stress"] for r in elig])):
    r["q_stress"] = round(q, 4)
n_sig_p = sum(1 for r in elig if r["p_stress"] <= 0.05)
n_sig_q = sum(1 for r in elig if r["q_stress"] == r["q_stress"] and r["q_stress"] <= 0.10)
print(f"[adr] stress proximity: {n_sig_p}/{len(elig)} at p<=0.05, {n_sig_q} survive FDR q<=0.10")

# validation: literature aneuploidy-selective compounds should recover (non-vetoed, FDR-significant)
known_pos = load_set(P["known_positives"]) if P.get("known_positives") and \
    os.path.exists(os.path.expanduser(P["known_positives"])) else set()
if known_pos:
    present = [r for r in rows if r["drug"].upper() in known_pos]
    recovered = [r for r in present if (not r["veto"])
                 and r["q_stress"] == r["q_stress"] and r["q_stress"] <= 0.10]
    lines = [f"positive controls present: {len(present)}/{len(known_pos)}",
             f"recovered (non-vetoed and FDR q<=0.10): {len(recovered)}"]
    for r in present:
        tag = "recovered" if r in recovered else ("vetoed" if r["veto"] else "not significant")
        lines.append(f"  {r['drug']:<18} p={r['p_stress']} q={r['q_stress']} -> {tag}")
    with open(os.path.join(OUT, "validation.txt"), "w") as vf:
        vf.write("positive-control recovery (Tang 2011 aneuploidy-selective compounds)\n")
        vf.write("\n".join(lines) + "\n")
    print("[adr] validation:", lines[0] + "; " + lines[1])

# rank non-vetoed by score, vetoed last
rows.sort(key=lambda r: (not r["veto"], r["score"]), reverse=True)

cols = ["rank", "drug", "recommendation", "score", "z_stress", "p_stress", "q_stress",
        "z_check", "dependency", "kinetic", "hub", "veto", "approved", "pediatric", "status",
        "targets", "reason"]
dest = os.path.join(OUT, "adr_ranking.tsv")
with open(dest, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for i, r in enumerate(rows, 1):
        fh.write("\t".join([str(i)] + [str(r[c]) for c in cols[1:]]) + "\n")

print(f"[adr] ranked {len(rows)} drugs -> {dest}")
for i, r in enumerate(rows[:15], 1):      # show only the top 15 in the console
    print(f"  #{i:<2} {r['drug']:<16} {r['recommendation']:<40} score={r['score']}")

# map: faint cloud of all scored drugs, labels only on the top hits 
sc = [r for r in rows if isinstance(r["closeness_stress"], float) and r["closeness_stress"] == r["closeness_stress"]]
if sc:
    top = [r for r in sc if not r["veto"]][:15]          # label only these
    fig, ax = plt.subplots(figsize=(8.5, 7))
    for r in sc:
        x, y = r["closeness_stress"], -r["z_check"]
        if r["veto"]:              col, mk, a, s = "#c0392b", "x", 0.5, 45
        elif r["dependency"]:      col, mk, a, s = "#15803d", "*", 1.0, 150
        elif r.get("kinetic"):     col, mk, a, s = "#2563eb", "D", 1.0, 70
        elif r["recommendation"].startswith("PRIORITIZE"): col, mk, a, s = "#15803d", "o", 0.9, 70
        else:                      col, mk, a, s = "#9aa0a6", "o", 0.25, 25
        ax.scatter(x, y, c=col, marker=mk, s=s, alpha=a, edgecolors="none", zorder=2)
    for r in top:                                        # labels on top hits only
        ax.annotate(r["drug"], (r["closeness_stress"], -r["z_check"]),
                    fontsize=7.5, xytext=(4, 3), textcoords="offset points")
    ax.axhline(1.64, color="#c0392b", lw=0.8, ls="--", alpha=0.5)
    ax.axvline(1.64, color="#15803d", lw=0.8, ls="--", alpha=0.5)
    ax.set_xlabel("closeness to the stress module (-z)  ->  exploit the consequence")
    ax.set_ylabel("closeness to the SAC module (-z)  ->  danger (vetoed)")
    ax.set_title("ADR-Net: ideal is lower-right (stress-close, SAC-far). "
                 "green = prioritize, blue = kinetic, red x = veto", fontsize=9.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "adr_map.png"), dpi=130)
    print(f"[adr] map -> {OUT}/adr_map.png")