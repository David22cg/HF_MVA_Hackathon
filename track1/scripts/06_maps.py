#   (1) interaction network-> string. The subgraph is induced on the
#       candidate genes + their highest-confidence neighbors, so it
#       reflects actual interactions and scales to any candidate set
#   (2) chromosome ideogram with the position of each candidate variant marked
#
# INPUT : STRING (optional, local): paths.string_links + paths.string_info download once (human, taxon 9606) from https://string-db.org/
#   (requires: pip install networkx matplotlib)
# =====================================================================
import sys, os, csv, gzip, yaml
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config.yaml"))
WORK = os.path.expanduser(cfg["paths"]["work"])
OUT  = os.path.expanduser(cfg["paths"]["out"]); os.makedirs(OUT, exist_ok=True)
CAND = os.path.join(WORK, "candidates.tsv")

STRING_LINKS = os.path.expanduser(cfg["paths"].get("string_links", "") or "")
STRING_INFO  = os.path.expanduser(cfg["paths"].get("string_info", "") or "")
STRING_MIN   = int(cfg.get("network", {}).get("string_min_score", 700))   # high-confidence edges
MAX_NEIGH    = int(cfg.get("network", {}).get("max_neighbors", 8))         # per candidate, caps the hairball


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


#candidate genes-loci
cand_genes, cand_loci = set(), []
if os.path.exists(CAND):
    for r in csv.DictReader(open(CAND), delimiter="\t"):
        cand_genes.add(r["symbol"].upper())
        try:
            cand_loci.append((r["symbol"], r["chrom_1"], int(r["pos_1"]), r["tier"]))
        except (ValueError, KeyError):
            pass

#intercation network fallback if string files are not available
CURATED_EDGES = [
    ("MAD1L1", "MAD2L1"), ("MAD2L1", "CDC20"), ("TRIP13", "MAD2L1"),
    ("MAD2L1BP", "MAD2L1"), ("MAD2L1BP", "TRIP13"), ("BUB3", "BUB1B"),
    ("BUB3", "BUB1"), ("BUB1", "BUB1B"), ("BUB1B", "CDC20"), ("TTK", "BUB1"),
    ("TTK", "MAD1L1"), ("CENPE", "BUB1B"), ("CEP57", "CENPE"), ("CDC20", "APC/C"),
]


def build_string_network():
    if not (STRING_LINKS and STRING_INFO and os.path.exists(STRING_LINKS) and os.path.exists(STRING_INFO)):
        return None, None
    prot2sym = {}
    with opener(STRING_INFO) as fh:
        next(fh, None)                                   # header
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                prot2sym[p[0]] = p[1].upper()
    #keep high-confidence edges incident to a candidate gene
    inc = defaultdict(list)                               # candidate -> [(neighbor, score)]
    cand_cand = set()                                     # edges between two candidates
    with opener(STRING_LINKS) as fh:
        next(fh, None)                                   # header: protein1 protein2 combined_score
        for line in fh:
            a = line.split()
            if len(a) < 3:
                continue
            try:
                score = int(a[2])
            except ValueError:
                continue
            if score < STRING_MIN:
                continue
            s1, s2 = prot2sym.get(a[0]), prot2sym.get(a[1])
            if not s1 or not s2 or s1 == s2:
                continue
            c1, c2 = s1 in cand_genes, s2 in cand_genes
            if c1 and c2:
                cand_cand.add((min(s1, s2), max(s1, s2)))
            elif c1:
                inc[s1].append((s2, score))
            elif c2:
                inc[s2].append((s1, score))
    #cap neighbors per candidate to MAX_NEIGH, and return the edges as a sorted list
    edges = set(cand_cand)
    for cand, lst in inc.items():
        for nb, _ in sorted(lst, key=lambda x: x[1], reverse=True)[:MAX_NEIGH]:
            edges.add((min(cand, nb), max(cand, nb)))
    return sorted(edges), f"STRING v12 (score>={STRING_MIN}, <={MAX_NEIGH} neighbors/gene)"


edges, src_label = build_string_network()
if edges is None:                                        # fallback
    edges, src_label = CURATED_EDGES, "curated SAC reference (STRING file not found)"
context = {"CDC20", "APC/C"}                         

try:
    import networkx as nx
    G = nx.Graph(); G.add_edges_from(edges)
    if G.number_of_nodes() == 0:
        print("[maps] no interactions to draw (no candidate gene found in STRING?)")
    else:
        pos = nx.spring_layout(G, seed=7, k=0.9)
        fig, ax = plt.subplots(figsize=(10, 8))
        for n in G.nodes():
            if n.upper() in cand_genes: fc, ec = "#c0392b", "black"   # variant in the proband
            elif n in context:          fc, ec = "#ecf0f1", "#95a5a6" # checkpoint target
            else:                        fc, ec = "#aed6f1", "#2980b9" # interacting partner (STRING)
            nx.draw_networkx_nodes(G, pos, nodelist=[n], node_color=fc,
                                   edgecolors=ec, node_size=2200, ax=ax)
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#7f8c8d", width=1.1, alpha=0.6)
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=8.5, font_weight="bold")
        ax.set_title(f"Protein-protein interaction network of the candidate genes\n"
                     f"red = candidate variant in the proband   |   source: {src_label}", fontsize=10)
        ax.axis("off")
        leg = [mpatches.Patch(color="#c0392b", label="candidate in the proband"),
               mpatches.Patch(color="#aed6f1", label="interacting partner")]
        ax.legend(handles=leg, loc="lower left", fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "pathway_network.png"), dpi=130)
        print(f"[maps] network ({G.number_of_nodes()} genes, {G.number_of_edges()} edges) "
              f"-> {OUT}/pathway_network.png   [{src_label}]")
except ImportError:
    print("[maps] networkx missing -> 'pip install networkx' for the network map")

# cartography
CHR_LEN = {"1":248956422,"2":242193529,"3":198295559,"4":190214555,"5":181538259,
    "6":170805979,"7":159345973,"8":145138636,"9":138394717,"10":133797422,
    "11":135086622,"12":133275309,"13":114364328,"14":107043718,"15":101991189,
    "16":90338345,"17":83257441,"18":80373285,"19":58617616,"20":64444167,
    "21":46709983,"22":50818468,"X":156040895,"Y":57227415}

def norm_chr(c): return str(c).replace("chr", "")

loci_by_chr = {}
for sym, chrom, p, tier in cand_loci:
    loci_by_chr.setdefault(norm_chr(chrom), []).append((sym, p, tier))

chroms = [c for c in CHR_LEN if c in loci_by_chr]
if chroms:
    maxlen = max(CHR_LEN[c] for c in chroms)
    fig2, ax2 = plt.subplots(figsize=(10, 0.8 * len(chroms) + 1.5))
    for i, c in enumerate(chroms):
        y = len(chroms) - i
        ax2.add_patch(mpatches.FancyBboxPatch((0, y - 0.18), CHR_LEN[c] / maxlen * 100, 0.36,
            boxstyle="round,pad=0.02,rounding_size=0.1", fc="#f2f3f4", ec="#34495e"))
        ax2.text(-2, y, f"chr{c}", ha="right", va="center", fontsize=9, fontweight="bold")
        TIER_COL = {"core": "#c0392b", "ext": "#c0392b", "phenocopy": "#e67e22", "acmg_sf": "#8e44ad"}
        for sym, p, tier in loci_by_chr[c]:
            x = p / maxlen * 100
            col = TIER_COL.get(tier, "#7f8c8d")
            ax2.plot([x], [y + 0.28], marker="v", color=col, ms=10)
            ax2.annotate(sym, (x, y + 0.45), ha="center", fontsize=8, color=col, fontweight="bold")
    ax2.set_xlim(-12, 106); ax2.set_ylim(0.3, len(chroms) + 1)
    ax2.set_title("Genomic cartography: position of the candidate variants (GRCh38)\n"
                  "red = MVA panel | orange = phenocopy | purple = ACMG secondary", fontsize=11)
    ax2.axis("off")
    fig2.tight_layout(); fig2.savefig(os.path.join(OUT, "genomic_cartography.png"), dpi=130)
    print(f"[maps] cartography -> {OUT}/genomic_cartography.png")
else:
    print("[maps] no candidate loci for the ideogram (empty candidates.tsv?)")