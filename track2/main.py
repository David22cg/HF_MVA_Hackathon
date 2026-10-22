#!/usr/bin/env python3
#  repurpose_adr.py       network (kinetic rescue + stress + dependency, SAC veto)
#  signature_reversal.py  transcriptomic reversion (CMap/L1000)
# Each .py score is z-normalised, then combined = w_network*z(net) + w_signature*z(rev).
# A drug vetoed always sinks to the bottom; 

import os, sys, csv, subprocess, argparse, yaml

BASE = os.path.dirname(os.path.abspath(__file__))


def zscores(values):
    # z-normalise a list of floats; returns {index: z}
    xs = [v for v in values if v == v]
    if len(xs) < 2:
        return {i: 0.0 for i in range(len(values))}
    mu = sum(xs) / len(xs)
    sd = (sum((x - mu) ** 2 for x in xs) / len(xs)) ** 0.5 or 1.0
    return {i: ((v - mu) / sd if v == v else 0.0) for i, v in enumerate(values)}


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description="Track 2")
    ap.add_argument("--config", default=os.path.join(BASE, "config_2.yaml"))
    ap.add_argument("--merge-only", action="store_true", help=" just merge out/*.tsv")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    OUT = os.path.expanduser(cfg["paths"].get("out", "./out")); os.makedirs(OUT, exist_ok=True)
    C = cfg.get("combine", {}); WN = float(C.get("w_network", 1.0)); WS = float(C.get("w_signature", 1.0))

    if not args.merge_only:
        for script in ("repurpose_adr.py", "signature_reversal.py"):
            print(f"[track2] running {script}")
            rc = subprocess.call([sys.executable, os.path.join(BASE, "scripts", script), args.config], cwd=BASE)
            if rc != 0:
                sys.exit(f"[track2] {script} failed (exit {rc}).")

    # key by a normalised name so the two arms match despite case/spacing
    # differences (DGIdb is upper-case, L1000 labels are lower-case).
    def norm(s):
        return (s or "").strip().lower()

    net = {norm(r["drug"]): r for r in read_tsv(os.path.join(OUT, "adr_ranking.tsv"))}
    sig = {norm(r["drug"]): r for r in read_tsv(os.path.join(OUT, "signature_reversal.tsv"))}
    if not net and not sig:
        sys.exit("[track2] no arm outputs in out/; run without --merge-only first.")

    both = set(net) & set(sig)
    print(f"[track2] drugs: {len(net)} in Arm A, {len(sig)} in Arm B, {len(both)} shared")
    if not both:
        print("  no shared names -> nothing can triangulate. Example names:")
        print("   Arm A:", list(net)[:5])
        print("   Arm B:", list(sig)[:5])

    drugs = list(dict.fromkeys(list(net) + list(sig)))      # union, order preserved

    def fnum(d, key):
        try:
            return float(d.get(key, "nan"))
        except (ValueError, TypeError):
            return float("nan")

    net_scores = [fnum(net.get(d, {}), "score") for d in drugs]
    rev_scores = [fnum(sig.get(d, {}), "reversal") for d in drugs]
    vetoed = [str(net.get(d, {}).get("veto", "")).lower() == "true" for d in drugs]
    # drop vetoed scores from the network z-norm so their -100 doesn't skew it
    net_for_z = [float("nan") if vetoed[i] else net_scores[i] for i in range(len(drugs))]
    zn = zscores(net_for_z); zs = zscores(rev_scores)

    rows = []
    for i, d in enumerate(drugs):
        nrec = net.get(d, {}).get("recommendation", "")
        scall = sig.get(d, {}).get("call", "")
        net_support = nrec.startswith("PRIORITIZE")
        sig_support = scall.startswith("REVERSES")
        sig_mimic   = scall.startswith("MIMICS")
        veto = vetoed[i]

        combined = WN * zn[i] + WS * zs[i]
        if veto:                                            # safety overrides everything
            support, final = "VETOED (SAC-acting)", "SAFETY-VETO"
            combined -= 100
        elif net_support and sig_support:                   # strongest
            support, final = "TRIANGULATED (both arms)", "PRIORITIZE"
        elif net_support and sig_mimic:                     #  temper
            support, final = "NETWORK-only (signature MIMICS - temper)", "CONSIDER (conflict)"
        elif net_support:
            support, final = "NETWORK-only", "CONSIDER"
        elif sig_support:
            support, final = "SIGNATURE-only", "CONSIDER"
        else:
            support, final = "weak/none", "LOW"

        disp = net.get(d, {}).get("drug") or sig.get(d, {}).get("drug") or d
        rows.append(dict(
            drug=disp, final_call=final, support=support, combined=round(combined, 3),
            network_rec=nrec or "-", network_score=net.get(d, {}).get("score", "-"),
            network_fdr=net.get(d, {}).get("q_stress", "-"),   # FDR of stress proximity
            reversal=sig.get(d, {}).get("reversal", "-"), signature_call=scall or "-",
            veto=veto, hub=net.get(d, {}).get("hub", "-"),
            targets=net.get(d, {}).get("targets", sig.get(d, {}).get("targets", "-")),
            reason=net.get(d, {}).get("reason", "") or sig.get(d, {}).get("reason", ""),
        ))

    # non-vetoed first, triangulated above the rest, then by combined score
    rows.sort(key=lambda r: (not r["veto"], r["support"].startswith("TRIANGULATED"), r["combined"]),
              reverse=True)

    cols = ["rank", "drug", "final_call", "support", "combined",
            "network_rec", "network_score", "network_fdr", "reversal", "signature_call",
            "hub", "veto", "targets", "reason"]
    dest = os.path.join(OUT, "repurpose_ranking.tsv")
    with open(dest, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(rows, 1):
            fh.write("\t".join([str(i)] + [str(r[c]) for c in cols[1:]]) + "\n")

    print(f"[track2] merged ranking ({len(rows)} drugs) -> {dest}")
    for i, r in enumerate(rows[:15], 1):      # top 15 only
        print(f"  #{i:<2} {r['drug']:<18} {r['final_call']:<14} {r['support']}")

    _figure(rows, OUT)


def _figure(rows, OUT):
    # scatter: x = network score, y = reversal.
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception:
        return
    pts = []
    for r in rows:
        try:
            x, y = float(r["network_score"]), float(r["reversal"])
        except (ValueError, TypeError):
            continue
        if x != x or y != y:
            continue
        if r["veto"]:
            x = max(x, -5)                       # keep vetoed visible at the left edge
        pts.append((r, x, y))
    if not pts:
        return
    label = [p for p in pts if not p[0]["veto"]][:15] + [p for p in pts if p[0]["veto"]][:5]
    fig, ax = plt.subplots(figsize=(9, 7))
    for r, x, y in pts:
        if r["veto"]:                               col, mk, a, s = "#c0392b", "x", 0.5, 45
        elif r["support"].startswith("TRIANGULATED"): col, mk, a, s = "#15803d", "*", 1.0, 160
        elif r["support"].startswith("SIGNATURE"):  col, mk, a, s = "#2563eb", "D", 0.9, 70
        elif r["final_call"].startswith("CONSIDER (conflict)"): col, mk, a, s = "#d97706", "s", 0.9, 60
        else:                                        col, mk, a, s = "#9aa0a6", "o", 0.3, 28
        ax.scatter(x, y, c=col, marker=mk, s=s, alpha=a, edgecolors="none", zorder=2)
    for r, x, y in label:
        ax.annotate(r["drug"], (x, y), fontsize=7.5, xytext=(4, 3), textcoords="offset points")
    ax.axhline(0, color="grey", lw=0.8); ax.axvline(0, color="grey", lw=0.8)
    ax.set_xlabel("Arm A network score  ->")
    ax.set_ylabel("Arm B signature reversal  (+ reverses / - mimics)  ->")
    ax.set_title("Track 2 merged: green star = both arms, blue = signature-only, "
                 "orange = conflict, red x = veto", fontsize=9.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "repurpose_map.png"), dpi=130)
    print(f"[track2] map -> {OUT}/repurpose_map.png")


if __name__ == "__main__":
    main()