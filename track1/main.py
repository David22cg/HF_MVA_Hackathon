#   python3 main.py                      # run the full pipeline
#   python3 main.py --from 3             # resume from stage 3 
#   python3 main.py --from 4 --to 7      # just the ranking->maps analysis stages
#   python3 main.py --config config.yaml --vcf /mnt/c/.../raw.vcf.gz
#   python3 main.py --no-conda           # run in the current env
#   python3 main.py --list               # list the stages and exit
# =====================================================================
import os, sys, subprocess, shutil, argparse, yaml

BASE = os.path.dirname(os.path.abspath(__file__))          # project root 
SCRIPTS = os.path.join(BASE, "scripts")

def stage_list(cfg_path, raw_vcf, work, out):
    def py(script, *a): return ["python3", os.path.join(SCRIPTS, script), *a]
    def sh(script):     return ["bash",    os.path.join(SCRIPTS, script)]
    return [
        (0, "inspect",      "mva", py("00_inspect.py", raw_vcf) if raw_vcf else None, None),
        (1, "prep",         "mva", sh("01_prep.sh"),                 os.path.join(work, "all.norm.vcf.gz")),
        (2, "annotate",     "vep", sh("02_annotate.sh"),            os.path.join(work, "annotated.vcf.gz")),
        (3, "panel_audit",  "mva", py("00b_panel_audit.py", cfg_path), None),
        (4, "filter_rank",  "mva", py("03_filter_rank.py", cfg_path), os.path.join(work, "candidates.tsv")),
        (5, "baf_map",      "mva", py("04_baf_map.py", cfg_path),     os.path.join(out,  "loh_scan.tsv")),
        (6, "make_csv",     "mva", py("05_make_csv.py", cfg_path),    None),
        (7, "maps",         "mva", py("06_maps.py", cfg_path),        None),
    ]

def main():
    ap = argparse.ArgumentParser(description="Track 1 MVA pipeline orchestrator")
    ap.add_argument("--config", default=os.path.join(BASE, "config.yaml"))
    ap.add_argument("--vcf", default=None, help="raw VCF for stage 0 (default: paths.vcf_raw in config)")
    ap.add_argument("--from", dest="frm", type=int, default=0, help="first stage id to run")
    ap.add_argument("--to", type=int, default=7, help="last stage id to run")
    ap.add_argument("--no-conda", action="store_true", help="run in the current env, no `conda run`")
    ap.add_argument("--list", action="store_true", help="list stages and exit")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    work = os.path.expanduser(cfg["paths"]["work"])
    out  = os.path.expanduser(cfg["paths"]["out"])
    raw_vcf = args.vcf or cfg["paths"].get("vcf_raw")
    stages = stage_list(args.config, raw_vcf, work, out)

    if args.list:
        print("Track 1 stages:")
        for sid, name, env, cmd, _ in stages:
            print(f"  {sid}  {name:<12} [env: {env}]")
        return

    use_conda = (not args.no_conda) and shutil.which("conda") is not None
    if not args.no_conda and not use_conda:
        print("[main] 'conda' not found on PATH -> running in the current env (make sure it has the deps).")

    print(f"[main] project root: {BASE}")
    print(f"[main] config: {args.config}   stages {args.frm}..{args.to}\n")

    for sid, name, env, cmd, produced in stages:
        if sid < args.frm or sid > args.to:
            continue
        if cmd is None:
            print(f"[main] stage {sid} ({name}) skipped (no raw VCF given; set paths.vcf_raw or pass --vcf)\n")
            continue
        full = (["conda", "run", "--no-capture-output", "-n", env] + cmd) if use_conda else cmd
        print(f"[main] STAGE {sid}: {name}   (env: {env})")
        rc = subprocess.call(full, cwd=BASE)                 # cwd=root so ./genes, ./hpo_terms.txt resolve
        if rc != 0:
            sys.exit(f"[main] stage {sid} ({name}) FAILED (exit {rc}). Fix it and resume with --from {sid}.")
        if produced and not os.path.exists(produced):
            sys.exit(f"[main] stage {sid} ({name}) finished but expected output is missing: {produced}")
        print(f"[main] stage {sid} ({name}) OK" + (f" -> {produced}" if produced else "") + "\n")

    print("done")

if __name__ == "__main__":
    main()