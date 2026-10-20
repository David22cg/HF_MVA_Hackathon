#!/usr/bin/env python3
# inspect a vcf.gz and report its data 
# USAGE:  python3 00_inspect.py <vcf.gz>
import sys, gzip, collections

# chr1 has a different length in each build 
CHR1_LEN = {248956422: "GRCh38", 249250621: "GRCh37/hg19"}
MAIN_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
# Common callers
KNOWN_CALLERS = ("gatk", "haplotypecaller", "dragen", "deepvariant", "strelka", "mutect")


#header
def parse_header(fh):
    h = dict(samples=[], callers=set(), info_ids=set(), chr1_len=None)
    for line in fh:
        if not line.startswith("#"):
            break                                   # first variant: end of header
        low = line.lower()
        # chr1 length declared in the contigs -> build
        if line.startswith("##contig") and ("id=chr1," in low or "id=1," in low):
            try:
                h["chr1_len"] = int(line.split("length=")[1].split(">")[0].split(",")[0])
            except Exception:
                pass
        elif line.startswith("##INFO=<ID="):        # catalogue of INFO fields (gnomAD, etc.)
            h["info_ids"].add(line.split("ID=")[1].split(",")[0])
        for t in KNOWN_CALLERS:
            if t in low:
                h["callers"].add(t)
        if line.startswith("#CHROM"):               # columns 10+ = sample names
            h["samples"] = line.rstrip("\n").split("\t")[9:]
    return h


# body one variant per line
def scan_variants(fh):
    s = dict(n_sites=0, n_multiallelic=0, n_indel_alleles=0,
             first_chrom=None, filters=collections.Counter(),
             format_fields=set(), per_chrom=collections.Counter())
    for line in fh:
        if line.startswith("#"):                    
            continue
        f = line.split("\t", 9)                    
        if len(f) < 9:
            continue
        chrom, ref, alt, flt, fmt = f[0], f[3], f[4], f[6], f[8]
        if s["first_chrom"] is None:
            s["first_chrom"] = chrom
        s["per_chrom"][chrom.replace("chr", "")] += 1
        s["filters"][flt] += 1
        s["format_fields"].update(fmt.split(":"))
        alts = alt.split(",")
        if len(alts) > 1:
            s["n_multiallelic"] += 1                # site with >1 ALT (must be split)
        for a in alts:
            if len(ref) != len(a):
                s["n_indel_alleles"] += 1           # different lengths = indel
        s["n_sites"] += 1
    return s


#report

def report(h, s):
    chr_prefix = str(s["first_chrom"]).startswith("chr")
    ad = "AD" in s["format_fields"]
    dp = "DP" in s["format_fields"]
    gt = "GT" in s["format_fields"]
    build = CHR1_LEN.get(h["chr1_len"], f"undetermined (chr1={h['chr1_len']})")

    print("=" * 66)
    print(f"Sample(s)            : {h['samples']}   -> use as proband_id")
    print(f"Caller               : {', '.join(sorted(h['callers'])) or 'not evident'}")
    print(f"First chromosome     : {s['first_chrom']}")
    print(f"'chr' prefix         : "
          + ("YES -> 01_prep strips it for annotation (05 re-adds it to the CSV)"
             if chr_prefix else "NO -> already Ensembl style, no renaming needed"))
    print(f"GT / AD / DP         : GT={gt}  AD={ad}  DP={dp}   "
          + ("(OK: VAF available -> mosaic rescue + BAF map)"
             if ad and dp else "( AD/DP MISSING: no VAF, no mosaic, no BAF)"))
    print(f"FILTER               : {dict(s['filters'].most_common())}")
    if set(s["filters"]) <= {"."}:
        print("   -> UNFILTERED VCF (FILTER='.'): do NOT use '-f PASS' or everything is deleted.")
    elif "PASS" in s["filters"]:
        nonpass = sum(v for k, v in s["filters"].items() if k not in ("PASS", "."))
        print(f"   -> {nonpass} non-PASS present: may hide MOSAICS (low VAF) ")
    print(f"Total variants       : {s['n_sites']}")
    print(f"Multiallelic sites   : {s['n_multiallelic']}  ")
    print(f"Indel alleles        : {s['n_indel_alleles']}   ")
    print(f"INFO tags            : {sorted(h['info_ids'])}")

    # per-chromosome integrity: no main chromosome should be empty
    missing = [c for c in MAIN_CHROMS if s["per_chrom"].get(c, 0) == 0]
    print("-" * 66)
    print("Variants per main chromosome (loss detection):")
    print("   " + "  ".join(f"{c}:{s['per_chrom'].get(c, 0)}" for c in MAIN_CHROMS))
    if missing:
        print(f" no variants in: {missing}  -> the VCF may come in truncated -> check")
    else:
        print("   OK: every main chromosome has variants ")
    print("=" * 66)


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 00_inspect.py <vcf.gz>")
    vcf = sys.argv[1]
    with gzip.open(vcf, "rt") as fh:
        h = parse_header(fh)
    with gzip.open(vcf, "rt") as fh:
        s = scan_variants(fh)
    report(h, s)


if __name__ == "__main__":
    main()