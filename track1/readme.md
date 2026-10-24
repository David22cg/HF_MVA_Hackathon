# Variant Prediction

Offline pipeline that predicts the causal variant for MVA

## Folder structure

```
└── track1/
    ├── main.py                 # one-command orchestrator (stages 0–7)
    ├── config.yaml             # all paths + thresholds (each one sourced)
    ├── hpo_terms.txt           # HPO symptoms (one HP:xxxxxxx per line)
    ├── scripts/
    │   ├── 00_inspect.py       # QC of the raw VCF
    │   ├── 00b_panel_audit.py  # panel-coverage sanity check
    │   ├── 01_prep.sh          # normalize (bcftools)                 [env: mva]
    │   ├── 02_annotate.sh      # VEP + AlphaMissense (+SpliceAI/LOFTEE)[env: vep]
    │   ├── 03_filter_rank.py   # filter + evidence-tiered ranking
    │   ├── 04_baf_map.py       # BAF/LOH map + mosaic-aneuploidy burden
    │   ├── 05_make_csv.py      # submission CSV
    │   ├── 06_maps.py          # STRING network + genomic cartography
    │   ├── make_spliceai.sh    # compute SpliceAI locally
    │   └── acmg_audit.py       # diagnose why no ACMG secondary finding appears
    ├── genes/
    │   ├── mva_g.txt           # MVA panel (core + SAC/kinetochore ext)
    │   ├── phenocopy.txt       # differential-diagnosis genes
    │   └── acmg_sf.txt         # ACMG secondary-findings genes
    ├── rankvar/
    │   └── compare_rankvar.py  # comparison of RankVar result and our ranking as benchmark
    ├── work/                   # (generated) intermediate files
    └── out/                    # (generated) submission CSV + figures
```

## Run

```bash
cd ../track1
python3 main.py                 # full pipeline, stages 0–7
# resume options:
python3 main.py --list          # show stages
python3 main.py --from 3        # skip prep/annotate if annotated.vcf.gz exists
python3 main.py --from 4 --to 7 # just ranking to maps
python3 main.py --no-conda      # run in the current env (don't switch via `conda run`)
```

`main.py` runs each stage in its conda env (`mva` / `vep`) via `conda run`, from
this folder so the config's relative paths (`./genes`, `./hpo_terms.txt`) resolve.
Edit the paths inside `01_prep.sh` / `02_annotate.sh` (raw VCF, FASTA, VEP cache)
and the `config.yaml` once before the first run.

## Independent comparison 

RankVar is run separately as an independent, phenotype-driven check. After the
pipeline produces `work/candidates.tsv`, compare the two rankings at gene and
variant level.

Installation and execution instructions of this tool here: 

```bash
git clone https://github.com/WGLab/RankVar.git
```
```bash
cd rankvar
python3 compare_rankvar.py rank_var.tsv ~/mva/work/candidates.tsv ../genes/mva_g.txt
```

## Setup

### 1) Software
```bash
# env for data wrangling + the python/bash scripts
conda create -n mva -c conda-forge -c bioconda \
  bcftools samtools htslib python=3.11 pandas numpy matplotlib pyyaml networkx -y

# env for VEP (keep channel order: conda-forge BEFORE bioconda)
conda create -n vep -c conda-forge -c bioconda "ensembl-vep>=116" -y

# VEP plugins repo (provides AlphaMissense.pm and SpliceAI.pm)
git clone https://github.com/Ensembl/VEP_plugins ~/.vep/Plugins
```

### 2) Reference data

**2a) FASTA**
```bash
cd ~/mva/ref
wget -c http://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -t Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz   # integrity check before gunzip
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
conda activate mva && samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

**2b) VEP cache**
```bash
cd ~/mva/cache
wget -c http://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
tar xzf homo_sapiens_vep_116_GRCh38.tar.gz          # extracts homo_sapiens/116_GRCh38/
```

**2c) ClinVar**
```bash
cd ~/mva/data
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

**2d) AlphaMissense plugin data**
```bash
cd ~/mva/data
wget -c https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
( zcat AlphaMissense_hg38.tsv.gz | head -n4 | grep '^#' ; \
  zcat AlphaMissense_hg38.tsv.gz | grep -v '^#' | sed 's/^chr//' ) \
  | bgzip > AlphaMissense_hg38.ens.tsv.gz
tabix -s1 -b2 -e2 -f AlphaMissense_hg38.ens.tsv.gz
```

**2e) SpliceAI plugin data**
```bash
conda activate vep
bash scripts/make_spliceai.sh      # -> ~/mva/data/spliceai.hg38.vcf.gz (+ .tbi)
```
> If the `spliceai` package complains: `np.fromstring` removed, patch
> `spliceai/utils.py` to `np.frombuffer(bytes(seq,"utf-8"), np.int8)`. If
> `pkg_resources` is missing, `pip install "setuptools<81"`. SpliceAI is optional
> and the pipeline still runs without it.

**2f) LOFTEE plugin data** (GRCh38, LoF confidence HC/LC)
```bash
cd ~/.vep/Plugins
git clone -b grch38 https://github.com/konradjk/loftee        # GRCh38 branch
cd ~/mva/data
wget -c https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz
wget -c https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.fai
wget -c https://personal.broadinstitute.org/konradk/loftee_data/GRCh38/human_ancestor.fa.gz.gzi
```
> GRCh38 LOFTEE also needs the GERP conservation bigwig and the conservation SQL.
> See the loftee GRCh38 README for the current URLs, then set `GERP_BIGWIG` and
> `CONSERVATION_FILE` in `02_annotate.sh`. LOFTEE is optional.

**2g) HPO gene to phenotype map**
```bash
cd ~/mva/data
wget -c https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt    # or from hpo.jax.org/data/download
```

**2h) STRING**
```bash
cd ~/mva/data
wget -c https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
wget -c https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
```

**2i) Patient HPO terms** — create `hpo_terms.txt`, one present-symptom HP id per line:
```
HP:XXXXXXX
...
```

## Outputs

Generated into `work/` (intermediate) and `out/`:

| File | Explanation |
|------|-------------|
| `out/<username>_<approach>.csv` | ranked causal-variant predictions with a calibrated probability per row |
| `out/loh_scan.tsv` | cryptic 2nd-allele (LOH) scan per candidate gene: heterozygosity-depletion vs chromosome density, depth ratio, allelic imbalance, and a flag (`SUSPECTED_BIALLELIC` / repetitive-gain / no-signal) |
| `out/aneuploidy_burden.tsv` | per-chromosome mosaic-aneuploidy burden, estimated mosaic fraction from the deviation of BAF away from 0.5 |
| `out/baf_map.png` | genome-wide B-allele-frequency map (visual of the mosaic aneuploidy) |
| `out/depth_karyotype.png` | per-chromosome read-depth karyotype |
| `out/pathway_network.png` | STRING PPI network of the candidate genes plus top neighbours |
| `out/genomic_cartography.png` | ideogram with the candidate variants placed on the chromosomes |
| `work/all.norm.vcf.gz` | normalized, left-aligned, split-multiallelic VCF |
| `work/annotated.vcf.gz` | VEP annotation: consequence, gnomAD AF, ClinVar, AlphaMissense (+SpliceAI/LOFTEE) |
| `work/candidates.tsv` | filtered and evidence-tier-ranked variant table before CSV formatting |
| `work/het_snps.tsv` | cached heterozygous-SNP table used by the LOH and BAF step |
