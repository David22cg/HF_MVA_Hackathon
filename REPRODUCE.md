# GenomEx

Two tracks, fully offline. 
Track 1 predicts the causal variant. Track 2 proposes repurposed drugs from it.

---

## Folder layout

The code repo:

```
mva_hackathon/
  track1/
    main.py                  # Track 1 orchestrator
    scripts/                 # Track 1 helper (acmg_audit.py)
    work/                    # intermediate files + candidates.tsv  (created on run)
    out/                     # submission CSV + figures             (created on run)
    config.yaml                # Track 1 settings and paths
    hpo_terms.txt              # patient HPO terms, one HP:xxxxxxx per line (you add this)
    genes/                     # panel, phenocopy and ACMG gene lists (included)
    scripts/                   # Track 1 stages 00..06 (included)
  track2/
    main.py                  # Track 2 orchestrator
    config_track2.yaml       # Track 2 settings and paths
    scripts/                 # Arm A, Arm B, build tables, bridge (included)
    modules/                 # gene modules: stress, kinetic, checkpoint, dependencies 
    data/                    # drug tables + disease signature 
    out/                     # repurpose_ranking.tsv + map          (created on run)
    rankvar/                   # optional independent check (included)
```

The reference data :

```
~/mva/
  ref/     # genome FASTA
  cache/   # VEP cache
  data/    # ClinVar, AlphaMissense, SpliceAI, LOFTEE, HPO, STRING
```

---

## 2. One-time downloads

Put them where `config.yaml` points (`~/mva/ref`, `~/mva/cache`, `~/mva/data`).

Track 1 needs: genome FASTA (Ensembl GRCh38), VEP cache (release 116), ClinVar, AlphaMissense, SpliceAI, LOFTEE, the HPO gene-to-phenotype file, and STRING (`9606.protein.links` + `.info`).

Track 2 needs: STRING (same files), the DGIdb interactions export, and the LINCS L1000 up/down GMT pair.

Full download commands are in `RUN.md`.

---

## 3. Before running 

Track 1 — `config.yaml`:
- `proband_id: ""`  
- `hf_username`  
- `paths.vcf_raw`  (and the other paths)
- add `hpo_terms.txt` with the present HPO terms

Track 2 — `track2/config_2.yaml`:
- `paths.track1_candidates: ../track1/work/candidates.tsv`
- `paths.drug_targets` and `paths.drug_signatures` point to your tables 
- `allow_demo: false`  
- `checkpoint.from_seed: false` 

---

## 4. Execution order

**4a. Track 1 — variant prediction**
```bash
python3 track1/main.py          # runs stages 0..7
```
Stages: 0 inspect, 1 normalize, 2 annotate (VEP), 3 panel audit, 4 filter and rank, 5 BAF/LOH map, 6 submission CSV, 7 network maps.
Resume after a fix with `python3 track1/main.py --from <stage>`.

Outputs:
- `track1/work/candidates.tsv`  (ranked candidates)
- `track1/out/<hf_username>_<approach>.csv`  (submission file)

**4b. Track 2 — build the real drug tables (one time)**
```bash
cd track2
python3 scripts/build_drug_targets.py interactions.tsv
python3 scripts/build_drug_signatures.py lincs_up.gmt lincs_down.gmt \
        --name-regex '^[^-]*-(.+)-[0-9.]+$'
```
This writes `data/drug_targets_real.tsv` and `data/drug_signatures_real.tsv`.
Point the two `paths.*` in `config_track2.yaml` to them.

**4c. Bridge — Track 1 result becomes the Track 2 seed**
```bash
cd track2
python3 from_track1.py          # reads ../track1/work/candidates.tsv, writes data/seed_genes.txt
```
Only gene symbols are read. No coordinates, genotypes or phenotype.

**4d. Track 2 — drug repurposing**
```bash
cd track2
python3 main.py                 # runs Arm A + Arm B, then merges
```
Output:
- `track2/out/repurpose_ranking.tsv`  (ranked drugs)
- `track2/out/repurpose_map.png`  (the two-method map)

**4e. Optional — independent check (RankVar)**
```bash
bash rankvar/run_rankvar.sh
python3 rankvar/compare_rankvar.py
```

---