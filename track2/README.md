# Track 2. Drug repurposing for Mosaic Variegated Aneuploidy

This directory contains an offline, scalable drug-repurposing model for Mosaic Variegated Aneuploidy. The model ranks approved and investigational drugs using two independent mechanism-based methods and a safety constraint derived from the disease biology. 
All outputs arehypotheses for experimental follow-up.

## Method

The causal variant weakens the spindle assembly checkpoint, which allows cells to
enter anaphase before the chromosomes are correctly attached and produces mosaic
aneuploidy. Because the defect is a loss of function in a structural checkpoint
protein, the model does not attempt to restore the protein. It targets the
downstream consequences of aneuploidy.

The first method measures network proximity on the STRING protein interaction
network reused from Track 1. For each drug it computes the distance from the drug
targets to a module of druggable aneuploidy-stress genes and compares it against a
degree-preserving random model, which yields a z-score and an empirical p-value.
Drugs are also rewarded when their targets slow the cell cycle, such as CDK4 and
CDK6 or mTOR inhibitors, because a longer cycle gives the error-correction
machinery more time without inhibiting the checkpoint, and when their targets are
synthetic-lethal dependencies of aneuploid cells such as KIF18A.

The second method measures transcriptomic reversion. A public aneuploidy-stress
expression signature defines which genes are increased and decreased in aneuploid
cells. Each drug is scored by how strongly its LINCS L1000 expression signature
reverses that disease signature. A drug that lowers the increased genes and raises
the decreased genes is a candidate. A drug that reinforces the signature is
flagged as unsuitable.

A safety constraint overrides both scores. Any drug that acts on the checkpoint
itself is vetoed, because further weakening mitotic fidelity in a patient whose checkpoint is already impaired would increase missegregation.

The two scores are z-normalised and combined. A drug supported by both methods is
labelled TRIANGULATED and carries the strongest evidence. Multiple-testing
correction is applied with the Benjamini-Hochberg procedure. The model also
reports how many literature compounds with established aneuploidy selectivity it
recovers, as an external validation.

The model is unsupervised and mechanism-based. Mosaic Variegated Aneuploidy affects fewer than fifty patients worldwide and has no labelled drug-outcome data, so a trained classifier would overfit. Network proximity and signature reversal require no training labels and every score is
explainable.

To apply the model to another disease, replace the module gene lists and the
disease signature. 

## Folder layout

```
track2/
  main.py                      orchestrator that runs both methods and merges them
  config_2.yaml           paths, weights and thresholds
  scripts/
    from_track1.py             builds the disease seed from the Track 1 result
    repurpose_adr.py           network-proximity method
    signature_reversal.py      transcriptomic-reversion method
    drug_targets.py      builds the real drug-target table from DGIdb
    drug_signatures.py   builds the real drug signatures from LINCS L1000
  modules/
    stress_module.txt          druggable aneuploidy-stress genes, rewarded
    kinetic_module.txt         cell-cycle-slowing genes, rewarded
    checkpoint_module.txt      checkpoint genes, vetoed
    aneuploidy_dependencies.txt synthetic-lethal dependencies, rewarded
  data/
    demo_drug_targets.tsv      small demo drug-target table
    drug_signatures.tsv        small synthetic demo signature table
    disease_signature.tsv      public aneuploidy-stress signature
    known_positives.txt        literature positive controls for validation
  out/                         generated rankings and figures
```

## Data sources

All inputs are public. STRING is the same file used in Track 1.

| Input | Source |
|-------|--------|
| STRING protein interactions | 9606.protein.links and info v12.0, combined score at least 700 |
| Drug targets | DGIdb interactions table from dgidb.org/downloads |
| Drug signatures | LINCS L1000 from Enrichr |
| Disease signature | literature aneuploidy-stress signature, editable |
| Gene modules | literature-curated, sources cited in each file header |

The demo tables in data are small and illustrative.

## Running the model

Install the dependencies in the same environment as Track 1.

```bash
pip install networkx numpy matplotlib pyyaml
cd ../track2
```

Build the real inputs and run.

```bash
python3 scripts/from_track1.py                 # writes data/seed_genes.txt from Track 1

python3 scripts/drug_targets.py interactions.tsv

base='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName='
wget -O lincs_up.gmt   "${base}LINCS_L1000_Chem_Pert_up"
wget -O lincs_down.gmt "${base}LINCS_L1000_Chem_Pert_down"
python3 scripts/build_drug_signatures.py lincs_up.gmt lincs_down.gmt \
  --name-regex '^[^-]*-(.+)-[0-9.]+$'

python3 main.py
```

Use python3 main.py --merge-only to re-merge existing results without recomputing. 
Change paths in config if needed.

## Outputs

| File | Content |
|------|---------|
| out/repurpose_ranking.tsv | final merged ranking with the final call, the support category, the combined score, the network false-discovery rate, the reversal score, the targets and the reason |
| out/adr_ranking.tsv | network-proximity results with the score, z-score, p-value, false-discovery rate and the kinetic, dependency and veto flags |
| out/signature_reversal.tsv | reversion results with the reversal score and the call of reverses, mimics or neutral |
| out/repurpose_map.png | merged map with the network score on the horizontal axis and the reversal score on the vertical axis |
| out/adr_map.png | network map with closeness to the stress module and closeness to the checkpoint |
| out/validation.txt | how many literature positive controls the model recovers |

In repurpose_ranking.tsv the support column is the main result. TRIANGULATED means
both methods agree. NETWORK-only and SIGNATURE-only mean a single method supports
the drug. VETOED means the drug acts on the checkpoint and is excluded for safety.


## Limitations

There is no direct drug-efficacy data for Mosaic Variegated Aneuploidy, so every
candidate is a hypothesis for follow-up. Network proximity captures drugs with
diffuse or metabolic mechanisms poorly, which is one reason the signature method is
added. The gene modules and the disease signature are literature-curated and
editable, and the results depend on them. STRING and DGIdb are biased toward
well-studied genes.