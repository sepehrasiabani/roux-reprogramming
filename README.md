# Predicting and Scoring Partial Reprogramming: A Gene-Regulatory-Network Reanalysis of a Pooled Yamanaka-Factor Screen

> **Status:** analysis complete; manuscript in preparation

A reanalysis of a pooled single-cell Yamanaka-factor screen (Roux et al., 2022) in young and aged mouse mesenchymal stromal cells, addressing two questions: whether a gene regulatory network (GRN) model can predict the transcriptional effect of partial reprogramming when validated against the measured screen, and whether simulated reprogramming moves cells toward a younger state along a data-derived ageing axis.

## Scientific question

Partial reprogramming — transient expression of the Yamanaka factors (Sox2, Oct4/Pou5f1, Klf4, Myc) — can reverse features of cellular ageing without erasing cell identity, but which factor combinations do so, and by what regulatory logic, is unresolved. Roux et al. (2022) measured a pooled screen of factor combinations in aging mouse mesenchymal stromal cells (MSCs). This analysis asks whether a GRN model built from that data can (1) predict the measured effect of each combination, and (2) be used to score whether reprogramming moves cells against a data-derived ageing axis.

## Key findings

- Predicted per-combination expression changes correlate with the measured screen (mean r = 0.13; 78% of combination × cell-type comparisons positive), significantly exceeding a gene-shuffled null (p = 2×10⁻¹¹) but not a combination-mismatched null (p = 0.47): the model recovers the genes that respond to reprogramming but does not resolve combination-specific effects.
- The factor combinations converge on a shared transcriptional response. Measured combination effects correlate r ≈ 0.83 with one another; a supervised classifier separated cell types at 86% accuracy but combinations at only 11% (chance 7%).
- All fifteen combinations scored positive on the ageing axis (movement toward an aged, activated state rather than a younger one), robust across overexpression doses and uniform across cell populations, driven predominantly by Myc.
- Ageing and myofibroblast activation are inseparable in these cells: an ageing axis constructed to exclude activation genes remained correlated with the activation-associated axis (r = 0.74) and produced the same result.
- A published senescence (SASP) panel did not track ageing in these cells; the data-derived ageing signature is a myofibroblast/activation program, corroborated by a compositional shift toward stressed and senescent populations with age.

## Approach

| Stage | Method | Notebook |
|---|---|---|
| Preprocessing and QC | Scanpy; Scrublet doublet removal | `notebooks/01_preprocessing.ipynb` |
| Clustering and annotation | Leiden clustering; marker-gene validation | `notebooks/02_clustering_annotation.ipynb` |
| GRN inference | CellOracle (mouse base GRN; per-cluster regression) | `notebooks/03_grn_construction.ipynb` |
| Perturbation and validation | In-silico simulation; correlation against the measured screen; null models | `notebooks/04_aging_axis_and_validation.ipynb` |
| Ageing axis and scoring | Data-derived ageing signature; gradient-based perturbation scoring | `notebooks/05_rejuvenation_scoring.ipynb` |
| Robustness | Activation-independent ageing axis; combination classifier | `notebooks/06_robustness_and_composition.ipynb` |

The cell populations are parallel identities rather than stages of a differentiation trajectory; analysis uses group contrasts (young versus aged, treated versus untreated) rather than trajectory inference. See `docs/decisions.md`.

## Reproducing this analysis

```bash
conda env create -f environment.yml
conda activate roux-reprogramming
pip install -e .
python scripts/download_data.py
```

GRN inference (`get_links`) is the longest step (~25 min); scoring across all fifteen combinations is ~50 min. The pipeline is CPU-bound.

## Repository layout

```
src/roux_reprogramming/   importable package
  io.py                   data loading and checkpoints
  preprocessing.py        QC, normalisation, feature selection
  propagation.py          GRN signal-propagation implementation
  scoring.py              ageing-axis construction and perturbation scoring
  classifier.py           combination classifier and cell-type control
  plotting.py             figure functions
notebooks/                analysis notebooks, one per stage
scripts/                  command-line entry points
tests/                    unit tests
figures/                  generated figures
data/                     gitignored; populated by scripts/download_data.py
docs/                     methods notes and decision log
```

## Data

- Source: Roux, A. E. et al. Diverse partial reprogramming strategies restore youthful gene expression and transiently suppress cell identity. *Cell Systems* 13, 574–587 (2022). doi:10.1016/j.cels.2022.05.002
- Accession: GEO GSE176206 (pooled Yamanaka-factor screen in young and aged mouse MSCs; expression only).
- Raw data is not committed; `scripts/download_data.py` retrieves it.

## Licence

MIT — see `LICENSE`.
