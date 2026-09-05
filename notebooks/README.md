# Notebooks

The analysis pipeline, one notebook per stage. Each notebook loads the checkpoint
written by the previous stage, so they run in order.

| Notebook | Contents | Output |
| `01_preprocessing.ipynb` | Loads the raw screen, reduces the object to its counts matrix, applies quality-control filtering, detects and removes doublets | `msc_clean_raw.h5ad` |
| `02_clustering_annotation.ipynb` | Normalisation, PCA, Leiden clustering, marker-gene annotation, removal of low-quality clusters | `msc_annotated.h5ad` |
| `03_grn_construction.ipynb` | Feature selection retaining the Yamanaka factors, CellOracle network inference against the mouse scATAC-seq base GRN | `oracle_roux_v2`, `links_roux_v2` |
| `04_aging_axis_and_validation.ipynb` | Data-derived ageing signature, comparison with a published senescence panel, simulation of all fifteen combinations, predicted-versus-measured validation with null models | `msc_aging_scored.h5ad`, `sim_results.pkl` |
| `05_rejuvenation_scoring.ipynb` | Ageing-axis scoring of each combination, dose-response control, per-cell-population breakdown, activation-independent axis test | `aging_axis_scores.pkl` |
| `06_robustness_and_composition.ipynb` | Combination classifier with cell-type positive control, compositional analysis of cell-type proportions | — |

Intermediate files are written to a project directory on Google Drive and are not
committed. `scripts/download_data.py` retrieves the raw data.