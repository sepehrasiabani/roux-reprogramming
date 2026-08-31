# Methods notes

Processing and analysis choices specific to this dataset, with the reasons they
matter. The Roux screen is raw single-cell output and requires full
preprocessing before network inference.

---

## Cell-level quality control

The screen data is filtered at the cell level before any downstream step:

- **Minimum genes detected** — removes near-empty droplets.
- **Mitochondrial fraction** — high mitochondrial content indicates dying or
  leaking cells; the fraction is computed and thresholded.
- **Doublet removal** (Scrublet) — droplets containing two cells produce blended
  transcriptomes that generate spurious gene co-variation and false regulatory
  edges during network inference.

Raw sequencer output contains doublets and dying cells; retaining them yields
unreliable networks and perturbation predictions.

## Scanpy function versions

Current Scanpy function names are used throughout:

- `sc.pp.normalize_total` for total-count normalisation.
- `sc.pp.highly_variable_genes` for feature selection.

## Data state

Value ranges are checked before normalisation to confirm the input is raw
counts, since normalisation and log-transformation must be applied exactly once.
`rank_genes_groups` operates on `adata.raw` by default; `use_raw=False` is set
where scoring or differential expression should run on the processed matrix.

## Base regulatory network

The screen provides expression data only, with no chromatin-accessibility assay.
The base GRN is therefore taken from CellOracle's built-in mouse scATAC atlas
rather than constructed from sample-specific accessibility, and reflects generic
rather than cell-specific regulatory potential.

## Feature selection and the Yamanaka factors

Highly-variable-gene selection excludes Myc and Oct4, which are not among the
most variable genes. Both are retained explicitly, as they are required as
perturbation targets. Networks were rebuilt on a feature set containing all four
factors.

## Ageing axis

The ageing axis is a continuous per-cell score built from a data-derived gene
signature and supplied to CellOracle's gradient framework in place of pseudotime.
The score is imposed from a defined gene set rather than derived from the data's
geometry. Its spatial coherence on the embedding was confirmed before use.

## Perturbation magnitude

Overexpression is capped at approximately twice each factor's observed maximum
expression; CellOracle rejects larger perturbations. A dose-response (2×, 1×,
0.5×) was run to separate directional effects from artefacts of extreme values.

## Scoring input

Gene-signature scoring uses normalised, log-transformed data, not scaled data.
Scaling standardises each gene independently and removes the cross-gene magnitude
information that signature scoring depends on. Scaling is applied for PCA,
clustering, and classification, where equal per-gene contribution is appropriate.
