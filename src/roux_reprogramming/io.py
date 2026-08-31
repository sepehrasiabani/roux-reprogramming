"""Data loading and checkpoint management.

Loads the Roux et al. (2022) MSC screen (GEO GSE176206) and the intermediate
checkpoints produced by the pipeline.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import scanpy as sc


def load_raw(path: str | Path) -> ad.AnnData:
    """Load the raw screen matrix and report its value range.

    Does not normalise or log; call :func:`preprocessing.normalize_log` on raw
    counts once, downstream.

    Parameters
    ----------
    path
        Path to the ``.h5ad`` (or gzipped ``.h5ad``) file.

    Returns
    -------
    AnnData
        Raw counts in ``.X``.
    """
    adata = sc.read_h5ad(path)
    _report_state(adata)
    return adata


def load_checkpoint(path: str | Path) -> ad.AnnData:
    """Load a pipeline checkpoint."""
    return sc.read_h5ad(path)


def merge_full_genes_with_annotations(
    full_gene_path: str | Path,
    annotated_path: str | Path,
) -> ad.AnnData:
    """Combine a full-gene checkpoint with cell-type annotations.

    The full-gene checkpoint retains all genes but no annotations; the annotated
    checkpoint carries ``cell_type`` and ``age`` on a gene subset. The annotated
    cells are a subset of the full-gene cells. Subsets the full-gene object to the
    annotated barcodes, in order, and transfers the labels.

    Returns
    -------
    AnnData
        Full gene set, restricted to annotated cells, with ``cell_type`` and
        ``age`` in ``.obs``.
    """
    full = sc.read_h5ad(full_gene_path)
    annot = sc.read_h5ad(annotated_path)

    if not annot.obs_names.isin(full.obs_names).all():
        raise ValueError(
            "Annotated cells are not a subset of the full-gene checkpoint."
        )

    merged = full[annot.obs_names].copy()
    merged.obs["cell_type"] = annot.obs["cell_type"].values
    merged.obs["age"] = annot.obs["age"].values
    return merged


def _report_state(adata: ad.AnnData) -> None:
    """Print value ranges to distinguish raw counts from normalised data."""
    x = adata.X
    xmax = float(x.max())
    print(f"loaded {adata.shape[0]} cells x {adata.shape[1]} genes | X max = {xmax:.2f}")
