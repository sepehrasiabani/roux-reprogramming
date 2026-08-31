"""Quality control, normalisation, and feature selection."""

from __future__ import annotations

from collections.abc import Sequence

import anndata as ad
import scanpy as sc

# Yamanaka factors in mouse gene symbols (Oct4 = Pou5f1).
YAMANAKA_FACTORS = ("Sox2", "Pou5f1", "Klf4", "Myc")


def normalize_log(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    """Normalise to a fixed total per cell and log1p-transform.

    Operates in place on ``.X`` and returns the object. Expects raw counts.
    """
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def select_hvg(
    adata: ad.AnnData,
    n_top_genes: int = 2000,
    force_keep: Sequence[str] = YAMANAKA_FACTORS,
) -> ad.AnnData:
    """Select highly variable genes and retain the specified genes of interest.

    The Yamanaka factors are added to the highly-variable set regardless of their
    variance ranking, as they are required as perturbation targets. Expects
    normalised, log-transformed input.

    Notes
    -----
    A ``force_keep`` gene absent from ``adata.var_names`` is skipped; the caller
    should confirm all factors are present, building on a full-gene checkpoint if
    necessary.
    """
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)

    missing = [g for g in force_keep if g not in adata.var_names]
    if missing:
        print(f"force_keep genes absent from this object: {missing}")
    for gene in force_keep:
        if gene in adata.var_names:
            adata.var.loc[gene, "highly_variable"] = True

    print(f"highly-variable genes: {int(adata.var['highly_variable'].sum())}")
    return adata


def prepare_for_grn(adata: ad.AnnData) -> ad.AnnData:
    """Subset to HVG, scale, and restore raw counts to ``.X`` for CellOracle.

    Stores the full normalised, log-transformed matrix in ``.raw``, subsets to the
    highly-variable genes, scales for the internal PCA and imputation steps, and
    restores raw counts to ``.X``. Requires a ``raw_count`` layer.
    """
    if "raw_count" not in adata.layers:
        raise ValueError("expected a 'raw_count' layer")

    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    adata.X = adata.layers["raw_count"].copy()
    return adata
