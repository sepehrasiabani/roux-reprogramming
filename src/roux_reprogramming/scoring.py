"""Ageing-axis construction and perturbation scoring.

Constructs a continuous ageing axis from a gene signature and scores each
simulated perturbation by its inner product with that axis. A negative inner
product indicates movement toward a younger state; a positive value, movement
toward an aged state.
"""

from __future__ import annotations

from collections.abc import Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.interpolate import griddata


def build_aging_score(
    adata: ad.AnnData,
    signature: Sequence[str],
    score_name: str = "aging_score",
) -> ad.AnnData:
    """Score each cell on an ageing gene signature.

    Uses ``sc.tl.score_genes``, which subtracts a matched random background gene
    set; the score is therefore centred near zero and may be negative. Expects
    normalised, log-transformed (not scaled) input.
    """
    present = [g for g in signature if g in adata.var_names]
    if len(present) < len(signature):
        missing = set(signature) - set(present)
        print(f"{len(missing)} signature genes absent: {sorted(missing)}")
    sc.tl.score_genes(adata, gene_list=present, score_name=score_name, use_raw=False)
    return adata


def validate_axis_separates_age(
    adata: ad.AnnData,
    score_key: str = "aging_score",
    age_key: str = "age",
) -> float:
    """Return (aged mean − young mean) for the score; positive tracks ageing."""
    young = adata.obs.loc[adata.obs[age_key] == "Young", score_key]
    aged = adata.obs.loc[adata.obs[age_key] == "Aged", score_key]
    diff = float(aged.mean() - young.mean())
    print(f"{score_key}: young={young.mean():+.3f} aged={aged.mean():+.3f} diff={diff:+.3f}")
    return diff


def perturbation_inner_product_per_cell(oracle, gradient) -> np.ndarray:
    """Per-cell inner product of the perturbation shift with the ageing gradient.

    Interpolates the gradient's reference flow to each cell's embedding position
    and takes the dot product with that cell's ``delta_embedding``. Requires that
    ``simulate_shift``, ``estimate_transition_prob``, and
    ``calculate_embedding_shift`` have been run on ``oracle``.

    Returns
    -------
    np.ndarray
        One inner-product value per cell (negative: toward younger).
    """
    shift = oracle.delta_embedding
    grid_points = gradient.gridpoints_coordinates
    grid_flow = gradient.ref_flow
    cell_coords = oracle.adata.obsm["X_umap"]

    flow_x = griddata(grid_points, grid_flow[:, 0], cell_coords, method="nearest")
    flow_y = griddata(grid_points, grid_flow[:, 1], cell_coords, method="nearest")
    aging_dir = np.stack([flow_x, flow_y], axis=1)

    return np.sum(shift * aging_dir, axis=1)


def score_by_cell_type(
    oracle,
    ip_per_cell: np.ndarray,
    cell_type_key: str = "cell_type",
) -> pd.Series:
    """Mean ageing-axis inner product within each cell type, sorted."""
    df = pd.DataFrame(
        {"cell_type": oracle.adata.obs[cell_type_key].values, "ip": ip_per_cell}
    )
    return df.groupby("cell_type", observed=True)["ip"].mean().sort_values()
