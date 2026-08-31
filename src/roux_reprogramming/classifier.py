"""Combination classifier with a cell-type control.

Tests whether the Yamanaka-factor combinations are separable in expression space.
A cell-type classifier trained through the same pipeline serves as a control: if
cell type is separable but combination is not, the combinations do not carry
distinguishing signal.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.model_selection import train_test_split


def build_features(adata: ad.AnnData, n_top_genes: int = 2000) -> np.ndarray:
    """Build a normalised, HVG-subset, scaled feature matrix.

    Expects normalised, log-transformed input. Scaling is applied so that all
    genes contribute on a comparable scale.
    """
    import scanpy as sc

    feat = adata.copy()
    sc.pp.highly_variable_genes(feat, n_top_genes=n_top_genes)
    feat = feat[:, feat.var["highly_variable"]].copy()
    sc.pp.scale(feat, max_value=10)
    return feat.X


def classify(
    X: np.ndarray,
    labels: np.ndarray,
    test_size: float = 0.25,
    n_estimators: int = 200,
    random_state: int = 0,
) -> dict[str, float]:
    """Train and test a random forest; return accuracy against chance.

    Returns accuracy, balanced accuracy, chance level (1/n_classes), and the
    number of classes, using a stratified held-out test split.
    """
    n_classes = len(np.unique(labels))
    x_tr, x_te, y_tr, y_te = train_test_split(
        X, labels, test_size=test_size, random_state=random_state, stratify=labels
    )
    clf = RandomForestClassifier(
        n_estimators=n_estimators, random_state=random_state, n_jobs=-1
    )
    clf.fit(x_tr, y_tr)
    pred = clf.predict(x_te)
    return {
        "accuracy": float(accuracy_score(y_te, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_te, pred)),
        "chance": 1.0 / n_classes,
        "n_classes": n_classes,
    }


def balance_classes(labels: np.ndarray, random_state: int = 0) -> np.ndarray:
    """Return indices sampling equal cells per class, down to the smallest class."""
    rng = np.random.RandomState(random_state)
    classes, counts = np.unique(labels, return_counts=True)
    per_class = counts.min()
    idx = []
    for cls in classes:
        cls_idx = np.where(labels == cls)[0]
        idx.extend(rng.choice(cls_idx, per_class, replace=False))
    return np.array(idx)
