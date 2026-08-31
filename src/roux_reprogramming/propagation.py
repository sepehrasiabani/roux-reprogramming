"""GRN signal propagation.

Implements the perturbation-propagation step of the CellOracle algorithm from
the published formulation. Correctness is verified against the library's output
in ``tests/test_propagation.py``.

Method
------
A fitted gene regulatory network gives, for each target gene, a linear model
predicting its expression from its regulators::

    expression[target] ~ sum_over_regulators( coefficient[regulator, target]
                                              * expression[regulator] )

Collecting the coefficients into a matrix ``W`` (regulators x targets) expresses
one round of regulatory influence as a matrix product. A perturbation propagates
outward in rounds: the perturbed gene shifts its direct targets, those targets
shift their own targets, and so on.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["propagate_perturbation", "simulate_knockout"]


def propagate_perturbation(
    delta: npt.NDArray[np.float64],
    coefficient_matrix: npt.NDArray[np.float64],
    clamped_indices: npt.NDArray[np.int_],
    clamped_values: npt.NDArray[np.float64],
    n_propagation: int = 3,
) -> npt.NDArray[np.float64]:
    """Propagate an expression perturbation through a regulatory network.

    Parameters
    ----------
    delta
        Initial expression change per gene, shape ``(n_cells, n_genes)``.
        Non-zero only at the perturbed genes on entry.
    coefficient_matrix
        Fitted GRN coefficients, shape ``(n_genes, n_genes)``, oriented so that
        ``coefficient_matrix[i, j]`` is the effect of regulator ``i`` on target
        ``j``. Positive values activate, negative values repress.
    clamped_indices
        Column indices of genes held fixed at their perturbed value throughout
        propagation (a knocked-out gene must not be revived by its own
        regulators).
    clamped_values
        Values to hold the clamped genes at, shape ``(n_cells, len(clamped_indices))``.
    n_propagation
        Number of propagation rounds. Each round lets influence travel one
        further step through the network. Defaults to 3, as in CellOracle.

    Returns
    -------
    numpy.ndarray
        Propagated expression change per gene, shape ``(n_cells, n_genes)``.

    Raises
    ------
    ValueError
        If array shapes are inconsistent or ``n_propagation`` is negative.

    Notes
    -----
    Re-clamping after every round is essential and easy to get wrong: without
    it the perturbed gene is pulled back toward its unperturbed value by its
    own upstream regulators, and the knockout silently undoes itself.
    """
    if n_propagation < 0:
        raise ValueError(f"n_propagation must be non-negative, got {n_propagation}")

    n_genes = delta.shape[1]
    if coefficient_matrix.shape != (n_genes, n_genes):
        raise ValueError(
            f"coefficient_matrix must be ({n_genes}, {n_genes}), "
            f"got {coefficient_matrix.shape}"
        )
    if clamped_values.shape != (delta.shape[0], len(clamped_indices)):
        raise ValueError(
            f"clamped_values must be ({delta.shape[0]}, {len(clamped_indices)}), "
            f"got {clamped_values.shape}"
        )

    current = delta.copy()
    for _ in range(n_propagation):
        current = current @ coefficient_matrix
        # Hold the perturbed genes fixed; see Notes.
        current[:, clamped_indices] = clamped_values
    return current


def simulate_knockout(
    expression: npt.NDArray[np.float64],
    coefficient_matrix: npt.NDArray[np.float64],
    gene_index: int,
    target_value: float = 0.0,
    n_propagation: int = 3,
) -> npt.NDArray[np.float64]:
    """Simulate knocking a single gene to a fixed value.

    Convenience wrapper around :func:`propagate_perturbation` for the common
    case of a single-gene knockout (``target_value=0.0``) or overexpression.

    Parameters
    ----------
    expression
        Unperturbed expression matrix, shape ``(n_cells, n_genes)``.
    coefficient_matrix
        Fitted GRN coefficients, shape ``(n_genes, n_genes)``.
    gene_index
        Column index of the gene to perturb.
    target_value
        Value to force the gene to. ``0.0`` is a knockout.
    n_propagation
        Number of propagation rounds.

    Returns
    -------
    numpy.ndarray
        Predicted perturbed expression, shape ``(n_cells, n_genes)``.

    Examples
    --------
    >>> import numpy as np
    >>> expression = np.array([[1.0, 2.0]])
    >>> coefficients = np.array([[0.0, 0.5], [0.0, 0.0]])
    >>> simulate_knockout(expression, coefficients, gene_index=0).shape
    (1, 2)
    """
    delta = np.zeros_like(expression)
    forced = np.full((expression.shape[0], 1), target_value)
    delta[:, [gene_index]] = forced - expression[:, [gene_index]]

    propagated = propagate_perturbation(
        delta=delta,
        coefficient_matrix=coefficient_matrix,
        clamped_indices=np.array([gene_index]),
        clamped_values=delta[:, [gene_index]],
        n_propagation=n_propagation,
    )
    return expression + propagated
