"""Unit tests for GRN signal propagation.

Fast and data-free so they run in CI on every push. The check against
CellOracle's output on reference data is marked ``slow`` and skipped by default.
"""

from __future__ import annotations

import numpy as np
import pytest

from roux_reprogramming.propagation import propagate_perturbation, simulate_knockout


@pytest.fixture
def simple_network() -> np.ndarray:
    """A three-gene chain: gene0 -> gene1 -> gene2, activating.

    Chosen so propagation depth is directly observable: gene2 can only be
    reached from gene0 in two rounds.
    """
    coefficients = np.zeros((3, 3))
    coefficients[0, 1] = 0.5  # gene0 activates gene1
    coefficients[1, 2] = 0.5  # gene1 activates gene2
    return coefficients


def test_knockout_sets_target_gene_to_zero(simple_network: np.ndarray) -> None:
    """A knockout must leave the perturbed gene at its forced value."""
    expression = np.array([[4.0, 4.0, 4.0]])
    result = simulate_knockout(expression, simple_network, gene_index=0)
    assert result[0, 0] == pytest.approx(0.0)


def test_propagation_reaches_downstream_genes(simple_network: np.ndarray) -> None:
    """Knocking out gene0 must reduce gene1, which activates it."""
    expression = np.array([[4.0, 4.0, 4.0]])
    result = simulate_knockout(expression, simple_network, gene_index=0)
    assert result[0, 1] < expression[0, 1]


def test_propagation_depth_is_respected(simple_network: np.ndarray) -> None:
    """With one round, influence cannot yet have reached the two-step target.

    This is the test that would catch an off-by-one in the propagation loop.
    """
    expression = np.array([[4.0, 4.0, 4.0]])
    one_round = simulate_knockout(
        expression, simple_network, gene_index=0, n_propagation=1
    )
    two_rounds = simulate_knockout(
        expression, simple_network, gene_index=0, n_propagation=2
    )
    assert one_round[0, 2] == pytest.approx(expression[0, 2])
    assert two_rounds[0, 2] < expression[0, 2]


def test_repression_increases_target(simple_network: np.ndarray) -> None:
    """Removing a repressor must raise its target (sign handling)."""
    coefficients = simple_network.copy()
    coefficients[0, 1] = -0.5  # gene0 now represses gene1
    expression = np.array([[4.0, 4.0, 4.0]])
    result = simulate_knockout(expression, coefficients, gene_index=0)
    assert result[0, 1] > expression[0, 1]


def test_zero_propagation_leaves_network_untouched(simple_network: np.ndarray) -> None:
    """With no propagation rounds, only the clamped gene changes."""
    expression = np.array([[4.0, 4.0, 4.0]])
    result = simulate_knockout(
        expression, simple_network, gene_index=0, n_propagation=0
    )
    np.testing.assert_allclose(result[0, 1:], expression[0, 1:])


def test_rejects_negative_propagation(simple_network: np.ndarray) -> None:
    expression = np.array([[4.0, 4.0, 4.0]])
    with pytest.raises(ValueError, match="non-negative"):
        simulate_knockout(expression, simple_network, gene_index=0, n_propagation=-1)


def test_rejects_mismatched_coefficient_matrix() -> None:
    with pytest.raises(ValueError, match="coefficient_matrix"):
        propagate_perturbation(
            delta=np.zeros((2, 3)),
            coefficient_matrix=np.zeros((4, 4)),
            clamped_indices=np.array([0]),
            clamped_values=np.zeros((2, 1)),
        )


@pytest.mark.slow
def test_matches_celloracle_on_tutorial_data() -> None:
    """Verify this implementation against CellOracle's own simulation.

    Skipped by default: requires cached reference data and is slow. Run with
    ``pytest -m slow``.
    """
    pytest.skip("TODO: implement once a reference Oracle object is cached locally.")
