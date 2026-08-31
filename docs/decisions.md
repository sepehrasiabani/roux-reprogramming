# Decision log

Design decisions taken during the analysis, with the alternatives considered and
the reasoning.

---

## No trajectory / pseudotime analysis
**Decision:** Treat cell types as parallel identities and use group contrasts (young vs. aged, treated vs. untreated); do not infer a developmental trajectory or pseudotime.
**Alternatives:** Trajectory inference (diffusion map, PAGA, DPT) with developmental-flow scoring.
**Reasoning:** The populations are distinct mesenchymal identities, not stages along one differentiation path. On the embedding, treated combinations are intermixed and age only partially separates; there is no continuous developmental axis to order cells along.

## Marker-gene annotation
**Decision:** Annotate Leiden clusters by canonical marker genes validated with `rank_genes_groups`.
**Alternatives:** Reference-based label transfer, or no annotation.
**Reasoning:** Cluster identity determines every per-cell-type result and must be correct. Marker validation refined two initial assignments. The dataset is pure mesenchyme, so a generic reference would add little.

## Continuous ageing score rather than the age label
**Decision:** Build a continuous per-cell ageing score from a gene signature rather than using the binary young/aged label.
**Alternatives:** Two-group contrasts on the binary label only.
**Reasoning:** Gradient construction requires a continuous quantity encoding direction and magnitude, which a binary label cannot provide; a continuous score also captures within-group heterogeneity in ageing.

## Data-derived ageing signature
**Decision:** Use a data-derived aged-vs-young signature (a myofibroblast/activation program) as the ageing axis.
**Alternatives:** A published SASP/senescence panel.
**Reasoning:** The published SASP panel anti-correlated with age in this dataset (r = −0.32); SASP genes are minimally expressed and the senescence core is mixed. The data-derived signature reflects the ageing program actually present in these cells.

## Retain the Yamanaka factors in the feature set
**Decision:** Retain Sox2, Oct4/Pou5f1, Klf4, and Myc in the HVG feature set regardless of variance ranking, and rebuild the GRN accordingly.
**Alternatives:** Use the default top-2000 HVG, which excluded Myc and Oct4.
**Reasoning:** The analysis concerns factor combinations; excluding two of the four factors would make the core question unanswerable. Genes required as perturbation targets are retained irrespective of variance.

## Overexpression capped at ~2× observed maximum; dose-response control
**Decision:** Set each factor's simulated overexpression to ~2× its observed maximum, and test 2×, 1×, and 0.5× as a dose-response.
**Alternatives:** A single fixed high value.
**Reasoning:** CellOracle rejects perturbations beyond ~2× the observed maximum. The dose-response distinguishes a directional effect from an artefact of extreme values; the pro-ageing direction held across doses.

## Predicted-vs-measured evaluated against null models
**Decision:** Evaluate the predicted-vs-measured correlation against a gene-shuffled null and a combination-mismatched null.
**Alternatives:** Report the raw correlation alone.
**Reasoning:** A correlation is uninterpretable without a null. The gene-shuffled null tests whether the correct genes are predicted; the combination-mismatched null tests whether predictions are combination-specific.

## Activation-independent ageing axis
**Decision:** Construct an ageing axis from ageing genes excluding myofibroblast/contractile markers, and re-score.
**Alternatives:** Report the entanglement as an unresolved caveat.
**Reasoning:** The pro-ageing result is interpretable only if "aged" is separable from "activated." The activation-independent axis still separated young from aged, correlated r = 0.74 with the original axis, and gave the same perturbation score — indicating the two processes are inseparable in these cells.

## Combination classifier with a cell-type control
**Decision:** Test combination separability with a supervised classifier paired with a cell-type classifier as a control.
**Alternatives:** Rely on the embedding and the combination-mismatched null alone.
**Reasoning:** A negative classifier result is credible only if the same pipeline succeeds on a task with real signal. Cell type classified at 86% and combination at 11%, isolating the failure to the combinations.

## Repository layout
**Decision:** A `src/` package with thin notebooks; CI running linting, type-checking, and tests.
**Alternatives:** Notebook-only.
**Reasoning:** Logic in `src/` is testable and importable; notebook-only logic is not.
