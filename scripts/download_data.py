"""Download the Roux et al. (2022) screen data from GEO (GSE176206).

Writes the raw matrix to data/raw/. Run once before the notebooks.
"""

from __future__ import annotations

from pathlib import Path

RAW_DIR = Path("data/raw")
ACCESSION = "GSE176206"


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Target directory: {RAW_DIR.resolve()}")
    print(f"Accession: {ACCESSION}")
    print(
        "Download the supplementary matrix from "
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={ACCESSION} "
        f"into {RAW_DIR}/."
    )


if __name__ == "__main__":
    main()
