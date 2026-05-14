#!/usr/bin/env python3
"""Export code cells from WMH_Analysis.ipynb into ordered step scripts.

This keeps the runnable scripts in sync with the original notebook without
manually copying code. The exporter is intentionally lightweight and only uses
Python's standard library.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

DEFAULT_NOTEBOOK = "WMH_Analysis.ipynb"
DEFAULT_OUTPUT_DIR = "scripts/notebook_steps"

STEP_SLUGS = {
    0: "flowchart",
    1: "remove_withdrawn_consent",
    2: "clean_final_dataset",
    3: "matching_ato_weighting",
    4: "wmh_transformation_diagnostics",
    5: "main_sensitivity_psm_wmh",
    6: "cognitive_secondary_outcomes",
    7: "wmh_cognitive_mediation",
    8: "mediation_plot_addon",
    9: "diagnostics_tables_5_cohorts",
    10: "additional_sensitivity_analysis",
    11: "stratify_cohorts",
    12: "stratified_ols_fdr",
    13: "interaction_ols_fdr",
    14: "annual_wmh_change",
    15: "duration_vs_wmh_loess",
    16: "print_package_versions",
}


def load_notebook(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def script_text(code_cell_index: int, source: str) -> str:
    title = source.strip().splitlines()[0] if source.strip() else f"Notebook code cell {code_cell_index}"
    header = (
        f"# Extracted from WMH_Analysis.ipynb, code cell {code_cell_index}.\n"
        f"# Notebook heading: {title}\n"
        "# Run this file from the repository root unless a local CONFIG section is edited.\n\n"
    )
    return header + source.rstrip() + "\n"


def export_steps(notebook_path: Path, output_dir: Path) -> list[Path]:
    notebook = load_notebook(notebook_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    code_cell_index = 0
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue

        source = "".join(cell.get("source", []))
        slug = STEP_SLUGS.get(code_cell_index, f"code_cell_{code_cell_index:02d}")
        path = output_dir / f"{code_cell_index + 1:02d}_{slug}.py"
        path.write_text(script_text(code_cell_index, source), encoding="utf-8")
        written.append(path)
        code_cell_index += 1

    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--notebook", default=DEFAULT_NOTEBOOK, help="Notebook to export")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Directory for exported scripts")
    args = parser.parse_args()

    notebook_path = Path(args.notebook)
    output_dir = Path(args.output_dir)
    written = export_steps(notebook_path, output_dir)

    for path in written:
        print(path)
    print(f"Exported {len(written)} code cells from {notebook_path} to {output_dir}")


if __name__ == "__main__":
    main()
