#!/usr/bin/env python3
# Copyright 2026 The Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Unified end-to-end Multiple Sequence Alignment and Visualization Pipeline.

Orchestrates pre-flight license verification, FASTA input validation, EBI
Clustal Omega alignment (or deterministic offline mock execution), Shannon
entropy and pairwise identity analysis, and standalone HTML5 dashboard
compilation.
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

try:
    from scripts.analyze_msa import analyze_alignment
    from scripts.generate_dashboard import (
        MANDATORY_CITATIONS,
        SCIENTIFIC_DISCLAIMER,
        compile_dashboard,
    )
    from scripts.msa_align import (
        MSAApiError,
        MSAValidationError,
        ensure_license_guard,
        run_alignment,
    )
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from analyze_msa import analyze_alignment  # type: ignore
    from generate_dashboard import (  # type: ignore
        MANDATORY_CITATIONS,
        SCIENTIFIC_DISCLAIMER,
        compile_dashboard,
    )
    from msa_align import (  # type: ignore
        MSAApiError,
        MSAValidationError,
        ensure_license_guard,
        run_alignment,
    )

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("run_pipeline")


def _format_rel_path(path: Path) -> str:
    """Format a path relative to the current working directory when possible."""
    try:
        return os.path.relpath(path, Path.cwd())
    except ValueError:
        return str(path)


def format_executive_summary(
    metrics: Dict[str, Any],
    dashboard_path: Path,
    aligned_fasta_path: Path,
    metrics_json_path: Path,
) -> str:
    """Build a structured plain-text executive summary with mandatory citations."""
    summary = metrics["summary"]
    columns = metrics["columns"]

    # Sort columns by conservation score descending, then position ascending
    top_conserved = sorted(
        columns,
        key=lambda c: (-float(c["conservation_score"]), int(c["position"])),
    )[:10]
    top_positions_str = ", ".join(
        f"Col {c['position']} ({c['modal_aa']}, score={c['conservation_score']:.2f})"
        for c in top_conserved
    )

    lines = [
        "========================================================================",
        "PROTEIN SEQUENCE MSA & VISUALIZATION PIPELINE - EXECUTIVE SUMMARY",
        "========================================================================",
        f"Number of Sequences         : {summary['num_sequences']}",
        f"Alignment Column Length (L) : {summary['alignment_length']}",
        f"Invariant Columns           : {summary['invariant_columns']} / {summary['alignment_length']} ({summary['conservation_index']:.2f}%)",
        f"Mean Conservation Score     : {summary['mean_conservation_score']:.4f}",
        f"Mean Pairwise Global ID     : {summary['mean_global_identity']:.2f}%",
        f"Mean Sequence Coverage ID   : {summary['mean_coverage_identity']:.2f}%",
        f"Mean Overlap Identity       : {summary['mean_overlap_identity']:.2f}%",
        f"Top Conserved Positions     : {top_positions_str}",
        "------------------------------------------------------------------------",
    ]

    mapped_features = metrics.get("features", {}).get("mapped_features", [])
    if mapped_features:
        lines.append("Functional Feature Projections:")
        for feat in mapped_features:
            if feat.get("status") == "mapped":
                lines.append(
                    f"  - {feat['name']}: Unaligned Pos {feat['unaligned_position']} "
                    f"({feat['reference_aa']}) -> Aligned Col {feat['alignment_column_position']} "
                    f"[{feat['conservation_state'].upper()}, conservation={feat['column_conservation']:.2f}]"
                )
            else:
                lines.append(
                    f"  - {feat.get('name')}: UNMAPPED ({feat.get('error', 'out of range')})"
                )
        lines.append(
            "------------------------------------------------------------------------"
        )

    lines.extend(
        [
            "Generated Artifacts:",
            f"  - Aligned FASTA  : {_format_rel_path(aligned_fasta_path)}",
            f"  - Metrics JSON   : {_format_rel_path(metrics_json_path)}",
            f"  - HTML Dashboard : {_format_rel_path(dashboard_path)}",
            "------------------------------------------------------------------------",
            "Academic Attribution (EMBL-EBI Clustal Omega):",
        ]
    )
    for cit in MANDATORY_CITATIONS:
        lines.append(
            f"  * {cit['label']}: {cit['authors']} {cit['title']} {cit['journal']} doi:{cit['doi']}"
        )
    lines.extend(
        [
            "------------------------------------------------------------------------",
            SCIENTIFIC_DISCLAIMER,
            "========================================================================",
        ]
    )
    return "\n".join(lines)


def execute_pipeline(
    input_fasta: Path,
    output_dir: Path,
    features: Optional[str] = None,
    title: Optional[str] = None,
    mock: bool = False,
    dry_run: bool = False,
) -> Dict[str, Any]:
    """Run the complete end-to-end alignment, analysis, and dashboard pipeline."""
    # Step 1: Pre-flight license guard
    license_file = ensure_license_guard()

    # If user passed a .html file path to -o instead of a directory, use its parent
    if output_dir.suffix.lower() in (".html", ".htm"):
        dashboard_path = output_dir
        target_dir = output_dir.parent
    else:
        target_dir = output_dir
        dashboard_path = target_dir / "dashboard.html"

    target_dir.mkdir(parents=True, exist_ok=True)
    aligned_fasta_path = target_dir / "aligned.fa"
    pim_path = target_dir / "matrix.pim"
    tree_path = target_dir / "tree.dnd"
    metrics_json_path = target_dir / "metrics.json"

    # Step 2 & 3: Validate input and execute alignment
    align_result = run_alignment(
        input_fasta=input_fasta,
        output_fasta=aligned_fasta_path,
        output_pim=pim_path,
        output_tree=tree_path,
        dry_run=dry_run,
        mock=mock,
    )
    if dry_run:
        return align_result

    # Step 4: Quantitative conservation and identity metrics
    metrics = analyze_alignment(
        alignment_path=aligned_fasta_path,
        features_arg=features,
        output_json_path=metrics_json_path,
    )

    # Step 5: Compile standalone interactive dashboard
    dashboard_title = title or f"MSA Analysis: {input_fasta.stem}"
    compile_dashboard(
        alignment_path=aligned_fasta_path,
        metrics_path=metrics_json_path,
        output_html_path=dashboard_path,
        pim_path=pim_path,
        tree_path=tree_path,
        title=dashboard_title,
    )

    # Step 6: Print formatted executive summary to stdout
    summary_text = format_executive_summary(
        metrics=metrics,
        dashboard_path=dashboard_path,
        aligned_fasta_path=aligned_fasta_path,
        metrics_json_path=metrics_json_path,
    )
    print(summary_text)

    return {
        "license_file": _format_rel_path(license_file),
        "aligned_fasta": _format_rel_path(aligned_fasta_path),
        "pim_file": _format_rel_path(pim_path),
        "tree_file": _format_rel_path(tree_path),
        "metrics_json": _format_rel_path(metrics_json_path),
        "dashboard_html": _format_rel_path(dashboard_path),
        "summary": metrics["summary"],
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Unified end-to-end protein Multiple Sequence Alignment (MSA), "
            "quantitative conservation profiling, and standalone HTML5 dashboard pipeline."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Path to input FASTA file containing >= 2 protein sequences.",
    )
    parser.add_argument(
        "-o",
        "--output",
        "--output-dir",
        dest="output_dir",
        type=Path,
        required=True,
        help="Output directory where aligned.fa, metrics.json, and dashboard.html will be saved.",
    )
    parser.add_argument(
        "--features",
        type=str,
        default=None,
        help="Optional JSON string or JSON file path mapping 1-based reference sequence features.",
    )
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Optional custom title for the generated HTML dashboard.",
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Run in deterministic offline mock mode without external network calls.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate input and print prepared EBI API payload without submitting.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        execute_pipeline(
            input_fasta=args.input,
            output_dir=args.output_dir,
            features=args.features,
            title=args.title,
            mock=args.mock,
            dry_run=args.dry_run,
        )
        return 0
    except (FileNotFoundError, MSAValidationError, MSAApiError) as exc:
        logger.error("%s", exc)
        return 1
    except Exception as exc:  # pragma: no cover
        logger.error("Unexpected error in run_pipeline: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
