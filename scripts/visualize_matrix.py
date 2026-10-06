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
"""Pairwise Percent Identity Matrix (PIM) heatmap and Newick dendrogram generator.

Parses Newick phylogenetic guide trees and Clustal Omega PIM files, computes
hierarchical cladogram SVG coordinates, and provides HTML5/SVG components for
interactive pairwise similarity exploration.
"""

from __future__ import annotations

import argparse
import html
import json
import logging
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("visualize_matrix")


def parse_pim_file(pim_path: Optional[Path]) -> Optional[Dict[str, Any]]:
    """Parse a Percent Identity Matrix from JSON or Clustal Omega text format."""
    if pim_path is None or not pim_path.is_file():
        return None
    raw = pim_path.read_text(encoding="utf-8").strip()
    if not raw:
        return None

    if raw.startswith("{"):
        try:
            data = json.loads(raw)
            if "matrix" in data and "sequence_ids" in data:
                return data
        except json.JSONDecodeError:
            pass

    # Parse Clustal text PIM format:
    # 1: sp|P69905|HBA_HUMAN 100.00 43.66 ...
    seq_ids: List[str] = []
    matrix: List[List[float]] = []
    for line in raw.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        cleaned = re.sub(r"^\d+\s*:\s*", "", s)
        tokens = cleaned.split()
        if len(tokens) < 2:
            continue
        seq_id = tokens[0]
        row_vals: List[float] = []
        for tok in tokens[1:]:
            try:
                row_vals.append(round(float(tok), 2))
            except ValueError:
                continue
        if row_vals:
            seq_ids.append(seq_id)
            matrix.append(row_vals)

    if seq_ids and matrix:
        return {
            "description": "Parsed Clustal Omega Percent Identity Matrix",
            "sequence_ids": seq_ids,
            "matrix": matrix,
            "raw_pim_text": raw,
        }
    return None


def parse_newick(newick_str: str) -> Dict[str, Any]:
    """Parse a Newick tree string into a nested dictionary hierarchy.

    Returns a dictionary node:
        {"name": str, "length": float, "children": List[node]}
    """
    cleaned = "".join(newick_str.strip().split())
    if cleaned.endswith(";"):
        cleaned = cleaned[:-1]
    if not cleaned:
        return {"name": "root", "length": 0.0, "children": []}

    pos = 0
    n = len(cleaned)

    def _parse_subtree() -> Dict[str, Any]:
        nonlocal pos
        children: List[Dict[str, Any]] = []
        if pos < n and cleaned[pos] == "(":
            pos += 1  # consume '('
            while True:
                children.append(_parse_subtree())
                if pos < n and cleaned[pos] == ",":
                    pos += 1
                elif pos < n and cleaned[pos] == ")":
                    pos += 1
                    break
                else:
                    break

        label_chars: List[str] = []
        while pos < n and cleaned[pos] not in (",", ")", ":", ";"):
            label_chars.append(cleaned[pos])
            pos += 1
        name = "".join(label_chars).strip()

        branch_len = 0.0
        if pos < n and cleaned[pos] == ":":
            pos += 1
            len_chars: List[str] = []
            while pos < n and cleaned[pos] not in (",", ")", ";"):
                len_chars.append(cleaned[pos])
                pos += 1
            try:
                branch_len = max(0.0, float("".join(len_chars)))
            except ValueError:
                branch_len = 0.0

        return {
            "name": name,
            "length": round(branch_len, 5),
            "children": children,
        }

    tree = _parse_subtree()
    if not tree.get("name") and tree.get("children"):
        tree["name"] = "root"
    return tree


def compute_tree_layout(root: Dict[str, Any]) -> Dict[str, Any]:
    """Assign cumulative X distances and leaf Y indices to a parsed Newick tree."""
    leaves: List[str] = []
    leaf_counter = 0

    def _assign(node: Dict[str, Any], parent_dist: float, is_root: bool = False) -> Tuple[float, float]:
        nonlocal leaf_counter
        own_len = 0.0 if is_root else float(node.get("length", 0.0))
        effective_len = own_len if own_len > 0.0 else (0.0 if is_root else 0.05)
        cum_x = round(parent_dist + effective_len, 5)
        node["x_dist"] = cum_x

        children = node.get("children", [])
        if not children:
            y_pos = float(leaf_counter)
            leaf_counter += 1
            node["y_pos"] = y_pos
            node["is_leaf"] = True
            leaves.append(str(node.get("name", f"Leaf_{leaf_counter}")))
            return y_pos, cum_x

        child_ys: List[float] = []
        max_x = cum_x
        for child in children:
            cy, cx = _assign(child, cum_x, is_root=False)
            child_ys.append(cy)
            if cx > max_x:
                max_x = cx

        y_pos = sum(child_ys) / len(child_ys)
        node["y_pos"] = round(y_pos, 4)
        node["is_leaf"] = False
        return y_pos, max_x

    _, max_dist = _assign(root, 0.0, is_root=True)
    return {
        "root": root,
        "leaves": leaves,
        "num_leaves": max(1, leaf_counter),
        "max_distance": max(0.01, round(max_dist, 5)),
    }


def render_tree_svg(
    newick_str: str,
    width: int = 640,
    row_height: int = 46,
) -> str:
    """Render an interactive SVG rectangular phylogram from a Newick string."""
    parsed = parse_newick(newick_str)
    layout = compute_tree_layout(parsed)
    num_leaves = int(layout["num_leaves"])
    max_dist = float(layout["max_distance"])

    top_pad = 28
    bottom_pad = 28
    left_pad = 28
    right_label_pad = 240
    height = max(180, top_pad + bottom_pad + max(1, num_leaves - 1) * row_height + 24)
    plot_width = max(160, width - left_pad - right_label_pad)

    def _scale_x(dist: float) -> float:
        return round(left_pad + (dist / max_dist) * plot_width, 2)

    def _scale_y(y_idx: float) -> float:
        if num_leaves <= 1:
            return round(height / 2.0, 2)
        return round(top_pad + y_idx * row_height, 2)

    elements: List[str] = [
        f'<svg id="phylotree-svg" xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 0 {width} {height}" width="100%" height="{height}" '
        f'class="phylotree-svg">',
    ]

    def _draw_node(node: Dict[str, Any]) -> None:
        nx = _scale_x(float(node["x_dist"]))
        ny = _scale_y(float(node["y_pos"]))
        children = node.get("children", [])

        for child in children:
            cx = _scale_x(float(child["x_dist"]))
            cy = _scale_y(float(child["y_pos"]))
            # Vertical connector from parent Y to child Y at parent X
            elements.append(
                f'<line x1="{nx}" y1="{ny}" x2="{nx}" y2="{cy}" '
                f'stroke="#64748B" stroke-width="2" />'
            )
            # Horizontal branch from parent X to child X at child Y
            elements.append(
                f'<line x1="{nx}" y1="{cy}" x2="{cx}" y2="{cy}" '
                f'stroke="#64748B" stroke-width="2" />'
            )
            branch_len = float(child.get("length", 0.0))
            if branch_len > 0:
                mid_x = round((nx + cx) / 2.0, 2)
                elements.append(
                    f'<text x="{mid_x}" y="{cy - 5}" text-anchor="middle" '
                    f'font-size="10" fill="#64748B" font-family="monospace">'
                    f"{branch_len:.3f}</text>"
                )
            _draw_node(child)

        if node.get("is_leaf"):
            leaf_name = str(node.get("name", "Leaf"))
            escaped_name = html.escape(leaf_name)
            elements.append(
                f'<g class="tree-leaf-node" data-seq-id="{escaped_name}" '
                f' style="cursor:pointer;">'
                f'<circle cx="{nx}" cy="{ny}" r="5" fill="#0284C7" stroke="#FFFFFF" stroke-width="1.5" />'
                f'<text x="{nx + 10}" y="{ny + 4}" font-size="12" font-weight="600" '
                f'font-family="monospace" fill="currentColor">{escaped_name}</text>'
                f"</g>"
            )
        else:
            elements.append(
                f'<circle cx="{nx}" cy="{ny}" r="3.5" fill="#475569" />'
            )

    _draw_node(layout["root"])
    elements.append("</svg>")
    return "\n".join(elements)


def get_identity_cell_color(pct: float) -> str:
    """Return hex color for a pairwise identity percentage (0..100%)."""
    clamped = max(0.0, min(100.0, pct))
    if clamped >= 90.0:
        return "#059669"  # High identity emerald
    if clamped >= 70.0:
        return "#10B981"
    if clamped >= 50.0:
        return "#0D9488"  # Moderate-high teal
    if clamped >= 35.0:
        return "#0284C7"  # Moderate sky blue
    if clamped >= 20.0:
        return "#2563EB"  # Low-moderate royal blue
    return "#1E3A8A"  # Distant homology deep blue


def generate_matrix_and_tree_component_html(tree_svg: str) -> str:
    """Return the HTML5 container markup for the PIM heatmap and phylogenetic cladogram."""
    return f"""
<div class="dual-grid-section">
  <section id="pim-heatmap-section" class="panel-card">
    <div class="panel-header">
      <div>
        <h2 class="panel-title">All-vs-All Pairwise Percent Identity Matrix (PIM)</h2>
        <p class="panel-subtitle">Interactive homology matrix across four mathematical denominator definitions</p>
      </div>
      <div class="toolbar-controls">
        <label class="control-label" for="pim-metric-select">Identity Metric
          <select id="pim-metric-select" class="control-input">
            <option value="global" selected>Global Identity (Matches / Alignment Length)</option>
            <option value="coverage">Sequence Coverage Identity (Matches / Min Length)</option>
            <option value="overlap">Overlap Identity (Excludes Terminal Gaps)</option>
            <option value="external">Clustal Omega PIM Endpoint</option>
          </select>
        </label>
        <button id="export-pim-csv-btn" type="button" class="action-btn action-btn-secondary">Export CSV</button>
        <button id="export-pim-json-btn" type="button" class="action-btn action-btn-secondary">Export JSON</button>
      </div>
    </div>

    <div class="pim-summary-strip">
      <div class="metric-pill">
        <span class="metric-pill-label">Multisequence Conservation Index:</span>
        <strong id="pim-conservation-index-val">--</strong>
      </div>
      <div class="metric-pill">
        <span class="metric-pill-label">Active Metric Mean Pairwise Identity:</span>
        <strong id="pim-active-mean-val">--</strong>
      </div>
    </div>

    <div id="pim-heatmap-container" class="pim-heatmap-container"></div>
    <div id="pim-tooltip" class="pim-tooltip hidden"></div>
  </section>

  <section id="phylotree-section" class="panel-card">
    <div class="panel-header">
      <div>
        <h2 class="panel-title">Phylogenetic Guide Tree Cladogram</h2>
        <p class="panel-subtitle">Hierarchical clustering tree (Newick format) with evolutionary branch lengths; hover or click leaf nodes to highlight sequence tracks</p>
      </div>
    </div>
    <div id="phylotree-container" class="phylotree-container">
      {tree_svg}
    </div>
    <details class="newick-raw-details">
      <summary>View Raw Newick Tree String</summary>
      <pre id="newick-raw-pre" class="newick-raw-pre"></pre>
    </details>
  </section>
</div>
""".strip()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate pairwise identity heatmap and Newick phylogenetic tree SVG/HTML."
    )
    parser.add_argument(
        "--tree",
        type=Path,
        required=True,
        help="Path to Newick guide tree (.dnd) file.",
    )
    parser.add_argument(
        "--pim",
        type=Path,
        default=None,
        help="Optional path to Clustal PIM or JSON file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output path (.svg or .html).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if not args.tree.is_file():
        logger.error("Newick tree file not found: %s", args.tree)
        return 1

    newick_str = args.tree.read_text(encoding="utf-8")
    tree_svg = render_tree_svg(newick_str)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.output.suffix.lower() == ".svg":
        args.output.write_text(tree_svg + "\n", encoding="utf-8")
    else:
        html_content = generate_matrix_and_tree_component_html(tree_svg)
        args.output.write_text(html_content + "\n", encoding="utf-8")

    logger.info("Saved matrix/tree visualization to %s", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
