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
"""Interactive multi-track sequence alignment visualizer module.

Provides biochemical color palettes (ClustalX, Zappo, Taylor, Hydrophobicity,
Monochromatic), SVG/HTML conservation bar track generators, and synchronized
dual-axis alignment viewer markup.
"""

from __future__ import annotations

import argparse
import html
import json
import logging
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("visualize_alignment")

COLOR_PALETTES: Dict[str, Dict[str, str]] = {
    "ClustalX": {
        "A": "#1F77B4",
        "I": "#1F77B4",
        "L": "#1F77B4",
        "M": "#1F77B4",
        "F": "#1F77B4",
        "W": "#1F77B4",
        "V": "#1F77B4",
        "N": "#2CA02C",
        "Q": "#2CA02C",
        "S": "#2CA02C",
        "T": "#2CA02C",
        "G": "#2CA02C",
        "D": "#E377C2",
        "E": "#E377C2",
        "K": "#D62728",
        "R": "#D62728",
        "C": "#FF7F0E",
        "P": "#FF7F0E",
        "H": "#17BECF",
        "Y": "#17BECF",
        "-": "#94A3B8",
        ".": "#94A3B8",
        "X": "#94A3B8",
    },
    "Zappo": {
        "I": "#FFA07A",
        "L": "#FFA07A",
        "V": "#FFA07A",
        "A": "#FFA07A",
        "M": "#FFA07A",
        "F": "#FFA07A",
        "W": "#FFA07A",
        "Y": "#FFA07A",
        "K": "#0000FF",
        "R": "#0000FF",
        "H": "#0000FF",
        "D": "#FF0000",
        "E": "#FF0000",
        "S": "#00FF00",
        "T": "#00FF00",
        "N": "#00FF00",
        "Q": "#00FF00",
        "P": "#FF00FF",
        "G": "#FF00FF",
        "C": "#FFFF00",
        "-": "#94A3B8",
        ".": "#94A3B8",
        "X": "#94A3B8",
    },
    "Taylor": {
        "A": "#CCFF00",
        "R": "#0000FF",
        "N": "#CC00FF",
        "D": "#FF0000",
        "C": "#FFFF00",
        "Q": "#FF00CC",
        "E": "#FF0066",
        "G": "#FF9900",
        "H": "#0066FF",
        "I": "#66FF00",
        "L": "#33FF00",
        "K": "#6600FF",
        "M": "#00FF00",
        "F": "#00FF66",
        "P": "#FFCC00",
        "S": "#FF3300",
        "T": "#FF6600",
        "W": "#00CCFF",
        "Y": "#00FFCC",
        "V": "#99FF00",
        "-": "#94A3B8",
        ".": "#94A3B8",
        "X": "#94A3B8",
    },
    "Hydrophobicity": {
        "I": "#D7191C",
        "V": "#DC2B24",
        "L": "#E3432F",
        "F": "#EF7546",
        "C": "#F48850",
        "M": "#F9AA60",
        "A": "#FBB265",
        "G": "#E6F598",
        "T": "#D1EC9C",
        "S": "#C7E89E",
        "W": "#BBE3A1",
        "Y": "#99D5A4",
        "P": "#80C7A5",
        "H": "#4393C3",
        "E": "#3682BA",
        "Q": "#3682BA",
        "D": "#3682BA",
        "N": "#3682BA",
        "K": "#2B6CB0",
        "R": "#1E4E8C",
        "-": "#94A3B8",
        ".": "#94A3B8",
        "X": "#94A3B8",
    },
    "Monochromatic": {
        aa: "#475569"
        for aa in "ACDEFGHIKLMNPQRSTVWY"
    }
    | {"-": "#94A3B8", ".": "#94A3B8", "X": "#94A3B8"},
}

KYTE_DOOLITTLE_SCORES: Dict[str, float] = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
    "-": 0.0,
}


def get_residue_color(residue: str, scheme: str = "ClustalX") -> str:
    """Return hex color code for an amino acid character under the requested scheme."""
    palette = COLOR_PALETTES.get(scheme, COLOR_PALETTES["ClustalX"])
    return palette.get(residue.upper(), "#94A3B8")


def get_conservation_bar_color(score: float) -> str:
    """Return color hex code for conservation bar height in [0.0, 1.0]."""
    if score >= 0.75:
        return "#16A34A"  # High conservation (green)
    if score >= 0.45:
        return "#D97706"  # Moderate conservation (amber/orange)
    return "#94A3B8"  # Variable (slate gray)


def extract_accession_and_label(seq_id: str) -> Dict[str, str]:
    """Parse UniProt-style 'sp|P69905|HBA_HUMAN' or generic FASTA identifiers."""
    parts = seq_id.split("|")
    if len(parts) >= 3:
        return {
            "db": parts[0],
            "accession": parts[1],
            "name": parts[2],
            "full_id": seq_id,
        }
    return {
        "db": "seq",
        "accession": seq_id[:12],
        "name": seq_id,
        "full_id": seq_id,
    }


def render_static_alignment_svg(
    sequences: List[Dict[str, Any]],
    columns: List[Dict[str, Any]],
    scheme: str = "ClustalX",
    cell_width: int = 16,
    cell_height: int = 22,
) -> str:
    """Generate a high-resolution standalone SVG of the alignment and conservation track."""
    label_width = 210
    top_margin = 76
    bottom_margin = 34
    align_len = len(columns)
    num_seqs = len(sequences)

    total_width = label_width + (align_len * cell_width) + 24
    total_height = top_margin + (num_seqs * cell_height) + bottom_margin

    svg_parts: List[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_width} {total_height}" '
        f'width="{total_width}" height="{total_height}" font-family="monospace" font-size="12">',
        f'<rect width="{total_width}" height="{total_height}" fill="#FFFFFF" />',
        '<text x="12" y="24" font-weight="bold" fill="#0F172A">Conservation</text>',
        '<text x="12" y="64" font-weight="bold" fill="#475569">Position</text>',
    ]

    # Render ruler and conservation bars
    for col in columns:
        idx = int(col["index"])
        pos = int(col["position"])
        score = float(col["conservation_score"])
        x = label_width + idx * cell_width
        bar_h = max(2, int(round(score * 34)))
        bar_y = 42 - bar_h
        bar_color = get_conservation_bar_color(score)
        svg_parts.append(
            f'<rect x="{x + 1}" y="{bar_y}" width="{cell_width - 2}" height="{bar_h}" '
            f'fill="{bar_color}" rx="1" />'
        )
        if pos == 1 or pos % 10 == 0:
            svg_parts.append(
                f'<text x="{x + 2}" y="64" font-size="10" fill="#475569">{pos}</text>'
            )

    # Render sequence rows
    for row_idx, seq_obj in enumerate(sequences):
        y = top_margin + row_idx * cell_height
        seq_id = html.escape(str(seq_obj["id"]))
        unaligned_len = int(seq_obj.get("unaligned_length", 0))
        svg_parts.append(
            f'<text x="12" y="{y + 15}" fill="#0F172A" font-weight="600">'
            f"{seq_id} ({unaligned_len} aa)</text>"
        )
        aligned_str = str(seq_obj["aligned_sequence"])
        for col_idx, aa in enumerate(aligned_str):
            x = label_width + col_idx * cell_width
            if aa not in ("-", "."):
                bg = get_residue_color(aa, scheme=scheme)
                svg_parts.append(
                    f'<rect x="{x}" y="{y + 1}" width="{cell_width - 1}" '
                    f'height="{cell_height - 2}" fill="{bg}" fill-opacity="0.28" rx="2" />'
                )
            text_color = "#0F172A" if aa not in ("-", ".") else "#94A3B8"
            svg_parts.append(
                f'<text x="{x + (cell_width // 2)}" y="{y + 15}" text-anchor="middle" '
                f'fill="{text_color}">{html.escape(aa)}</text>'
            )

    # Consensus row
    cons_y = top_margin + num_seqs * cell_height + 6
    svg_parts.append(
        f'<text x="12" y="{cons_y + 14}" font-weight="bold" fill="#0F172A">Consensus</text>'
    )
    for col in columns:
        idx = int(col["index"])
        cons_aa = html.escape(str(col["consensus_aa"]))
        x = label_width + idx * cell_width
        svg_parts.append(
            f'<text x="{x + (cell_width // 2)}" y="{cons_y + 14}" text-anchor="middle" '
            f'font-weight="bold" fill="#1E293B">{cons_aa}</text>'
        )

    svg_parts.append("</svg>")
    return "\n".join(svg_parts)


def generate_alignment_component_html() -> str:
    """Return the HTML5 container markup for the interactive alignment viewer."""
    return """
<section id="alignment-viewer-section" class="panel-card">
  <div class="panel-header">
    <div>
      <h2 class="panel-title">Interactive Multi-Track Sequence Alignment Viewer</h2>
      <p class="panel-subtitle">Synchronized dual-axis alignment matrix with per-column Shannon entropy conservation bars and dynamic consensus sequence</p>
    </div>
    <div class="toolbar-controls">
      <label class="control-label" for="color-scheme-select">Color Scheme
        <select id="color-scheme-select" class="control-input">
          <option value="ClustalX" selected>ClustalX (Biochemical)</option>
          <option value="Zappo">Zappo (Physicochemical)</option>
          <option value="Taylor">Taylor (Rainbow)</option>
          <option value="Hydrophobicity">Hydrophobicity (Kyte-Doolittle)</option>
          <option value="Monochromatic">Monochromatic</option>
        </select>
      </label>
      <label class="control-label" for="consensus-threshold-slider">Consensus Threshold: <span id="consensus-threshold-val">80%</span>
        <input id="consensus-threshold-slider" type="range" min="50" max="100" step="5" value="80" class="control-slider" />
      </label>
      <label class="control-label" for="font-size-slider">Cell Size: <span id="font-size-val">18px</span>
        <input id="font-size-slider" type="range" min="14" max="26" step="2" value="18" class="control-slider" />
      </label>
      <label class="control-label" for="jump-position-input">Jump to Position
        <input id="jump-position-input" type="number" min="1" placeholder="Col #" class="control-input control-input-sm" />
      </label>
      <label class="control-label" for="motif-search-input">Find Motif (Regex)
        <input id="motif-search-input" type="text" placeholder="e.g. H..H or G..G.G" class="control-input" />
      </label>
      <button id="clear-highlights-btn" type="button" class="action-btn action-btn-secondary">Clear Filter</button>
    </div>
  </div>

  <div id="motif-match-banner" class="info-banner hidden"></div>

  <div class="msa-viewport-wrapper">
    <div id="msa-fixed-headers" class="msa-fixed-headers">
      <div class="msa-header-cell msa-ruler-label">Position Ruler</div>
      <div class="msa-header-cell msa-conservation-label">Conservation Score</div>
      <div id="msa-sequence-labels" class="msa-sequence-labels"></div>
      <div class="msa-header-cell msa-consensus-label">Consensus Sequence</div>
    </div>
    <div id="msa-scroll-viewport" class="msa-scroll-viewport">
      <div id="msa-ruler-track" class="msa-ruler-track"></div>
      <div id="msa-conservation-track" class="msa-conservation-track"></div>
      <div id="msa-sequence-tracks" class="msa-sequence-tracks"></div>
      <div id="msa-consensus-track" class="msa-consensus-track"></div>
    </div>
  </div>

  <div id="column-inspector-drawer" class="inspector-drawer">
    <div class="inspector-header">
      <h3 id="inspector-title" class="inspector-title">Column Inspector (Click any alignment column to inspect residue breakdown)</h3>
      <span id="inspector-badge" class="status-badge">Ready</span>
    </div>
    <div id="inspector-body" class="inspector-body">
      <p class="text-muted">Select a column in the alignment viewer or conservation bar chart above to inspect exact Shannon entropy (bits), normalized conservation index, gap frequency, and amino acid distribution.</p>
    </div>
  </div>
</section>
""".strip()


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate standalone alignment SVG or HTML fragment from MSA metrics."
    )
    parser.add_argument(
        "--metrics",
        type=Path,
        required=True,
        help="Path to computed metrics.json from analyze_msa.py.",
    )
    parser.add_argument(
        "--scheme",
        type=str,
        default="ClustalX",
        choices=list(COLOR_PALETTES.keys()),
        help="Biochemical color scheme (default: ClustalX).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output file path (.svg or .html).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if not args.metrics.is_file():
        logger.error("Metrics JSON file not found: %s", args.metrics)
        return 1

    metrics_data = json.loads(args.metrics.read_text(encoding="utf-8"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.suffix.lower() == ".svg":
        svg_text = render_static_alignment_svg(
            sequences=metrics_data["sequences"],
            columns=metrics_data["columns"],
            scheme=args.scheme,
        )
        args.output.write_text(svg_text, encoding="utf-8")
    else:
        html_fragment = generate_alignment_component_html()
        args.output.write_text(html_fragment + "\n", encoding="utf-8")

    logger.info("Wrote alignment visualization artifact to %s", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
