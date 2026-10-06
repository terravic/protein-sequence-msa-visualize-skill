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
"""Command-line smoke tests verifying every CLI script entrypoint."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _run_cli(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-B", *args],
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )


def test_msa_align_dry_run_smoke() -> None:
    proc = _run_cli(
        [
            "scripts/msa_align.py",
            "sample_data/sample_globins.fasta",
            "--dry-run",
        ]
    )
    assert proc.returncode == 0, proc.stderr
    payload = json.loads(proc.stdout)
    assert payload["num_sequences"] == 5
    assert "sp|P69905|HBA_HUMAN" in payload["sequence_ids"]


def test_msa_align_mock_smoke(tmp_path: Path) -> None:
    out_fa = tmp_path / "aligned.fa"
    out_pim = tmp_path / "matrix.pim"
    out_tree = tmp_path / "tree.dnd"
    proc = _run_cli(
        [
            "scripts/msa_align.py",
            "sample_data/sample_globins.fasta",
            "-o",
            str(out_fa),
            "--output-pim",
            str(out_pim),
            "--output-tree",
            str(out_tree),
            "--mock",
        ]
    )
    assert proc.returncode == 0, proc.stderr
    assert out_fa.is_file() and out_fa.stat().st_size > 0
    assert out_pim.is_file() and out_pim.stat().st_size > 0
    assert out_tree.is_file() and out_tree.stat().st_size > 0


def test_analyze_msa_smoke(tmp_path: Path) -> None:
    out_json = tmp_path / "metrics.json"
    proc = _run_cli(
        [
            "scripts/analyze_msa.py",
            "--alignment",
            "sample_data/sample_alignment_result.fa",
            "--features",
            '{"reference_id": "sp|P69905|HBA_HUMAN", "features": [{"name": "Distal His", "position": 59}]}',
            "--output-json",
            str(out_json),
        ]
    )
    assert proc.returncode == 0, proc.stderr
    assert out_json.is_file()
    metrics = json.loads(out_json.read_text(encoding="utf-8"))
    assert metrics["summary"]["num_sequences"] == 5
    assert metrics["summary"]["alignment_length"] == 165
    assert metrics["features"]["mapped_features"][0]["status"] == "mapped"


def test_visualize_modules_smoke(tmp_path: Path) -> None:
    metrics_json = tmp_path / "metrics.json"
    _run_cli(
        [
            "scripts/analyze_msa.py",
            "--alignment",
            "sample_data/sample_alignment_result.fa",
            "--output-json",
            str(metrics_json),
        ]
    )
    aln_svg = tmp_path / "alignment.svg"
    proc_aln = _run_cli(
        [
            "scripts/visualize_alignment.py",
            "--metrics",
            str(metrics_json),
            "--scheme",
            "ClustalX",
            "-o",
            str(aln_svg),
        ]
    )
    assert proc_aln.returncode == 0, proc_aln.stderr
    assert aln_svg.is_file() and "<svg" in aln_svg.read_text(encoding="utf-8")

    tree_svg = tmp_path / "tree.svg"
    proc_tree = _run_cli(
        [
            "scripts/visualize_matrix.py",
            "--tree",
            "sample_data/sample_phylotree.dnd",
            "--pim",
            "sample_data/sample_pim.json",
            "-o",
            str(tree_svg),
        ]
    )
    assert proc_tree.returncode == 0, proc_tree.stderr
    assert tree_svg.is_file() and "<svg" in tree_svg.read_text(encoding="utf-8")


def test_generate_dashboard_smoke(tmp_path: Path) -> None:
    metrics_json = tmp_path / "metrics.json"
    _run_cli(
        [
            "scripts/analyze_msa.py",
            "--alignment",
            "sample_data/sample_alignment_result.fa",
            "--output-json",
            str(metrics_json),
        ]
    )
    dashboard_html = tmp_path / "dashboard.html"
    proc = _run_cli(
        [
            "scripts/generate_dashboard.py",
            "--alignment",
            "sample_data/sample_alignment_result.fa",
            "--metrics",
            str(metrics_json),
            "--pim",
            "sample_data/sample_pim.json",
            "--tree",
            "sample_data/sample_phylotree.dnd",
            "-o",
            str(dashboard_html),
        ]
    )
    assert proc.returncode == 0, proc.stderr
    assert dashboard_html.is_file()
    html_text = dashboard_html.read_text(encoding="utf-8")
    assert 'id="msa-data"' in html_text
    assert 'id="homology-3d-viewport"' in html_text
    assert "Sievers, F. et al." in html_text
    assert "Madeira, F. et al." in html_text
