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
"""End-to-end system, compliance, and synthetic data verification tests."""

from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
import sys

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SKILL_SUBDIR = PROJECT_ROOT / "skills" / "protein_sequence_msa_visualize"


def _run_cli(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-B", *args],
        cwd=str(PROJECT_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )


def test_end_to_end_pipeline_globins_with_features(tmp_path: Path) -> None:
    out_dir = tmp_path / "globin_e2e"
    features_arg = json.dumps(
        {
            "reference_id": "sp|P69905|HBA_HUMAN",
            "features": [
                {"name": "Distal Histidine (E7)", "position": 59, "type": "active_site"},
                {"name": "Proximal Histidine (F8)", "position": 88, "type": "binding_site"},
            ],
        }
    )
    proc = _run_cli(
        [
            "scripts/run_pipeline.py",
            "sample_data/sample_globins.fasta",
            "--features",
            features_arg,
            "-o",
            str(out_dir),
            "--mock",
        ]
    )
    assert proc.returncode == 0, proc.stderr
    assert "EXECUTIVE SUMMARY" in proc.stdout
    assert "Distal Histidine (E7): Unaligned Pos 59 (H) -> Aligned Col 66 [INVARIANT" in proc.stdout
    assert "Proximal Histidine (F8): Unaligned Pos 88 (H) -> Aligned Col 95 [INVARIANT" in proc.stdout

    for filename in ("aligned.fa", "matrix.pim", "tree.dnd", "metrics.json", "dashboard.html"):
        artifact = out_dir / filename
        assert artifact.is_file() and artifact.stat().st_size > 0

    html_content = (out_dir / "dashboard.html").read_text(encoding="utf-8")
    assert "#1F77B4" in html_content  # ClustalX hydrophobic blue
    assert "10.1038/msb.2011.75" in html_content
    assert "10.1093/nar/gkae241" in html_content
    assert "10.1093/nar/gkq313" in html_content
    assert "AS-IS Scientific Disclaimer" in html_content
    assert "input-data-profile-section" in html_content
    assert "Input Sequences &amp; Residue Composition" in html_content
    assert "homology-3d-viewport" in html_content
    assert "3D Sequence Homology Network" in html_content


def test_root_dashboard_and_ui_requirements() -> None:
    root_dashboard = PROJECT_ROOT / "dashboard.html"
    assert root_dashboard.is_file() and root_dashboard.stat().st_size > 0
    html_content = root_dashboard.read_text(encoding="utf-8")

    btn_match = re.search(
        r'<button[^>]*id="toggle-theme-btn"[^>]*>(.*?)</button>',
        html_content,
        flags=re.DOTALL,
    )
    assert btn_match is not None, "Missing #toggle-theme-btn in dashboard.html"
    btn_inner = btn_match.group(1)
    assert "theme-icon-moon" in btn_inner
    assert "theme-icon-sun" in btn_inner
    stripped_text = re.sub(r"<[^>]+>", "", btn_inner).strip()
    assert stripped_text == ""

    assert 'id="input-data-profile-section"' in html_content
    assert 'id="homology-3d-viewport"' in html_content
    assert "3D Sequence Homology Network" in html_content
    assert "3D Conservation Landscape" in html_content


def test_end_to_end_pipeline_ubiquitin_and_skill_subdir_parity(tmp_path: Path) -> None:
    out_dir = tmp_path / "skill_subdir_e2e"
    proc = _run_cli(
        [
            "skills/protein_sequence_msa_visualize/scripts/run_pipeline.py",
            "skills/protein_sequence_msa_visualize/sample_data/sample_ubiquitin.fasta",
            "-o",
            str(out_dir),
            "--mock",
        ]
    )
    assert proc.returncode == 0, proc.stderr
    metrics = json.loads((out_dir / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["summary"]["num_sequences"] == 4
    assert metrics["summary"]["alignment_length"] >= 76
    assert (out_dir / "dashboard.html").is_file()


def test_pipeline_error_handling_on_invalid_input(tmp_path: Path) -> None:
    bad_fasta = tmp_path / "only_one.fasta"
    bad_fasta.write_text(">seq1\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH\n", encoding="utf-8")
    proc = _run_cli(
        [
            "scripts/run_pipeline.py",
            str(bad_fasta),
            "-o",
            str(tmp_path / "fail_out"),
            "--mock",
        ]
    )
    assert proc.returncode != 0


def test_legal_and_documentation_compliance() -> None:
    license_guard = PROJECT_ROOT / ".licenses" / "protein_sequence_msa_visualize_LICENSE.txt"
    assert license_guard.is_file()
    skill_licenses = (PROJECT_ROOT / "SKILL_LICENSES.md").read_text(encoding="utf-8")
    assert "protein-sequence-msa-visualize-skill | https://www.ebi.ac.uk/about/terms-of-use/ and https://www.ebi.ac.uk/jdispatcher/msa/clustalo" in skill_licenses

    script_dirs = [
        PROJECT_ROOT / "scripts",
        SKILL_SUBDIR / "scripts",
    ]
    for sdir in script_dirs:
        for py_file in sdir.glob("*.py"):
            content = py_file.read_text(encoding="utf-8")
            assert "Licensed under the Apache License, Version 2.0" in content, f"Missing Apache header in {py_file}"

    doc_files = [
        PROJECT_ROOT / "README.md",
        PROJECT_ROOT / "SKILL.md",
        PROJECT_ROOT / "SKILL_LICENSES.md",
        PROJECT_ROOT / "documents" / "prd.md",
        PROJECT_ROOT / "references" / "color_schemes.md",
        PROJECT_ROOT / "references" / "dashboard_schema.md",
        SKILL_SUBDIR / "SKILL.md",
    ]
    for doc in doc_files:
        assert doc.is_file(), f"Missing documentation file: {doc}"

    overview_img = PROJECT_ROOT / "assets" / "skill_overview.png"
    assert overview_img.is_file() and overview_img.stat().st_size > 0
    readme_text = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
    assert "assets/skill_overview.png" in readme_text
