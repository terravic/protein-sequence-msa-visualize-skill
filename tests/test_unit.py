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
"""Unit tests for quantitative MSA metrics, parsers, visualizers, and validation guards."""

from __future__ import annotations

import json
import math
from pathlib import Path
import pytest

from scripts.analyze_msa import (
    LOG2_20,
    analyze_alignment,
    compute_column_metrics,
    compute_identity_matrices,
    load_features_spec,
    project_functional_features,
    validate_aligned_records,
)
from scripts.msa_align import (
    MAX_FILE_SIZE_BYTES,
    MSAApiError,
    MSAValidationError,
    PoliteHttpClient,
    build_newick_from_pim,
    ensure_license_guard,
    execute_mock_alignment,
    format_fasta,
    parse_fasta_records,
    submit_and_poll_ebi,
    validate_input_fasta,
)
from scripts.visualize_alignment import (
    COLOR_PALETTES,
    extract_accession_and_label,
    get_conservation_bar_color,
    get_residue_color,
    render_static_alignment_svg,
)
from scripts.visualize_matrix import (
    compute_tree_layout,
    get_identity_cell_color,
    parse_newick,
    parse_pim_file,
    render_tree_svg,
)


def test_parse_fasta_records_valid() -> None:
    fasta = ">seq1 First sequence\nACDEFG\nHIKLMN\n>seq2 Second sequence\nACDEFGHIKLMN\n"
    records = parse_fasta_records(fasta)
    assert len(records) == 2
    assert records[0] == ("seq1", "seq1 First sequence", "ACDEFGHIKLMN")
    assert records[1] == ("seq2", "seq2 Second sequence", "ACDEFGHIKLMN")


def test_parse_fasta_records_malformed() -> None:
    with pytest.raises(MSAValidationError, match="precedes header"):
        parse_fasta_records("ACDEFG\n>seq1\nACDEFG")

    with pytest.raises(MSAValidationError, match="missing sequence identifier"):
        parse_fasta_records(">   \nACDEFG\n")

    with pytest.raises(MSAValidationError, match="Empty sequence"):
        parse_fasta_records(">seq1\n>seq2\nACDEFG\n")


def test_validate_input_fasta_guards(tmp_path: Path) -> None:
    # Missing file
    with pytest.raises(FileNotFoundError):
        validate_input_fasta(tmp_path / "nonexistent.fasta")

    # Single sequence (< 2)
    single_fa = tmp_path / "single.fasta"
    single_fa.write_text(">seq1\nACDEFGHIKLMNPQRSTVWY\n", encoding="utf-8")
    with pytest.raises(MSAValidationError, match="at least 2 sequences"):
        validate_input_fasta(single_fa)

    # Gap-only sequence
    gap_fa = tmp_path / "gaps.fasta"
    gap_fa.write_text(">seq1\n------\n>seq2\nACDEFG\n", encoding="utf-8")
    with pytest.raises(MSAValidationError, match="no amino acid residues"):
        validate_input_fasta(gap_fa)

    # Oversized file (> 4 MB)
    large_fa = tmp_path / "oversized.fasta"
    large_fa.write_bytes(b">s1\n" + (b"A" * (MAX_FILE_SIZE_BYTES + 100)))
    with pytest.raises(MSAValidationError, match="exceeds maximum limit"):
        validate_input_fasta(large_fa)


def test_validate_aligned_records_length_mismatch() -> None:
    records = [
        ("seq1", "seq1", "ACDEFG"),
        ("seq2", "seq2", "ACDEF"),
    ]
    with pytest.raises(MSAValidationError, match="Alignment length mismatch"):
        validate_aligned_records(records)


def test_shannon_entropy_and_conservation_math() -> None:
    # Column 0: Invariant 'H' across 4 sequences -> H=0.0, Conservation=1.0, Consensus='H'
    # Column 1: 3 'A', 1 'V' (75% A) -> Consensus='a' (>=50% and <80%)
    # Column 2: 2 'A', 2 'V' (50% A) -> H=1.0 bit, Conservation = 1 - 1/log2(20)
    # Column 3: 4 distinct hydrophobic ('I', 'L', 'V', 'F') -> 25% each -> H=2.0 bits, Consensus='+'
    # Column 4: 4 unrelated ('A', 'D', 'K', 'C') -> 25% each -> Consensus='.'
    # Column 5: All gaps '-' -> H=log2(20), Conservation=0.0, GapRate=1.0, Consensus='-'
    records = [
        ("s1", "s1", "HAAIA-"),
        ("s2", "s2", "HAALD-"),
        ("s3", "s3", "HAVVK-"),
        ("s4", "s4", "HVVFC-"),
    ]
    cols, consensus, invariant_count = compute_column_metrics(records)

    assert invariant_count == 1
    assert consensus == "Haa+.-"

    # Column 0 (Invariant H)
    assert cols[0]["shannon_entropy"] == 0.0
    assert cols[0]["conservation_score"] == 1.0
    assert cols[0]["is_invariant"] is True
    assert cols[0]["consensus_aa"] == "H"

    # Column 2 (2 A, 2 V -> 1.0 bit entropy)
    assert math.isclose(cols[2]["shannon_entropy"], 1.0, abs_tol=1e-4)
    expected_cons_col2 = round(1.0 - (1.0 / LOG2_20), 4)
    assert math.isclose(cols[2]["conservation_score"], expected_cons_col2, abs_tol=1e-4)

    # Column 5 (All gaps)
    assert cols[5]["gap_rate"] == 1.0
    assert cols[5]["conservation_score"] == 0.0
    assert cols[5]["consensus_aa"] == "-"


def test_four_identity_metrics_definitions() -> None:
    # Sequence A: --ACDEFG-- (unaligned len = 6, aligned len L = 10, nongap bounds 2..7)
    # Sequence B: MKACDWFG-- (unaligned len = 8, aligned len L = 10, nongap bounds 0..7)
    # Overlap region between A and B: start=max(2,0)=2, end=min(7,7)=7 -> 6 cols (terminal_gaps = 4)
    # Matches between A and B in cols 2..7: A, C, D, F, G = 5 matches (E vs W mismatch)
    records = [
        ("seqA", "seqA", "--ACDEFG--"),
        ("seqB", "seqB", "MKACDWFG--"),
    ]
    id_mats = compute_identity_matrices(records)

    # Global Identity: 5 / 10 = 50.0%
    assert id_mats["global"][0][1] == 50.0
    # Coverage Identity: 5 / min(6, 8) = 5 / 6 = 83.33%
    assert id_mats["coverage"][0][1] == 83.33
    # Overlap Identity: 5 / (10 - 4) = 5 / 6 = 83.33%
    assert id_mats["overlap"][0][1] == 83.33


def test_feature_coordinate_mapping_and_edge_cases() -> None:
    # seq1 has leading gap and internal gap: "-M-ACDEF"
    # Unaligned pos 1 = 'M' -> col index 1 (position 2)
    # Unaligned pos 2 = 'A' -> col index 3 (position 4)
    # Unaligned pos 4 = 'D' -> col index 5 (position 6)
    records = [
        ("sp|P12345|REF_SEQ", "sp|P12345|REF_SEQ", "-M-ACDEF"),
        ("sp|P67890|HOMOLOG", "sp|P67890|HOMOLOG", "MMKACEEF"),
    ]
    cols, _, _ = compute_column_metrics(records)
    spec = load_features_spec(
        json.dumps(
            {
                "reference_id": "P12345",
                "features": [
                    {"name": "Initiator Met", "position": 1, "type": "site"},
                    {"name": "Catalytic Asp", "position": 4, "type": "active_site"},
                    {"name": "Out of Bounds", "position": 99, "type": "invalid"},
                ],
            }
        )
    )
    projected = project_functional_features(records, cols, spec)
    assert projected["resolved_reference_id"] == "sp|P12345|REF_SEQ"
    mapped = projected["mapped_features"]
    assert len(mapped) == 3

    # Feature 1: Met at unaligned pos 1 -> col index 1 (pos 2), invariant M
    assert mapped[0]["status"] == "mapped"
    assert mapped[0]["alignment_column_index"] == 1
    assert mapped[0]["alignment_column_position"] == 2
    assert mapped[0]["reference_aa"] == "M"
    assert mapped[0]["is_invariant"] is True
    assert mapped[0]["conservation_state"] == "invariant"

    # Feature 2: Asp at unaligned pos 4 -> col index 5 (pos 6), D vs E
    assert mapped[1]["status"] == "mapped"
    assert mapped[1]["alignment_column_index"] == 5
    assert mapped[1]["alignment_column_position"] == 6
    assert mapped[1]["reference_aa"] == "D"
    assert mapped[1]["residues_across_msa"]["sp|P67890|HOMOLOG"] == "E"

    # Feature 3: Out of bounds position 99 -> gracefully marked unmapped
    assert mapped[2]["status"] == "unmapped"
    assert mapped[2]["alignment_column_index"] is None


def test_newick_parser_and_tree_svg() -> None:
    newick = "(((A:0.15,B:0.20):0.05,C:0.30):0.10,(D:0.25,E:0.35):0.08);"
    tree = parse_newick(newick)
    layout = compute_tree_layout(tree)
    assert layout["num_leaves"] == 5
    assert layout["leaves"] == ["A", "B", "C", "D", "E"]
    svg = render_tree_svg(newick)
    assert "<svg" in svg
    assert 'data-seq-id="A"' in svg
    assert 'data-seq-id="E"' in svg


def test_color_schemes_and_alignment_svg() -> None:
    for scheme in ("ClustalX", "Zappo", "Taylor", "Hydrophobicity", "Monochromatic"):
        assert scheme in COLOR_PALETTES
        color = get_residue_color("H", scheme=scheme)
        assert color.startswith("#") and len(color) == 7

    assert get_conservation_bar_color(0.90) == "#16A34A"
    assert get_conservation_bar_color(0.55) == "#D97706"
    assert get_conservation_bar_color(0.15) == "#94A3B8"

    info = extract_accession_and_label("sp|P69905|HBA_HUMAN")
    assert info["accession"] == "P69905"
    assert info["name"] == "HBA_HUMAN"

    assert get_identity_cell_color(95.0) == "#059669"
    assert get_identity_cell_color(10.0) == "#1E3A8A"


def test_license_guard_creation(tmp_path: Path) -> None:
    lic_path = ensure_license_guard(workspace_root=tmp_path)
    assert lic_path.is_file()
    text = lic_path.read_text(encoding="utf-8")
    assert "https://www.ebi.ac.uk/about/terms-of-use/" in text
    assert "Timestamp:" in text


def test_mock_ebi_submission_and_polling() -> None:
    class MockResponse:
        def __init__(self, status_code: int, text: str) -> None:
            self.status_code = status_code
            self.text = text

        def raise_for_status(self) -> None:
            if self.status_code >= 400:
                raise MSAApiError(f"HTTP {self.status_code}")

    class MockHttpClient:
        def __init__(self) -> None:
            self.status_calls = 0

        def post(self, url: str, data: dict) -> MockResponse:
            assert url.endswith("/run")
            assert data["email"] == "test@example.org"
            return MockResponse(200, "clustalo-R20261005-123456")

        def get(self, url: str) -> MockResponse:
            if "/status/" in url:
                self.status_calls += 1
                if self.status_calls == 1:
                    return MockResponse(200, "RUNNING")
                return MockResponse(200, "FINISHED")
            if url.endswith("/fa"):
                return MockResponse(200, ">s1\nACDEFG\n>s2\nACDEFG\n")
            if url.endswith("/pim"):
                return MockResponse(200, "1: s1 100.00 100.00\n2: s2 100.00 100.00\n")
            if url.endswith("/phylotree"):
                return MockResponse(200, "(s1:0.01,s2:0.01);")
            return MockResponse(404, "Not found")

    client = MockHttpClient()
    fa, pim, tree = submit_and_poll_ebi(
        raw_fasta=">s1\nACDEFG\n>s2\nACDEFG\n",
        email="test@example.org",
        http_client=client,
        poll_interval=0.01,
    )
    assert client.status_calls == 2
    assert ">s1" in fa
    assert "100.00" in pim
    assert "s1:0.01" in tree
