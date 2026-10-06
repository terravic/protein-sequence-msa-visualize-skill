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
"""Quantitative Multiple Sequence Alignment (MSA) metrics and feature mapper.

Computes per-column Shannon entropy, normalized conservation scores, gap rates,
consensus sequences, four sequence identity definitions (Sequence Coverage,
Global Identity, Overlap Identity, Conservation Index), and projects 1-based
unaligned reference sequence feature coordinates onto aligned MSA columns.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import logging
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

try:
    from scripts.msa_align import MSAValidationError, parse_fasta_records
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from msa_align import MSAValidationError, parse_fasta_records  # type: ignore

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("analyze_msa")

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
GAP_CHARS = frozenset("-. ")
LOG2_20 = math.log2(20.0)

BIOCHEMICAL_GROUPS: Dict[str, frozenset[str]] = {
    "hydrophobic": frozenset("AVLIMFWY"),
    "polar": frozenset("STNQ"),
    "acidic": frozenset("DE"),
    "basic": frozenset("KRH"),
    "special": frozenset("CGP"),
}


def validate_aligned_records(
    records: List[Tuple[str, str, str]]
) -> Tuple[int, int, List[Tuple[str, str, str]]]:
    """Verify at least 2 sequences and uniform aligned length L across all records."""
    if len(records) < 2:
        raise MSAValidationError(
            f"Aligned FASTA must contain at least 2 sequences; found {len(records)}."
        )
    normalized: List[Tuple[str, str, str]] = []
    expected_len = len(records[0][2])
    if expected_len == 0:
        raise MSAValidationError("Aligned sequences cannot have length 0.")

    for seq_id, desc, seq in records:
        if len(seq) != expected_len:
            raise MSAValidationError(
                f"Alignment length mismatch: sequence '{seq_id}' has length {len(seq)}, "
                f"expected {expected_len} (matching '{records[0][0]}')."
            )
        normalized.append((seq_id, desc, seq.upper()))

    return len(normalized), expected_len, normalized


def _has_physicochemical_similarity(aa_counts: Counter[str], total_nongap: int) -> bool:
    """Return True if >= 50% of non-gap residues fall into a single biochemical class."""
    if total_nongap <= 0:
        return False
    for group_set in BIOCHEMICAL_GROUPS.values():
        group_count = sum(cnt for aa, cnt in aa_counts.items() if aa in group_set)
        if (group_count / total_nongap) >= 0.50:
            return True
    return False


def compute_column_metrics(
    aligned_records: List[Tuple[str, str, str]]
) -> Tuple[List[Dict[str, Any]], str, int]:
    """Compute per-column Shannon entropy, conservation, gap rate, and consensus."""
    num_seqs = len(aligned_records)
    align_len = len(aligned_records[0][2])
    columns: List[Dict[str, Any]] = []
    consensus_chars: List[str] = []
    invariant_count = 0

    for col_idx in range(align_len):
        col_chars = [seq[col_idx] for _, _, seq in aligned_records]
        gap_count = sum(1 for c in col_chars if c in GAP_CHARS)
        nongap_chars = [c for c in col_chars if c not in GAP_CHARS]
        std_chars = [c for c in nongap_chars if c in STANDARD_AMINO_ACIDS]

        aa_counts: Counter[str] = Counter(nongap_chars)
        std_counts: Counter[str] = Counter(std_chars)
        num_nongap = len(nongap_chars)
        num_std = len(std_chars)
        gap_rate = round(gap_count / num_seqs, 4)

        if num_std == 0:
            entropy = round(LOG2_20, 4)
            raw_conservation = 0.0
            conservation = 0.0
        else:
            raw_h = 0.0
            for count in std_counts.values():
                if count > 0:
                    p_a = count / num_std
                    raw_h -= p_a * math.log2(p_a)
            entropy = round(max(0.0, raw_h), 4)
            raw_cons = max(0.0, min(1.0, 1.0 - (raw_h / LOG2_20)))
            raw_conservation = round(raw_cons, 4)
            conservation = round(raw_cons * (1.0 - (gap_count / num_seqs)), 4)

        if num_nongap == 0:
            modal_aa = "-"
            modal_freq = 0.0
            nongap_modal_freq = 0.0
            consensus_aa = "-"
            is_invariant = False
        else:
            # Deterministic tie-breaking: highest count, then alphabetical residue
            modal_aa, modal_count = sorted(
                aa_counts.items(), key=lambda kv: (-kv[1], kv[0])
            )[0]
            modal_freq = round(modal_count / num_seqs, 4)
            nongap_modal_freq = round(modal_count / num_nongap, 4)
            is_invariant = (gap_count == 0) and (modal_count == num_seqs)
            if is_invariant:
                invariant_count += 1

            if modal_freq >= 0.80:
                consensus_aa = modal_aa.upper()
            elif modal_freq >= 0.50:
                consensus_aa = modal_aa.lower()
            elif _has_physicochemical_similarity(aa_counts, num_nongap):
                consensus_aa = "+"
            else:
                consensus_aa = "."

        consensus_chars.append(consensus_aa)
        columns.append(
            {
                "index": col_idx,
                "position": col_idx + 1,
                "shannon_entropy": entropy,
                "conservation_score": conservation,
                "raw_conservation_score": raw_conservation,
                "gap_rate": gap_rate,
                "consensus_aa": consensus_aa,
                "modal_aa": modal_aa,
                "modal_frequency": modal_freq,
                "nongap_modal_frequency": nongap_modal_freq,
                "is_invariant": is_invariant,
                "aa_counts": dict(sorted(aa_counts.items())),
            }
        )

    return columns, "".join(consensus_chars), invariant_count


def _nongap_bounds(seq: str) -> Tuple[int, int, int]:
    """Return (first_nongap_idx, last_nongap_idx, unaligned_length) for aligned sequence."""
    first = -1
    last = -1
    count = 0
    for idx, ch in enumerate(seq):
        if ch not in GAP_CHARS:
            if first == -1:
                first = idx
            last = idx
            count += 1
    return first, last, count


def compute_identity_matrices(
    aligned_records: List[Tuple[str, str, str]]
) -> Dict[str, Any]:
    """Compute Global, Sequence Coverage, and Overlap identity matrices."""
    num_seqs = len(aligned_records)
    align_len = len(aligned_records[0][2])
    seq_ids = [sid for sid, _, _ in aligned_records]
    seqs = [seq for _, _, seq in aligned_records]
    bounds = [_nongap_bounds(seq) for seq in seqs]

    global_mat = [[0.0] * num_seqs for _ in range(num_seqs)]
    coverage_mat = [[0.0] * num_seqs for _ in range(num_seqs)]
    overlap_mat = [[0.0] * num_seqs for _ in range(num_seqs)]
    pairwise_details: List[Dict[str, Any]] = []

    for i in range(num_seqs):
        first_i, last_i, len_i = bounds[i]
        for j in range(num_seqs):
            first_j, last_j, len_j = bounds[j]
            matches = 0
            internal_gaps = 0

            overlap_start = max(first_i, first_j) if (first_i >= 0 and first_j >= 0) else 0
            overlap_end = min(last_i, last_j) if (last_i >= 0 and last_j >= 0) else -1

            if overlap_end >= overlap_start:
                terminal_gaps = overlap_start + (align_len - 1 - overlap_end)
                overlap_len = align_len - terminal_gaps
            else:
                terminal_gaps = align_len
                overlap_len = 0

            for col in range(align_len):
                ca = seqs[i][col]
                cb = seqs[j][col]
                a_gap = ca in GAP_CHARS
                b_gap = cb in GAP_CHARS
                if not a_gap and not b_gap and ca == cb:
                    matches += 1
                if overlap_start <= col <= overlap_end and (a_gap or b_gap):
                    internal_gaps += 1

            if i == j:
                global_id = 100.0
                coverage_id = 100.0
                overlap_id = 100.0
            else:
                global_id = round(100.0 * matches / max(1, align_len), 2)
                min_len = min(len_i, len_j)
                coverage_id = (
                    round(100.0 * matches / min_len, 2) if min_len > 0 else 0.0
                )
                overlap_id = (
                    round(100.0 * matches / overlap_len, 2)
                    if overlap_len > 0
                    else 0.0
                )

            global_mat[i][j] = global_id
            coverage_mat[i][j] = coverage_id
            overlap_mat[i][j] = overlap_id

            if i < j:
                pairwise_details.append(
                    {
                        "seq_a": seq_ids[i],
                        "seq_b": seq_ids[j],
                        "matches": matches,
                        "aligned_length": align_len,
                        "min_unaligned_length": min(len_i, len_j),
                        "overlap_length": overlap_len,
                        "terminal_gaps": terminal_gaps,
                        "internal_gaps": internal_gaps,
                        "global_identity": global_id,
                        "coverage_identity": coverage_id,
                        "overlap_identity": overlap_id,
                    }
                )

    return {
        "sequence_ids": seq_ids,
        "global": global_mat,
        "coverage": coverage_mat,
        "overlap": overlap_mat,
        "pairwise_details": pairwise_details,
    }


def _resolve_reference_sequence(
    reference_id: str, aligned_records: List[Tuple[str, str, str]]
) -> Optional[Tuple[str, str, str]]:
    """Find reference sequence by exact ID, UniProt token, or case-insensitive substring."""
    ref_clean = reference_id.strip()
    if not ref_clean:
        return aligned_records[0] if aligned_records else None

    for rec in aligned_records:
        if rec[0] == ref_clean:
            return rec
    for rec in aligned_records:
        tokens = rec[0].split("|")
        if ref_clean in tokens:
            return rec
    ref_lower = ref_clean.lower()
    for rec in aligned_records:
        if ref_lower in rec[0].lower() or ref_lower in rec[1].lower():
            return rec
    return None


def load_features_spec(features_arg: Optional[str]) -> Optional[Dict[str, Any]]:
    """Parse --features from either a JSON string or a JSON file path."""
    if not features_arg:
        return None
    stripped = features_arg.strip()
    if not stripped:
        return None

    candidate_path = Path(stripped)
    if not stripped.startswith("{") and candidate_path.is_file():
        raw_json = candidate_path.read_text(encoding="utf-8")
    else:
        raw_json = stripped

    try:
        parsed = json.loads(raw_json)
    except json.JSONDecodeError as exc:
        raise MSAValidationError(
            f"Invalid JSON provided to --features: {exc}"
        ) from exc

    if not isinstance(parsed, dict):
        raise MSAValidationError("--features JSON root must be an object.")
    return parsed


def project_functional_features(
    aligned_records: List[Tuple[str, str, str]],
    columns: List[Dict[str, Any]],
    features_spec: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Map 1-based unaligned residue positions on a reference sequence to MSA columns."""
    if not features_spec:
        return {
            "reference_id": None,
            "resolved_reference_id": None,
            "mapped_features": [],
        }

    ref_query = str(features_spec.get("reference_id", "")).strip()
    raw_features = features_spec.get("features", [])
    if not isinstance(raw_features, list):
        raise MSAValidationError("'features' key in --features must be a list.")

    resolved = _resolve_reference_sequence(ref_query, aligned_records)
    if resolved is None:
        return {
            "reference_id": ref_query,
            "resolved_reference_id": None,
            "mapped_features": [
                {
                    "name": str(feat.get("name", "Unnamed Feature")),
                    "type": str(feat.get("type", "functional_site")),
                    "unaligned_position": feat.get("position"),
                    "alignment_column_index": None,
                    "alignment_column_position": None,
                    "status": "unmapped",
                    "error": f"Reference sequence '{ref_query}' not found in alignment.",
                }
                for feat in raw_features
                if isinstance(feat, dict)
            ],
        }

    ref_id, _, ref_aligned_seq = resolved
    # Build 1-based unaligned coordinate -> 0-based MSA column index lookup
    unaligned_to_col: Dict[int, int] = {}
    unaligned_counter = 0
    for col_idx, ch in enumerate(ref_aligned_seq):
        if ch not in GAP_CHARS:
            unaligned_counter += 1
            unaligned_to_col[unaligned_counter] = col_idx

    mapped_list: List[Dict[str, Any]] = []
    for feat in raw_features:
        if not isinstance(feat, dict):
            continue
        name = str(feat.get("name", "Functional Site"))
        feat_type = str(feat.get("type", "functional_site"))
        desc = str(feat.get("description", ""))
        pos = feat.get("position")

        if not isinstance(pos, int) or pos < 1 or pos not in unaligned_to_col:
            mapped_list.append(
                {
                    "name": name,
                    "type": feat_type,
                    "description": desc,
                    "unaligned_position": pos,
                    "alignment_column_index": None,
                    "alignment_column_position": None,
                    "reference_aa": None,
                    "column_conservation": None,
                    "column_entropy": None,
                    "is_invariant": False,
                    "conservation_state": "unmapped",
                    "status": "unmapped",
                    "error": (
                        f"Position {pos} is out of bounds for reference '{ref_id}' "
                        f"(unaligned length {unaligned_counter})."
                    ),
                    "residues_across_msa": {},
                }
            )
            continue

        col_idx = unaligned_to_col[pos]
        col_meta = columns[col_idx]
        ref_aa = ref_aligned_seq[col_idx]
        residues_across = {
            sid: seq[col_idx] for sid, _, seq in aligned_records
        }
        unique_nongap = {
            aa for aa in residues_across.values() if aa not in GAP_CHARS
        }
        has_gaps = any(aa in GAP_CHARS for aa in residues_across.values())

        if not has_gaps and len(unique_nongap) == 1:
            cons_state = "invariant"
        elif col_meta["conservation_score"] >= 0.65:
            cons_state = "conserved"
        else:
            cons_state = "divergent"

        mapped_list.append(
            {
                "name": name,
                "type": feat_type,
                "description": desc,
                "unaligned_position": pos,
                "alignment_column_index": col_idx,
                "alignment_column_position": col_idx + 1,
                "reference_aa": ref_aa,
                "column_conservation": col_meta["conservation_score"],
                "column_entropy": col_meta["shannon_entropy"],
                "is_invariant": col_meta["is_invariant"],
                "conservation_state": cons_state,
                "status": "mapped",
                "residues_across_msa": residues_across,
            }
        )

    return {
        "reference_id": ref_query or ref_id,
        "resolved_reference_id": ref_id,
        "mapped_features": mapped_list,
    }


def analyze_alignment(
    alignment_path: Path,
    features_arg: Optional[str] = None,
    output_json_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run full quantitative MSA analysis and return structured dictionary."""
    if not alignment_path.is_file():
        raise FileNotFoundError(f"Aligned FASTA file not found: {alignment_path}")

    raw_text = alignment_path.read_text(encoding="utf-8")
    raw_records = parse_fasta_records(raw_text)
    num_seqs, align_len, aligned_records = validate_aligned_records(raw_records)

    columns, consensus_seq, invariant_cols = compute_column_metrics(aligned_records)
    identity_matrices = compute_identity_matrices(aligned_records)
    features_spec = load_features_spec(features_arg)
    projected_features = project_functional_features(
        aligned_records, columns, features_spec
    )

    pairwise = identity_matrices["pairwise_details"]
    num_pairs = max(1, len(pairwise))
    mean_global = round(
        sum(p["global_identity"] for p in pairwise) / num_pairs, 2
    )
    mean_coverage = round(
        sum(p["coverage_identity"] for p in pairwise) / num_pairs, 2
    )
    mean_overlap = round(
        sum(p["overlap_identity"] for p in pairwise) / num_pairs, 2
    )
    mean_cons = round(
        sum(c["conservation_score"] for c in columns) / max(1, align_len), 4
    )
    mean_gap = round(
        sum(c["gap_rate"] for c in columns) / max(1, align_len), 4
    )
    cons_index_pct = round(100.0 * invariant_cols / max(1, align_len), 2)

    sequences_meta = []
    for sid, desc, seq in aligned_records:
        gaps = sum(1 for ch in seq if ch in GAP_CHARS)
        sequences_meta.append(
            {
                "id": sid,
                "description": desc,
                "aligned_sequence": seq,
                "unaligned_length": align_len - gaps,
                "aligned_length": align_len,
                "gap_count": gaps,
            }
        )

    result = {
        "summary": {
            "num_sequences": num_seqs,
            "alignment_length": align_len,
            "invariant_columns": invariant_cols,
            "conservation_index": cons_index_pct,
            "conservation_index_ratio": round(invariant_cols / max(1, align_len), 4),
            "mean_conservation_score": mean_cons,
            "mean_gap_rate": mean_gap,
            "mean_global_identity": mean_global,
            "mean_coverage_identity": mean_coverage,
            "mean_overlap_identity": mean_overlap,
            "consensus_sequence": consensus_seq,
        },
        "sequences": sequences_meta,
        "columns": columns,
        "identity_matrices": identity_matrices,
        "features": projected_features,
    }

    if output_json_path is not None:
        output_json_path.parent.mkdir(parents=True, exist_ok=True)
        output_json_path.write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8"
        )
        logger.info("Saved MSA analytical metrics to %s", output_json_path)

    return result


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute quantitative conservation, Shannon entropy, pairwise identity "
            "metrics, and functional feature coordinate projections from an aligned FASTA."
        )
    )
    parser.add_argument(
        "--alignment",
        type=Path,
        required=True,
        help="Path to aligned FASTA file.",
    )
    parser.add_argument(
        "--features",
        type=str,
        default=None,
        help="Optional JSON string or JSON file path specifying functional sites on a reference sequence.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        required=True,
        help="Path to save computed analytical metrics in JSON format.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        analyze_alignment(
            alignment_path=args.alignment,
            features_arg=args.features,
            output_json_path=args.output_json,
        )
        return 0
    except (FileNotFoundError, MSAValidationError) as exc:
        logger.error("%s", exc)
        return 1
    except Exception as exc:  # pragma: no cover
        logger.error("Unexpected error in analyze_msa: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
