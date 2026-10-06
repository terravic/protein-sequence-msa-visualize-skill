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
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "polite-http",
#     "python-dotenv",
# ]
# ///
"""Enhanced EMBL-EBI Clustal Omega REST API client and deterministic mock engine.

Supports 1.0 QPS rate-limited job submission, status polling, multi-format
artifact retrieval (FASTA alignment, Percent Identity Matrix, Newick guide tree),
workspace license verification, dry-run payload inspection, and offline mock
execution.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional, Tuple
import urllib.error
import urllib.parse
import urllib.request

try:
    from dotenv import load_dotenv as _dotenv_load
except ImportError:
    _dotenv_load = None

try:
    from polite_http.http_client import HttpClient as _PoliteHttpClientExternal
except ImportError:
    _PoliteHttpClientExternal = None

try:
    import requests as _requests
except ImportError:
    _requests = None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("msa_align")

EBI_CLUSTALO_BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
MAX_FILE_SIZE_BYTES = 4 * 1024 * 1024  # 4 MB
MIN_SEQUENCES = 2
MAX_SEQUENCES = 4000
POLL_INTERVAL_SECONDS = 10.0
MAX_POLL_TIMEOUT_SECONDS = 15 * 60  # 15 minutes
DEFAULT_QPS = 1.0

LICENSE_NOTIFICATION_TEXT = (
    "Notification: Please check the EMBL-EBI terms of use at "
    "https://www.ebi.ac.uk/about/terms-of-use/ and Clustal Omega service "
    "documentation at https://www.ebi.ac.uk/jdispatcher/msa/clustalo. "
    "Always provide a valid contact email via USER_EMAIL to comply with "
    "EBI Job Dispatcher policies."
)


class MSAValidationError(ValueError):
    """Raised when input FASTA or environment parameters fail validation."""


class MSAApiError(RuntimeError):
    """Raised when the EMBL-EBI Clustal Omega REST API returns an error."""


class _FallbackHttpResponse:
    """Normalized HTTP response wrapper for fallback client."""

    def __init__(self, status_code: int, text: str) -> None:
        self.status_code = status_code
        self.text = text

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise MSAApiError(
                f"HTTP {self.status_code} error from EBI endpoint: {self.text[:300]}"
            )


class PoliteHttpClient:
    """Rate-limited HTTP client enforcing 1.0 QPS politeness for EMBL-EBI."""

    def __init__(self, qps: float = DEFAULT_QPS, timeout: float = 30.0) -> None:
        if qps <= 0:
            raise ValueError("qps must be positive.")
        self.qps = qps
        self.min_interval = 1.0 / qps
        self.timeout = timeout
        self._last_request_ts: float = 0.0
        self._external_client = (
            _PoliteHttpClientExternal(qps=qps)
            if _PoliteHttpClientExternal is not None
            else None
        )

    def _throttle(self) -> None:
        now = time.monotonic()
        elapsed = now - self._last_request_ts
        if self._last_request_ts > 0 and elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last_request_ts = time.monotonic()

    def post(self, url: str, data: Dict[str, Any]) -> Any:
        if self._external_client is not None and hasattr(_PoliteHttpClientExternal, "post"):
            return self._external_client.post(url, data=data)
        self._throttle()
        if _requests is not None:
            resp = _requests.post(url, data=data, timeout=self.timeout)
            return _FallbackHttpResponse(resp.status_code, resp.text)
        encoded = urllib.parse.urlencode(data).encode("utf-8")
        req = urllib.request.Request(url, data=encoded, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                body = resp.read().decode("utf-8", errors="replace")
                return _FallbackHttpResponse(resp.status, body)
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            return _FallbackHttpResponse(exc.code, body)

    def get(self, url: str) -> Any:
        if self._external_client is not None and hasattr(_PoliteHttpClientExternal, "get"):
            return self._external_client.get(url)
        self._throttle()
        if _requests is not None:
            resp = _requests.get(url, timeout=self.timeout)
            return _FallbackHttpResponse(resp.status_code, resp.text)
        req = urllib.request.Request(url, method="GET")
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                body = resp.read().decode("utf-8", errors="replace")
                return _FallbackHttpResponse(resp.status, body)
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            return _FallbackHttpResponse(exc.code, body)


def resolve_skill_root() -> Path:
    """Return the root directory of the skill containing scripts/ and sample_data/."""
    return Path(__file__).resolve().parent.parent


def resolve_workspace_root() -> Path:
    """Locate the top-level repository workspace root."""
    skill_root = resolve_skill_root()
    # If located inside <workspace>/skills/protein_sequence_msa_visualize/
    if skill_root.parent.name == "skills":
        return skill_root.parent.parent
    return skill_root


def ensure_license_guard(workspace_root: Optional[Path] = None) -> Path:
    """Verify or create the EMBL-EBI license notification file in .licenses/."""
    root = workspace_root or resolve_workspace_root()
    licenses_dir = root / ".licenses"
    primary_license = licenses_dir / "protein_sequence_msa_visualize_LICENSE.txt"
    fallback_license = licenses_dir / "protein_sequence_msa_LICENSE.txt"

    if primary_license.exists():
        return primary_license
    if fallback_license.exists():
        return fallback_license

    licenses_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    content = f"{LICENSE_NOTIFICATION_TEXT}\nTimestamp: {timestamp}\n"
    primary_license.write_text(content, encoding="utf-8")
    logger.info("Notice: %s", LICENSE_NOTIFICATION_TEXT)
    logger.info("Created license verification file at %s", primary_license)
    return primary_license


def _parse_env_file(env_path: Path) -> None:
    """Fallback parser for .env files when python-dotenv is not installed."""
    if not env_path.is_file():
        return
    try:
        for raw_line in env_path.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            if line.startswith("export "):
                line = line[len("export ") :].strip()
            key, val = line.split("=", 1)
            key = key.strip()
            val = val.strip().strip("'").strip('"')
            if key and key not in os.environ:
                os.environ[key] = val
    except OSError:
        pass


def load_user_email() -> Optional[str]:
    """Load USER_EMAIL from environment variables, ~/.env, or workspace .env."""
    home_env = Path.home() / ".env"
    local_env = resolve_workspace_root() / ".env"
    if _dotenv_load is not None:
        if home_env.is_file():
            _dotenv_load(home_env, override=False)
        if local_env.is_file():
            _dotenv_load(local_env, override=False)
    else:
        _parse_env_file(home_env)
        _parse_env_file(local_env)

    email = os.environ.get("USER_EMAIL", "").strip()
    return email if email else None


def parse_fasta_records(fasta_text: str) -> List[Tuple[str, str, str]]:
    """Parse FASTA text into list of (seq_id, description, sequence).

    Raises:
        MSAValidationError: If FASTA syntax is malformed or contains empty sequences.
    """
    records: List[Tuple[str, str, str]] = []
    current_id: Optional[str] = None
    current_desc: str = ""
    current_chunks: List[str] = []

    for line_num, raw_line in enumerate(fasta_text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if current_id is not None:
                seq_str = "".join(current_chunks)
                if not seq_str:
                    raise MSAValidationError(
                        f"Empty sequence encountered for FASTA record '{current_id}'."
                    )
                records.append((current_id, current_desc, seq_str))
            header = line[1:].strip()
            if not header:
                raise MSAValidationError(
                    f"Malformed FASTA header on line {line_num}: missing sequence identifier."
                )
            parts = header.split(None, 1)
            current_id = parts[0]
            current_desc = header
            current_chunks = []
        else:
            if current_id is None:
                raise MSAValidationError(
                    f"Malformed FASTA input on line {line_num}: sequence data precedes header '>'."
                )
            cleaned = "".join(ch for ch in line if not ch.isspace())
            current_chunks.append(cleaned)

    if current_id is not None:
        seq_str = "".join(current_chunks)
        if not seq_str:
            raise MSAValidationError(
                f"Empty sequence encountered for FASTA record '{current_id}'."
            )
        records.append((current_id, current_desc, seq_str))

    return records


def validate_input_fasta(input_path: Path) -> Tuple[str, List[Tuple[str, str, str]]]:
    """Validate FASTA file existence, byte size <= 4 MB, and sequence count 2..4000."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA file not found: {input_path}")
    if not input_path.is_file():
        raise MSAValidationError(f"Input path is not a regular file: {input_path}")

    file_size = input_path.stat().st_size
    if file_size > MAX_FILE_SIZE_BYTES:
        raise MSAValidationError(
            f"Input file size ({file_size} bytes) exceeds maximum limit of "
            f"{MAX_FILE_SIZE_BYTES} bytes (4 MB)."
        )
    if file_size == 0:
        raise MSAValidationError(f"Input FASTA file is empty: {input_path}")

    raw_text = input_path.read_text(encoding="utf-8")
    records = parse_fasta_records(raw_text)

    num_seqs = len(records)
    if num_seqs < MIN_SEQUENCES:
        raise MSAValidationError(
            f"Input FASTA contains {num_seqs} sequence(s); at least {MIN_SEQUENCES} "
            "sequences are required for multiple sequence alignment."
        )
    if num_seqs > MAX_SEQUENCES:
        raise MSAValidationError(
            f"Input FASTA contains {num_seqs} sequences, exceeding the maximum limit "
            f"of {MAX_SEQUENCES} sequences."
        )

    for seq_id, _, seq in records:
        ungapped = seq.replace("-", "").replace(".", "")
        if not ungapped:
            raise MSAValidationError(
                f"Sequence '{seq_id}' contains no amino acid residues."
            )

    return raw_text, records


def format_fasta(records: List[Tuple[str, str, str]], line_width: int = 60) -> str:
    """Format (seq_id, description, aligned_seq) records into standard wrapped FASTA."""
    lines: List[str] = []
    for seq_id, desc, seq in records:
        header = desc if desc else seq_id
        lines.append(f">{header}")
        for i in range(0, len(seq), line_width):
            lines.append(seq[i : i + line_width])
    return "\n".join(lines) + "\n"


def _pairwise_needleman_wunsch(seq_a: str, seq_b: str) -> Tuple[str, str]:
    """Lightweight global pairwise alignment for moderate-length protein sequences."""
    m, n = len(seq_a), len(seq_b)
    if m * n > 250_000:
        # Fast center-padded fallback for very long sequences
        target_len = max(m, n)
        return seq_a.ljust(target_len, "-"), seq_b.ljust(target_len, "-")

    match_score = 3
    mismatch_score = -1
    gap_penalty = -2

    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = i * gap_penalty
    for j in range(1, n + 1):
        dp[0][j] = j * gap_penalty

    for i in range(1, m + 1):
        ca = seq_a[i - 1]
        row_prev = dp[i - 1]
        row_curr = dp[i]
        for j in range(1, n + 1):
            score = match_score if ca == seq_b[j - 1] else mismatch_score
            diag = row_prev[j - 1] + score
            up = row_prev[j] + gap_penalty
            left = row_curr[j - 1] + gap_penalty
            row_curr[j] = max(diag, up, left)

    aligned_a: List[str] = []
    aligned_b: List[str] = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = match_score if seq_a[i - 1] == seq_b[j - 1] else mismatch_score
            if dp[i][j] == dp[i - 1][j - 1] + score:
                aligned_a.append(seq_a[i - 1])
                aligned_b.append(seq_b[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty:
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
        else:
            aligned_a.append("-")
            aligned_b.append(seq_b[j - 1])
            j -= 1

    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b))


def _synthesize_msa(records: List[Tuple[str, str, str]]) -> List[Tuple[str, str, str]]:
    """Deterministically align arbitrary sequences for offline mock mode."""
    stripped = [
        (sid, desc, seq.replace("-", "").replace(".", "").upper())
        for sid, desc, seq in records
    ]
    # If all inputs already have identical aligned length with gaps, preserve them
    orig_lengths = {len(seq) for _, _, seq in records}
    if len(orig_lengths) == 1 and any("-" in seq for _, _, seq in records):
        return [(sid, desc, seq.upper()) for sid, desc, seq in records]

    # Center-star alignment anchored on the longest sequence
    anchor_idx = max(range(len(stripped)), key=lambda idx: len(stripped[idx][2]))
    anchor_seq = stripped[anchor_idx][2]

    # Track gap insertions before each position of anchor_seq (0..len(anchor_seq))
    max_insertions = [0] * (len(anchor_seq) + 1)
    per_seq_insertions: List[Tuple[List[str], List[str]]] = []

    for idx, (_, _, seq) in enumerate(stripped):
        if idx == anchor_idx:
            per_seq_insertions.append((list(anchor_seq), [""] * (len(anchor_seq) + 1)))
            continue
        aln_anchor, aln_other = _pairwise_needleman_wunsch(anchor_seq, seq)
        anchor_chars: List[str] = []
        ins_buckets: List[List[str]] = [[] for _ in range(len(anchor_seq) + 1)]
        anchor_pos = 0
        for ca, cb in zip(aln_anchor, aln_other):
            if ca == "-":
                ins_buckets[anchor_pos].append(cb)
            else:
                anchor_chars.append(cb)
                anchor_pos += 1
        ins_strs = ["".join(bucket) for bucket in ins_buckets]
        for pos, s in enumerate(ins_strs):
            if len(s) > max_insertions[pos]:
                max_insertions[pos] = len(s)
        per_seq_insertions.append((anchor_chars, ins_strs))

    aligned_records: List[Tuple[str, str, str]] = []
    for (sid, desc, _), (matched_chars, ins_strs) in zip(stripped, per_seq_insertions):
        pieces: List[str] = []
        for pos in range(len(anchor_seq)):
            pieces.append(ins_strs[pos].ljust(max_insertions[pos], "-"))
            pieces.append(matched_chars[pos])
        pieces.append(ins_strs[len(anchor_seq)].ljust(max_insertions[len(anchor_seq)], "-"))
        aligned_records.append((sid, desc, "".join(pieces)))

    return aligned_records


def _compute_pim_from_aligned(
    aligned_records: List[Tuple[str, str, str]]
) -> Dict[str, Any]:
    """Build PIM dictionary and Clustal-formatted PIM text from aligned records."""
    seq_ids = [sid for sid, _, _ in aligned_records]
    seqs = [seq for _, _, seq in aligned_records]
    n = len(seqs)
    matrix: List[List[float]] = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                matrix[i][j] = 100.0
                continue
            matches = 0
            compared = 0
            for ca, cb in zip(seqs[i], seqs[j]):
                if ca not in ("-", ".") and cb not in ("-", "."):
                    compared += 1
                    if ca == cb:
                        matches += 1
            pct = round((100.0 * matches / compared), 2) if compared > 0 else 0.0
            matrix[i][j] = pct

    pim_lines = ["# Percent Identity Matrix - generated by Clustal Omega pipeline"]
    for idx, sid in enumerate(seq_ids):
        row_str = " ".join(f"{val:6.2f}" for val in matrix[idx])
        pim_lines.append(f"{idx + 1}: {sid:<22} {row_str}")

    return {
        "description": "Percent Identity Matrix (PIM)",
        "sequence_ids": seq_ids,
        "matrix": matrix,
        "raw_pim_text": "\n".join(pim_lines) + "\n",
    }


def build_newick_from_pim(seq_ids: List[str], matrix: List[List[float]]) -> str:
    """Construct a valid hierarchical Newick guide tree using UPGMA clustering."""
    if not seq_ids:
        return "();"
    if len(seq_ids) == 1:
        return f"({seq_ids[0]}:0.01000);"

    clusters: Dict[int, Tuple[str, float, int]] = {
        i: (sid, 0.0, 1) for i, sid in enumerate(seq_ids)
    }
    dists: Dict[Tuple[int, int], float] = {}
    n = len(seq_ids)
    for i in range(n):
        for j in range(i + 1, n):
            ident = max(0.0, min(100.0, matrix[i][j]))
            dist = max(0.002, round(1.0 - (ident / 100.0), 5))
            dists[(i, j)] = dist
            dists[(j, i)] = dist

    next_id = n
    while len(clusters) > 1:
        active_ids = sorted(clusters.keys())
        best_pair = (active_ids[0], active_ids[1])
        best_dist = dists.get(best_pair, 0.5)
        for idx_a in range(len(active_ids)):
            for idx_b in range(idx_a + 1, len(active_ids)):
                pair = (active_ids[idx_a], active_ids[idx_b])
                d = dists.get(pair, 0.5)
                if d < best_dist:
                    best_dist = d
                    best_pair = pair

        a, b = best_pair
        newick_a, height_a, size_a = clusters.pop(a)
        newick_b, height_b, size_b = clusters.pop(b)
        new_height = best_dist / 2.0
        branch_a = max(0.001, new_height - height_a)
        branch_b = max(0.001, new_height - height_b)
        merged_newick = f"({newick_a}:{branch_a:.5f},{newick_b}:{branch_b:.5f})"

        for other_id, (_, _, _) in clusters.items():
            d_ao = dists.get((a, other_id), 0.5)
            d_bo = dists.get((b, other_id), 0.5)
            d_new = (d_ao * size_a + d_bo * size_b) / (size_a + size_b)
            dists[(next_id, other_id)] = d_new
            dists[(other_id, next_id)] = d_new

        clusters[next_id] = (merged_newick, new_height, size_a + size_b)
        next_id += 1

    final_newick = next(iter(clusters.values()))[0]
    return f"{final_newick};\n"


def execute_mock_alignment(
    records: List[Tuple[str, str, str]],
    sample_data_dir: Optional[Path] = None,
) -> Tuple[str, str, str]:
    """Run deterministic offline alignment using fixtures or synthetic alignment.

    Returns:
        Tuple of (aligned_fasta_str, pim_output_str, newick_tree_str).
    """
    sample_dir = sample_data_dir or (resolve_skill_root() / "sample_data")
    input_ids = [sid for sid, _, _ in records]
    globin_ids = [
        "sp|P69905|HBA_HUMAN",
        "sp|P68871|HBB_HUMAN",
        "sp|P02144|MYG_HUMAN",
        "sp|P02238|LGB1_SOYBN",
        "sp|Q9NPG2|NGB_HUMAN",
    ]

    if input_ids == globin_ids:
        fa_path = sample_dir / "sample_alignment_result.fa"
        pim_path = sample_dir / "sample_pim.json"
        tree_path = sample_dir / "sample_phylotree.dnd"
        if fa_path.is_file() and pim_path.is_file() and tree_path.is_file():
            logger.info("Loaded pre-computed globin fixtures from %s", sample_dir)
            return (
                fa_path.read_text(encoding="utf-8"),
                pim_path.read_text(encoding="utf-8"),
                tree_path.read_text(encoding="utf-8"),
            )

    logger.info(
        "Generating deterministic offline mock alignment for %d sequences.",
        len(records),
    )
    aligned_records = _synthesize_msa(records)
    aligned_fasta = format_fasta(aligned_records)
    pim_dict = _compute_pim_from_aligned(aligned_records)
    pim_json = json.dumps(pim_dict, indent=2) + "\n"
    newick_tree = build_newick_from_pim(pim_dict["sequence_ids"], pim_dict["matrix"])
    return aligned_fasta, pim_json, newick_tree


def submit_and_poll_ebi(
    raw_fasta: str,
    email: str,
    title: str = "protein_sequence_msa_visualize",
    http_client: Optional[Any] = None,
    poll_interval: float = POLL_INTERVAL_SECONDS,
    timeout_seconds: float = MAX_POLL_TIMEOUT_SECONDS,
) -> Tuple[str, str, str]:
    """Submit MSA job to EMBL-EBI Clustal Omega API, poll status, and fetch results."""
    client = http_client or PoliteHttpClient(qps=DEFAULT_QPS)
    run_url = f"{EBI_CLUSTALO_BASE_URL}/run"
    payload = {
        "email": email,
        "title": title,
        "sequence": raw_fasta,
        "stype": "protein",
        "outfmt": "fa",
    }

    logger.info("Submitting Clustal Omega job to %s (qps=%.1f)...", run_url, DEFAULT_QPS)
    run_resp = client.post(run_url, data=payload)
    if hasattr(run_resp, "raise_for_status"):
        run_resp.raise_for_status()
    job_id = run_resp.text.strip()
    if not job_id:
        raise MSAApiError("EBI Clustal Omega returned an empty job_id.")
    logger.info("Submitted EBI Clustal Omega job_id: %s", job_id)

    status_url = f"{EBI_CLUSTALO_BASE_URL}/status/{job_id}"
    start_ts = time.monotonic()
    while True:
        status_resp = client.get(status_url)
        if hasattr(status_resp, "raise_for_status"):
            status_resp.raise_for_status()
        status = status_resp.text.strip().upper()
        logger.info("Job %s status: %s", job_id, status)

        if status == "FINISHED":
            break
        if status in ("ERROR", "FAILURE", "NOT_FOUND"):
            raise MSAApiError(
                f"EBI Clustal Omega job {job_id} terminated with status: {status}"
            )
        if (time.monotonic() - start_ts) >= timeout_seconds:
            raise MSAApiError(
                f"EBI Clustal Omega job {job_id} timed out after {timeout_seconds} seconds."
            )
        time.sleep(poll_interval)

    fa_resp = client.get(f"{EBI_CLUSTALO_BASE_URL}/result/{job_id}/fa")
    if hasattr(fa_resp, "raise_for_status"):
        fa_resp.raise_for_status()
    aligned_fasta = fa_resp.text

    pim_text = ""
    try:
        pim_resp = client.get(f"{EBI_CLUSTALO_BASE_URL}/result/{job_id}/pim")
        if getattr(pim_resp, "status_code", 200) < 400:
            pim_text = pim_resp.text
    except Exception as exc:  # pragma: no cover - defensive fallback
        logger.warning("Could not retrieve optional PIM endpoint: %s", exc)

    tree_text = ""
    try:
        tree_resp = client.get(f"{EBI_CLUSTALO_BASE_URL}/result/{job_id}/phylotree")
        if getattr(tree_resp, "status_code", 200) < 400:
            tree_text = tree_resp.text
    except Exception as exc:  # pragma: no cover - defensive fallback
        logger.warning("Could not retrieve optional phylotree endpoint: %s", exc)

    if not pim_text or not tree_text:
        aligned_records = parse_fasta_records(aligned_fasta)
        computed_pim = _compute_pim_from_aligned(aligned_records)
        if not pim_text:
            pim_text = json.dumps(computed_pim, indent=2) + "\n"
        if not tree_text:
            tree_text = build_newick_from_pim(
                computed_pim["sequence_ids"], computed_pim["matrix"]
            )

    return aligned_fasta, pim_text, tree_text


def run_alignment(
    input_fasta: Path,
    output_fasta: Path,
    output_pim: Optional[Path] = None,
    output_tree: Optional[Path] = None,
    dry_run: bool = False,
    mock: bool = False,
    http_client: Optional[Any] = None,
) -> Dict[str, Any]:
    """Validate input, enforce license guard, and execute MSA alignment."""
    ensure_license_guard()
    raw_fasta, records = validate_input_fasta(input_fasta)
    email = load_user_email()

    if dry_run:
        payload_preview = {
            "endpoint": f"{EBI_CLUSTALO_BASE_URL}/run",
            "email": email or "<USER_EMAIL_NOT_SET>",
            "title": f"clustalo_{input_fasta.stem}",
            "num_sequences": len(records),
            "sequence_ids": [sid for sid, _, _ in records],
            "file_size_bytes": input_fasta.stat().st_size,
            "sequence_preview": raw_fasta[:240] + ("..." if len(raw_fasta) > 240 else ""),
        }
        print(json.dumps(payload_preview, indent=2))
        return payload_preview

    if mock:
        aligned_fasta, pim_data, tree_data = execute_mock_alignment(records)
    else:
        if not email:
            raise MSAValidationError(
                "USER_EMAIL environment variable is required for live EMBL-EBI Clustal Omega "
                "requests. Set USER_EMAIL in your environment or ~/.env, or pass --mock for "
                "offline execution."
            )
        aligned_fasta, pim_data, tree_data = submit_and_poll_ebi(
            raw_fasta=raw_fasta,
            email=email,
            title=f"clustalo_{input_fasta.stem}",
            http_client=http_client,
        )

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_fasta.write_text(aligned_fasta, encoding="utf-8")
    logger.info("Saved aligned FASTA to %s", output_fasta)

    if output_pim is not None:
        output_pim.parent.mkdir(parents=True, exist_ok=True)
        output_pim.write_text(pim_data, encoding="utf-8")
        logger.info("Saved Percent Identity Matrix to %s", output_pim)

    if output_tree is not None:
        output_tree.parent.mkdir(parents=True, exist_ok=True)
        output_tree.write_text(tree_data, encoding="utf-8")
        logger.info("Saved Newick phylogenetic guide tree to %s", output_tree)

    return {
        "input_fasta": str(input_fasta),
        "output_fasta": str(output_fasta),
        "output_pim": str(output_pim) if output_pim else None,
        "output_tree": str(output_tree) if output_tree else None,
        "num_sequences": len(records),
        "mock": mock,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Align multiple protein sequences via EMBL-EBI Clustal Omega REST API "
            "or deterministic offline mock engine."
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
        type=Path,
        required=False,
        help="Path where the output alignment file in FASTA format will be saved.",
    )
    parser.add_argument(
        "--output-pim",
        type=Path,
        default=None,
        help="Optional path to save the Percent Identity Matrix (PIM).",
    )
    parser.add_argument(
        "--output-tree",
        type=Path,
        default=None,
        help="Optional path to save the Newick phylogenetic guide tree.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the prepared POST payload and exit without contacting the server.",
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Bypass network calls and return pre-computed or deterministic synthetic alignments.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if not args.dry_run and args.output is None:
        parser.error("the following arguments are required: -o/--output")

    try:
        run_alignment(
            input_fasta=args.input,
            output_fasta=args.output or Path("dry_run_aligned.fa"),
            output_pim=args.output_pim,
            output_tree=args.output_tree,
            dry_run=args.dry_run,
            mock=args.mock,
        )
        return 0
    except (FileNotFoundError, MSAValidationError, MSAApiError) as exc:
        logger.error("%s", exc)
        return 1
    except Exception as exc:  # pragma: no cover
        logger.error("Unexpected error in msa_align: %s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
