# Product Requirements Document: Protein Sequence MSA Visualize Skill

## 1. Document Metadata

| Field | Value |
| :--- | :--- |
| **Skill Identifier** | `protein-sequence-msa-visualize-skill` |
| **Package Location** | Root directory and `skills/protein_sequence_msa_visualize/` |
| **License** | Apache License 2.0 |
| **Upstream Service** | EMBL-EBI Clustal Omega REST API (`https://www.ebi.ac.uk/Tools/services/rest/clustalo`) |
| **Runtime Requirement** | Python 3.10 or higher |
| **Status** | Active |

---

## 2. User Stories and Scope

### 2.1 Target Audience

1. **Computational Biologists and Bioinformaticians (Technical Users)**
   - Require automated, rate-limited command-line pipelines to align homologous protein families, compute quantitative per-column Shannon entropy and pairwise identity metrics across multiple mathematical definitions, project catalytic or post-translational modification (PTM) coordinates onto alignment columns, and export machine-readable JSON/CSV artifacts.
2. **Structural Biologists, Protein Engineers, and Life Science Researchers (Domain Users)**
   - Require an interactive, self-contained visual interface (`dashboard.html`) to inspect input sequence composition, conserved motifs, biochemical coloring schemes (ClustalX, Zappo, Taylor, Hydrophobicity), 3D sequence homology relationships, all-vs-all percent identity matrices (PIM), and phylogenetic guide trees in a standard web browser.
3. **Automated Agent Workflows (System Users)**
   - Require deterministic CLI execution with structured stdout summaries, clear exit codes, offline mock execution (`--mock`) for sandboxed testing, and compliance with EMBL-EBI usage policies (`qps=1.0`, `USER_EMAIL` tracking, and license verification guards).

### 2.2 User Stories

- **US-01 (End-to-End Alignment and Visualization)**: As a protein engineer, I want to provide a multi-sequence FASTA file and receive both the aligned FASTA and a standalone HTML dashboard so that I can inspect conserved domains and pairwise sequence identities immediately in a browser.
- **US-02 (Multi-Metric Identity Comparison)**: As a computational biologist, I want pairwise sequence identity calculated across Sequence Coverage, Global Alignment Length, Overlap (excluding terminal gaps), and Multisequence Conservation Index so that length disparities or N/C-terminal extensions do not distort homology assessments.
- **US-03 (Functional Residue Projection)**: As a biochemist studying an enzyme family, I want to specify 1-based unaligned residue coordinates on a reference sequence (such as a catalytic triad or heme-binding histidine) and map them automatically to aligned MSA columns so that I can determine whether active-site residues are invariant or mutated across homologs.
- **US-04 (Interactive Biochemical Color and 3D Homology Inspection)**: As a structural biologist, I want to toggle between ClustalX, Zappo, Taylor, and Kyte-Doolittle Hydrophobicity color palettes, inspect a 3D homology network and conservation landscape, and search for peptide motifs via regular expressions so that physicochemical conservation patterns are clear.
- **US-05 (Sandboxed Offline Execution)**: As a system developer running in a network-restricted environment, I want a deterministic `--mock` mode that executes the complete analytical and visualization pipeline using local fixtures or synthetic alignment generation without external network egress.
- **US-06 (Service Governance)**: As a compliance administrator, I want automatic verification of EMBL-EBI terms of use, enforced 1.0 QPS rate limiting, academic citations, and synthetic or public reference sequences in sample datasets.

### 2.3 Scope Definition

| Capability Category | In Scope | Out of Scope |
| :--- | :--- | :--- |
| **Sequence Input** | Protein FASTA files containing 2 to 4,000 sequences up to 4 MB total size. | Nucleotide (DNA/RNA) genome-scale alignments; raw sequencing reads (FASTQ/BAM). |
| **Alignment Engine** | EMBL-EBI Clustal Omega REST API (`fa`, `pim`, `phylotree` endpoints) and deterministic local mock alignment engine (`--mock`). | Local compilation or bundling of C/C++ Clustal Omega binaries; structural 3D superposition (PDB/mmCIF). |
| **Quantitative Metrics** | Input sequence composition and GRAVY hydropathy, per-column Shannon entropy, normalized conservation score, gap rate, consensus sequence generation, and four identity definitions (Coverage, Global, Overlap, Conservation Index). | Bayesian phylogenetic tree inference (e.g., MrBayes, BEAST) or bootstrap resampling across thousands of trees. |
| **Feature Mapping** | Coordinate projection from 1-based unaligned reference sequence positions to 1-based and 0-based MSA column indices with cross-sequence residue breakdown. | Automated de novo domain prediction against external HMM databases (e.g., Pfam/InterProScan API calls). |
| **Visualization** | Standalone HTML5/SVG/JS dashboard (`dashboard.html`) with input sequence composition table, 3D sequence homology network and conservation landscape, dual-axis scrolling alignment viewer, conservation bar chart, interactive PIM heatmap, Newick cladogram, motif search, and FASTA/CSV/SVG export. | Multi-user server-hosted database persistence or collaborative real-time editing backends. |

---

## 3. Technical Architecture

### 3.1 System Sequence and Data Flow Diagram

```text
+-----------------------------------------------------------------------------------+
|                        AGENT HARNESS / CLI USER INVOCATION                        |
|        python3 scripts/run_pipeline.py <INPUT_FASTA> -o <OUT_DIR> [OPTIONS]       |
+-----------------------------------------+-----------------------------------------+
                                          |
                                          v
+-----------------------------------------------------------------------------------+
| 1. PRE-FLIGHT GOVERNANCE & INPUT VALIDATION (run_pipeline.py / msa_align.py)      |
|    - Check/create .licenses/protein_sequence_msa_visualize_LICENSE.txt            |
|    - Validate file exists, size <= 4 MB, sequence count 2 <= N <= 4000            |
|    - Verify non-empty amino acid sequences and unique/resolvable headers          |
+--------------------+-----------------------------------------+--------------------+
                     |                                         |
         [Live Mode: Default]                         [Offline Mode: --mock]
                     |                                         |
                     v                                         v
+----------------------------------------+   +--------------------------------------+
| 2A. EMBL-EBI CLUSTAL OMEGA REST API    |   | 2B. DETERMINISTIC MOCK ENGINE        |
|     (scripts/msa_align.py)             |   |     (scripts/msa_align.py --mock)    |
|     - Load USER_EMAIL from env / .env  |   |     - Match known fixtures in        |
|     - Enforce PoliteHttpClient qps=1.0 |   |       sample_data/ (globins, etc.)   |
|     - POST /clustalo/run -> job_id     |   |     - Or synthesize deterministic    |
|     - Poll /status/{job_id} (10s step) |   |       gapped alignment, pairwise     |
|     - Fetch /result/{job_id}/fa        |   |       PIM, and neighbor-joining      |
|     - Fetch /result/{job_id}/pim       |   |       Newick guide tree              |
|     - Fetch /result/{job_id}/phylotree |   +------------------+-------------------+
+--------------------+-------------------+                      |
                     |                                          |
                     +--------------------+---------------------+
                                          |
                                          v
                     Artifacts Written to Output Directory:
                     - aligned.fa (Aligned FASTA)
                     - matrix.pim (Percent Identity Matrix)
                     - tree.dnd   (Newick Phylogenetic Guide Tree)
                                          |
                                          v
+-----------------------------------------------------------------------------------+
| 3. QUANTITATIVE ANALYTICAL ENGINE (scripts/analyze_msa.py)                        |
|    - Verify uniform alignment length L across all N sequences                     |
|    - Compute per-column amino acid distribution, gap rate, Shannon entropy H(j),  |
|      and normalized conservation index C(j) = 1 - H(j)/log2(20)                   |
|    - Build consensus sequence (>=80% upper, >=50% lower, <50% '+' or '.')         |
|    - Compute N x N pairwise identity matrices:                                    |
|      * Global Identity:   Matches(A,B) / L                                        |
|      * Coverage Identity: Matches(A,B) / min(Len(A), Len(B))                      |
|      * Overlap Identity:  Matches(A,B) / (L - TerminalGaps(A,B))                  |
|      * Conservation Index: InvariantColumns / L                                   |
|    - Project 1-based unaligned reference features (--features) to MSA columns     |
|    - Emit structured JSON -> metrics.json                                         |
+-----------------------------------------+-----------------------------------------+
                                          |
                                          v
+-----------------------------------------------------------------------------------+
| 4. VISUALIZATION & DASHBOARD COMPILER                                             |
|    - scripts/visualize_alignment.py:                                              |
|      Renders multi-track viewer, ruler, conservation bars, and palette definitions|
|      (ClustalX, Zappo, Taylor, Kyte-Doolittle Hydrophobicity)                     |
|    - scripts/visualize_matrix.py:                                                 |
|      Renders interactive N x N PIM heatmap and parses Newick tree into SVG        |
|      cladogram with branch length scaling and leaf hover hooks                    |
|    - scripts/generate_dashboard.py:                                               |
|      Computes input sequence composition profiles, embeds JSON payload in         |
|      <script id="msa-data" type="application/json">, injects standalone CSS       |
|      variables (light/dark), renders 3D homology network (#homology-3d-viewport), |
|      binds interactive controls (motif regex, jump-to-column, threshold slider,   |
|      SVG/CSV/FASTA exporters), and appends academic citations -> dashboard.html   |
+-----------------------------------------+-----------------------------------------+
                                          |
                                          v
+-----------------------------------------------------------------------------------+
| 5. EXECUTIVE SUMMARY OUTPUT (stdout)                                              |
|    - Sequence count N, aligned length L, conservation index, mean global identity |
|    - Top conserved columns and mapped functional feature status                   |
|    - Relative paths to generated aligned.fa, metrics.json, and dashboard.html     |
+-----------------------------------------------------------------------------------+
```

### 3.2 Module Inventory

1. `scripts/msa_align.py`: CLI client for EMBL-EBI Clustal Omega with 1.0 QPS rate limiting, input validation, `.licenses/` verification, `--dry-run` payload inspection, multi-endpoint retrieval (`fa`, `pim`, `phylotree`), and `--mock` offline execution.
2. `scripts/analyze_msa.py`: Quantitative calculation engine producing `metrics.json` containing column-level entropy, conservation scores, gap rates, consensus sequence, four identity metric definitions, and projected feature coordinates.
3. `scripts/visualize_alignment.py`: Component generator for the synchronized dual-axis alignment track, conservation bar chart, biochemical color palettes, column inspector, and motif highlighter.
4. `scripts/visualize_matrix.py`: Component generator for the all-vs-all pairwise identity heatmap and Newick phylogenetic tree parser/SVG cladogram renderer.
5. `scripts/generate_dashboard.py`: Standalone HTML5 dashboard assembler combining input sequence composition profiles, data payloads, 3D sequence homology network and conservation landscape, visual modules, theme controls, export utilities, scientific disclaimers, and literature citations into `dashboard.html`.
6. `scripts/run_pipeline.py`: Unified end-to-end orchestrator chaining license check, alignment, quantitative analysis, dashboard compilation, and executive summary reporting.

---

## 4. Functional Requirements

### 4.1 Input Requirements and Validation Rules

- **FR-01 (FASTA Format Parsing)**: The parser must accept standard single-line and multi-line FASTA files where each record begins with `>` followed by an identifier and optional description, and subsequent lines contain amino acid characters (plus `-` or `.` in pre-aligned files).
- **FR-02 (File Size Limit)**: Input FASTA files must not exceed `4,194,304` bytes (4 MB). Files exceeding 4 MB must be rejected prior to network submission with a descriptive error message and non-zero exit code.
- **FR-03 (Sequence Cardinality Bounds)**: Input FASTA files must contain at least 2 sequences and at most 4,000 sequences ($2 \le N \le 4000$). Inputs outside this range must be rejected immediately.
- **FR-04 (License Guard Enforcement)**: Before submitting API jobs or generating outputs, the pipeline must verify the existence of `.licenses/protein_sequence_msa_visualize_LICENSE.txt` (or `.licenses/protein_sequence_msa_LICENSE.txt`) in the project root. If absent, the tool must print the EMBL-EBI terms notice to stderr/stdout and create `.licenses/protein_sequence_msa_visualize_LICENSE.txt` with an ISO 8601 UTC timestamp.
- **FR-05 (Email Registration in Live Mode)**: When running against the live EMBL-EBI REST API (neither `--mock` nor `--dry-run`), `USER_EMAIL` must be loaded from environment variables or `.env`. If missing or empty, execution must halt with instructions on setting `USER_EMAIL`.

### 4.2 Mathematical and Algorithmic Processing Stages

- **FR-06 (Column Shannon Entropy and Conservation Index)**:
  For an alignment of $N$ sequences and uniform aligned length $L$, at each column $j \in \{0, \dots, L-1\}$:
  - Let $A_{std} = \{\text{A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V}\}$ be the 20 standard amino acids.
  - Let $N_{\text{non-gap}}(j)$ be the count of standard amino acid residues at column $j$, and $N_{\text{gap}}(j)$ be the count of gap characters (`-`, `.`).
  - Relative frequency for each $a \in A_{std}$:
    $$p_j(a) = \frac{\text{count}_j(a)}{N_{\text{non-gap}}(j)} \quad \text{for } N_{\text{non-gap}}(j) > 0$$
  - Shannon entropy (in bits):
    $$H(j) = -\sum_{a \in A_{std}, p_j(a) > 0} p_j(a) \log_2 p_j(a)$$
    If $N_{\text{non-gap}}(j) = 0$, $H(j) = \log_2(20)$ and conservation is $0.0$. If all non-gap residues are identical, $H(j) = 0.0$.
  - Normalized residue conservation score:
    $$C_{\text{raw}}(j) = \max\left(0.0, \min\left(1.0, 1.0 - \frac{H(j)}{\log_2(20)}\right)\right)$$
  - Gap penalty weighting: Both raw non-gap conservation $C_{\text{raw}}(j)$ and gap-adjusted conservation $C(j) = C_{\text{raw}}(j) \times (1 - \text{GapRate}(j))$ are recorded, where $\text{GapRate}(j) = \frac{N_{\text{gap}}(j)}{N}$. When a column has zero gaps ($N_{\text{gap}}(j) = 0$), $C(j) = C_{\text{raw}}(j) = 1 - \frac{H(j)}{\log_2(20)}$.

- **FR-07 (Consensus Sequence Rules)**:
  At each column $j$:
  - Determine the modal amino acid $a^*$ and its frequency $f^* = \frac{\text{count}_j(a^*)}{N}$ (and non-gap frequency $f_{\text{ng}}^* = \frac{\text{count}_j(a^*)}{N_{\text{non-gap}}(j)}$).
  - If column consists entirely of gaps ($N_{\text{non-gap}}(j) = 0$), emit `-`.
  - If $f^* \ge 0.80$, emit uppercase $a^*$ (e.g., `H`).
  - Else if $f^* \ge 0.50$, emit lowercase $a^*$ (e.g., `h`).
  - Else if all non-gap residues belong to the same conserved physicochemical group (hydrophobic `AVLIMFWY`, polar `STNQ`, acidic `DE`, basic `KRH`) or $\ge 50\%$ share a physicochemical group, emit `+`.
  - Otherwise, emit `.`.

- **FR-08 (Pairwise and Multisequence Identity Metrics)**:
  For every pair of aligned sequences $(S_u, S_v)$ of aligned length $L$:
  - Let $\text{Matches}(S_u, S_v)$ be the number of columns $j \in [0, L-1]$ where $S_u[j] = S_v[j]$ and $S_u[j] \notin \{\text{'-'}, \text{'.'}\}$.
  - Let $\text{Len}(S_u)$ and $\text{Len}(S_v)$ be the unaligned lengths (count of non-gap characters) of $S_u$ and $S_v$.
  - Let $\text{TerminalGaps}(S_u, S_v)$ be the count of columns prior to $\max(\text{first\_nongap}(S_u), \text{first\_nongap}(S_v))$ plus columns after $\min(\text{last\_nongap}(S_u), \text{last\_nongap}(S_v))$.
  - **Global Identity (%)**:
    $$\text{ID}_{\text{global}}(S_u, S_v) = 100.0 \times \frac{\text{Matches}(S_u, S_v)}{L}$$
  - **Sequence Coverage Identity (%)**:
    $$\text{ID}_{\text{coverage}}(S_u, S_v) = 100.0 \times \frac{\text{Matches}(S_u, S_v)}{\min(\text{Len}(S_u), \text{Len}(S_v))}$$
  - **Overlap Identity (%)**:
    $$\text{ID}_{\text{overlap}}(S_u, S_v) = 100.0 \times \frac{\text{Matches}(S_u, S_v)}{\max(1, L - \text{TerminalGaps}(S_u, S_v))}$$
  - **Multisequence Conservation Index (%)**:
    $$\text{CI}_{\text{msa}} = 100.0 \times \frac{\text{Count of Strictly Invariant Non-Gap Columns}}{L}$$

- **FR-09 (Functional Feature Coordinate Projection)**:
  - Input format: JSON string or path to JSON file containing `reference_id` and a list of `features`, each with `name`, 1-based unaligned `position`, and optional `type` / `description`.
  - Reference matching: Match `reference_id` by exact sequence ID, accession substring (e.g., `P69905` matching `sp|P69905|HBA_HUMAN`), or prefix match.
  - Coordinate projection: Traverse the aligned sequence characters $j = 0 \dots L-1$, incrementing unaligned counter $k$ whenever $S_{\text{ref}}[j] \notin \{\text{'-'}, \text{'.'}\}$. When $k == \text{position}$, record 0-based column index $j$ and 1-based column number $j + 1$.
  - Cross-sequence profile: Extract the residue at column $j$ for every sequence in the alignment, determine whether the site is invariant, conserved, or divergent, and include the breakdown in `metrics.json` and `dashboard.html`.

### 4.3 Output Artifacts

- **FR-10 (`aligned.fa`)**: Standard multi-sequence FASTA alignment where all sequences have uniform length $L$.
- **FR-11 (`matrix.pim` / `sample_pim.json`)**: Pairwise identity matrix output from EBI or mock engine.
- **FR-12 (`tree.dnd`)**: Newick-formatted hierarchical clustering tree with branch lengths.
- **FR-13 (`metrics.json`)**: Complete analytical JSON conforming to `references/dashboard_schema.md` and `assets/metrics_schema.json`.
- **FR-14 (`dashboard.html` & Interactive UI Artifact Delivery)**: Standalone single-file HTML5 dashboard structured for native HTML/CSS/JS rendering inside an `<iframe>` preview pane and local browser:
  - Self-contained CSS with light and dark theme variables (`--bg-primary`, `--bg-secondary`, `--text-main`, `--text-muted`, `--border-color`, `--accent-color`).
  - **Light and Dark Theme Toggle (`#toggle-theme-btn`)**: Header button containing inline SVG Sun and Moon vector graphics for switching between light and dark themes.
  - Embedded `<script id="msa-data" type="application/json">` payload containing processed user sequence data.
  - **Input Sequences & Residue Composition Panel (`#input-data-profile-section`)**: Displays raw unaligned input sequence lengths, alignment gap additions, GRAVY hydropathy scores, and stacked biochemical class bars (Hydrophobic, Polar, Basic, Acidic, Special).
  - **3D Sequence Homology Network & Conservation Landscape (`#homology-3d-viewport`)**: Interactive 3D SVG viewport (`viewBox="0 0 680 390"`) pre-populated at compile time via `render_initial_3d_network_svg()` after 60 force-directed 3D relaxation steps so the network renders immediately on load, and paired with client-side 3D perspective projection, orbital reference rings, drag-to-orbit rotation, scroll-to-zoom, identity threshold filtering, node click selection, and a mode selector for the 3D Conservation Landscape ribbon.
  - **Interactive Multi-Track Sequence Alignment Viewer**: ClustalX, Zappo, Taylor, Hydrophobicity, and Monochromatic schemes, per-position Conservation Bar Chart, dynamic Consensus Track with adjustable threshold slider (50% to 100%), regex motif search (`Find Motif`), and column jump (`Jump to Column`).
  - **Interactive All-vs-All PIM Heatmap & SVG Phylogenetic Guide Tree Cladogram**: Symmetric $N \times N$ matrix with metric selector (Global, Coverage, Overlap, Clustal PIM) and Newick tree with bidirectional leaf-node highlighting.
  - **Functional Feature & Motif Coordinate Mapper**: Cards and live coordinate projector mapping 1-based unaligned reference residues onto aligned MSA columns.
  - Client-side export buttons for Aligned FASTA (`.fa`), Alignment SVG (`.svg`), and Identity Matrix CSV/JSON (`.csv`, `.json`).
  - Academic citations (Sievers et al. 2011; Madeira et al. 2024; Goujon et al. 2010) and the AS-IS scientific disclaimer.

---

## 5. Non-Functional Requirements

- **NFR-01 (API Politeness and Rate Limiting)**: All HTTP requests to EMBL-EBI must pass through a rate limiter enforcing at most 1.0 queries per second (`qps=1.0`, minimum 1.0-second interval between requests) and a 10-second polling interval with a 15-minute maximum job timeout.
- **NFR-02 (Performance Latency)**:
  - Offline mock pipeline (`--mock`) on standard protein families (5 to 50 sequences, $L \le 1000$) must complete end-to-end (alignment, metrics calculation, HTML generation) in under 2.0 seconds on standard CPU hardware.
  - Quantitative analysis (`analyze_msa.py`) must process 500 aligned sequences of length 500 in under 5.0 seconds using $O(N \cdot L)$ column operations and optimized pairwise comparisons.
- **NFR-03 (Data Privacy and Sequence Confidentiality)**:
  - Only public UniProt reference sequences or synthetic constructs are stored in `sample_data/`, `assets/`, or `tests/`.
  - `USER_EMAIL` is read exclusively from runtime environment variables or local `.env` files and is never hardcoded or written into generated HTML dashboards.
- **NFR-04 (Standalone File Portability)**: `dashboard.html` must open directly via the `file://` protocol in any modern browser without requiring a local HTTP server, CORS workarounds, or external API calls.
- **NFR-05 (Dependency Resilience)**: While `polite-http` and `python-dotenv` are primary dependencies, `scripts/msa_align.py` includes standard-library fallbacks for `.env` parsing and 1.0 QPS rate-limited HTTP transport so CLI execution succeeds in minimal Python 3.10+ environments.

---

## 6. Error Handling and Edge Cases

| Scenario ID | Input / System Condition | Detection Stage | Expected Behavior and Recovery | Exit Code |
| :--- | :--- | :--- | :--- | :--- |
| **ERR-01** | Input FASTA file path does not exist or is a directory | `msa_align.py` / `run_pipeline.py` pre-flight check | Log explicit `FileNotFoundError` message naming the missing path; terminate before any API or file write operations. | `1` |
| **ERR-02** | Input FASTA file size exceeds 4 MB (`> 4,194,304` bytes) | `msa_align.py` validation guard | Reject input with message stating actual byte size and 4 MB EBI Clustal Omega limit. | `1` |
| **ERR-03** | Input FASTA has fewer than 2 sequences ($N < 2$) | `msa_align.py` / `analyze_msa.py` parser | Reject input with message indicating at least 2 sequences are required for multiple sequence alignment. | `1` |
| **ERR-04** | Input FASTA has more than 4,000 sequences ($N > 4000$) | `msa_align.py` validation guard | Reject input with message stating sequence count and 4,000 sequence maximum limit. | `1` |
| **ERR-05** | FASTA contains empty sequence record or non-sequence text | FASTA parser in `msa_align.py` and `analyze_msa.py` | Reject malformed FASTA record identifying the offending header or line number. | `1` |
| **ERR-06** | `USER_EMAIL` missing when running live EBI mode | `msa_align.py` credential loader | Abort before network call; instruct user to set `USER_EMAIL` in environment or `.env` or pass `--mock` / `--dry-run`. | `1` |
| **ERR-07** | `.licenses/protein_sequence_msa_visualize_LICENSE.txt` missing | Pre-flight license guard | Automatically print EMBL-EBI terms notification and write `.licenses/protein_sequence_msa_visualize_LICENSE.txt` with UTC timestamp; continue execution. | `0` |
| **ERR-08** | Aligned FASTA passed to `analyze_msa.py` has unequal sequence lengths | `analyze_msa.py` length uniformity check | Raise validation error detailing sequence IDs with mismatched lengths; abort metric computation. | `1` |
| **ERR-09** | Alignment column contains 100% gaps (`N_non_gap == 0`) | `analyze_msa.py` entropy calculator | Guard against division by zero; assign `entropy = log2(20)`, `conservation = 0.0`, `gap_rate = 1.0`, and consensus character `'-'`. | `0` |
| **ERR-10** | Two sequences have zero overlapping non-gap positions (`Len == 0` or `L == TerminalGaps`) | `analyze_msa.py` pairwise identity calculator | Guard denominator with `max(1, denom)` or return `0.0%` identity without raising `ZeroDivisionError`. | `0` |
| **ERR-11** | `--features` specifies unknown `reference_id` or out-of-bounds residue position | `analyze_msa.py` feature mapper | Record feature status as `"unmapped"` with a clear diagnostic reason in `metrics.json` (or raise error if JSON syntax is invalid) without crashing the pipeline. | `0` (unmapped recorded) / `1` (malformed JSON) |
| **ERR-12** | EBI Job Dispatcher returns `ERROR`, `FAILURE`, `NOT_FOUND`, or times out after 15 minutes | `msa_align.py` polling loop | Log EBI job ID and terminal status; exit cleanly with descriptive runtime error message. | `1` |
| **ERR-13** | Optional `--pim` or `--tree` files are omitted when calling `generate_dashboard.py` | `generate_dashboard.py` input loader | Fall back to pairwise identity matrices in `metrics.json` and auto-generate a hierarchical distance guide tree from the identity matrix so all dashboard panels render. | `0` |
