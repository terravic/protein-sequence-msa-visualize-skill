---
name: protein-sequence-msa-visualize-skill
description: >
  Performs multiple sequence alignment (MSA) of protein sequences using EBI
  Clustal Omega and generates interactive standalone HTML dashboards with
  multi-track sequence viewers, dynamic conservation bars, all-vs-all percent
  identity heatmaps, and phylogenetic guide trees. Use when you need to align
  multiple proteins, assess residue conservation, map functional catalytic sites,
  or evaluate evolutionary similarity. Supports up to 4000 sequences and 4 MB.
---

# Protein Sequence MSA Visualize

## Overview
Provides automated multiple sequence alignment (MSA) of protein sequences via the EMBL-EBI Clustal Omega API, coupled with an analytical pipeline that calculates Shannon entropy conservation, pairwise identity matrices (Global, Coverage, Overlap), and consensus sequences. Compiles results into an interactive, standalone HTML dashboard (`dashboard.html`) featuring input sequence composition profiling, a 3D Sequence Homology Network and Conservation Landscape, multi-track sequence viewers, biochemical color schemes (ClustalX, Zappo, Taylor, Hydrophobicity), conservation bar charts, and Newick guide tree dendrograms.

## Prerequisites
1. **Python Environment**: Python 3.10+ with `polite-http` and `python-dotenv` installed.
2. **Contact Email (`USER_EMAIL`)**: EBI Job Dispatcher requires a contact email for job tracking. Ensure `USER_EMAIL` is set in your environment or `.env`.
3. **License Guard**: Verify that `.licenses/protein_sequence_msa_visualize_LICENSE.txt` exists in the project root. If absent, display the licensing terms notice and create the file with an ISO timestamp before querying.

## Core Rules
- **Execute Utility Scripts**: Always invoke the provided Python scripts in `scripts/`. Do not perform ad-hoc `curl` or unthrottled requests against the EBI API.
- **Rate Limiting Compliance**: Enforce the 1.0 QPS limit to respect EMBL-EBI service quotas.
- **Scientific Attribution**: Every report and dashboard must state that the alignment was performed using **EBI Clustal Omega**, citing Sievers et al. (2011) and Madeira et al. (2024).
- **Computed Data Only**: Report only calculated metrics and alignment columns directly computed from the alignment output.
- **Interactive UI Delivery**: Process the user's input sequence data first using the Python pipeline (`scripts/run_pipeline.py` or `scripts/generate_dashboard.py`) and save the resulting standalone HTML file to `dashboard.html` so it can be rendered directly in an HTML preview pane or web browser.

## Interactive Dashboard Architecture
The generated `dashboard.html` is structured for native HTML/CSS/JS rendering inside an `<iframe>` preview pane or web browser:
- **Light and Dark Theme Toggle**: Header toggle button (`#toggle-theme-btn`) displaying inline SVG Sun and Moon vector graphics for switching between light and dark themes.
- **Input Sequences & Residue Composition**: Compares raw unaligned input sequence lengths, alignment gap insertions, GRAVY hydropathy indices, and biochemical class distributions (Hydrophobic, Polar, Basic, Acidic, Special).
- **3D Sequence Homology Network**: Interactive 3D SVG viewport (`#homology-3d-viewport`) rendering sequences as 3D nodes connected by homology edges scaled by pairwise identity, with mode switching to a **3D Conservation Landscape** ribbon across alignment columns.
- **Interactive Multi-Track Sequence Alignment Viewer**: Dual-axis scrollable alignment matrix with ClustalX, Zappo, Taylor, Hydrophobicity, and Monochromatic palettes, per-position Shannon entropy bars, regex motif search, and column click inspection.
- **Pairwise Identity Matrix & Phylogenetic Guide Tree**: Symmetric $N \times N$ heatmap across four identity definitions paired with an interactive SVG Newick cladogram.
- **Functional Feature & Motif Coordinate Mapper**: Projects 1-based unaligned reference sequence coordinates onto aligned MSA columns.

## Workflows

### 1. Unified End-to-End Pipeline
Run full alignment, analytical profiling, and dashboard generation:
```bash
python3 scripts/run_pipeline.py <INPUT_FASTA> -o <OUTPUT_DIR>
```

Pipeline with known catalytic or functional feature mapping:
```bash
python3 scripts/run_pipeline.py <INPUT_FASTA> \
    --features '{"reference_id": "UniProt_ID", "features": [{"name": "Catalytic Dyad", "position": 25}]}' \
    -o <OUTPUT_DIR>
```

Offline synthetic mock mode (for sandboxed testing without external network access):
```bash
python3 scripts/run_pipeline.py sample_data/sample_globins.fasta -o ./output/mock_globins/ --mock
```

### 2. Modular Subcommands
Perform individual steps modularly:
```bash
# Step 1: Align sequences via EBI Clustal Omega
python3 scripts/msa_align.py input.fasta -o output/aligned.fa --output-pim output/matrix.pim --output-tree output/tree.dnd

# Step 2: Compute quantitative conservation and identity metrics
python3 scripts/analyze_msa.py --alignment output/aligned.fa --output-json output/metrics.json

# Step 3: Generate standalone HTML dashboard
python3 scripts/generate_dashboard.py --alignment output/aligned.fa --metrics output/metrics.json --pim output/matrix.pim --tree output/tree.dnd -o output/dashboard.html
```
