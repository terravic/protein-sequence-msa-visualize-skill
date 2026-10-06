# Dashboard and Analytical Metrics JSON Schema Specification

This document defines the JSON contract produced by `scripts/analyze_msa.py` (`metrics.json`) and embedded inside `<script id="msa-data" type="application/json">` by `scripts/generate_dashboard.py`.

---

## 1. Top-Level `metrics.json` Contract

```json
{
  "summary": {
    "num_sequences": 5,
    "alignment_length": 165,
    "invariant_columns": 12,
    "conservation_index": 7.27,
    "mean_conservation_score": 0.4812,
    "mean_gap_rate": 0.1152,
    "mean_global_identity": 31.88,
    "mean_coverage_identity": 35.42,
    "mean_overlap_identity": 33.15,
    "consensus_sequence": "mV-LSpADkTnVkAAWGkVG..."
  },
  "sequences": [
    {
      "id": "sp|P69905|HBA_HUMAN",
      "description": "sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens",
      "aligned_sequence": "MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR------",
      "unaligned_length": 142,
      "aligned_length": 165,
      "gap_count": 23
    }
  ],
  "columns": [
    {
      "index": 0,
      "position": 1,
      "shannon_entropy": 0.7219,
      "conservation_score": 0.833,
      "raw_conservation_score": 0.833,
      "gap_rate": 0.0,
      "consensus_aa": "m",
      "modal_aa": "M",
      "modal_frequency": 0.8,
      "nongap_modal_frequency": 0.8,
      "is_invariant": false,
      "aa_counts": {
        "M": 4,
        "V": 1
      }
    }
  ],
  "identity_matrices": {
    "sequence_ids": [
      "sp|P69905|HBA_HUMAN",
      "sp|P68871|HBB_HUMAN"
    ],
    "global": [
      [100.0, 41.21],
      [41.21, 100.0]
    ],
    "coverage": [
      [100.0, 47.89],
      [47.89, 100.0]
    ],
    "overlap": [
      [100.0, 43.59],
      [43.59, 100.0]
    ],
    "pairwise_details": [
      {
        "seq_a": "sp|P69905|HBA_HUMAN",
        "seq_b": "sp|P68871|HBB_HUMAN",
        "matches": 68,
        "aligned_length": 165,
        "min_unaligned_length": 142,
        "overlap_length": 156,
        "terminal_gaps": 9,
        "internal_gaps": 16,
        "global_identity": 41.21,
        "coverage_identity": 47.89,
        "overlap_identity": 43.59
      }
    ]
  },
  "features": {
    "reference_id": "sp|P69905|HBA_HUMAN",
    "resolved_reference_id": "sp|P69905|HBA_HUMAN",
    "mapped_features": [
      {
        "name": "Distal Histidine (E7)",
        "type": "active_site",
        "unaligned_position": 59,
        "alignment_column_index": 64,
        "alignment_column_position": 65,
        "reference_aa": "H",
        "column_conservation": 1.0,
        "column_entropy": 0.0,
        "is_invariant": true,
        "status": "mapped",
        "residues_across_msa": {
          "sp|P69905|HBA_HUMAN": "H",
          "sp|P68871|HBB_HUMAN": "H",
          "sp|P02144|MYG_HUMAN": "H",
          "sp|P02238|LGB1_SOYBN": "H",
          "sp|Q9NPG2|NGB_HUMAN": "H"
        }
      }
    ]
  }
}
```

---

## 2. Embedded Dashboard Data Payload (`#msa-data`)

When `scripts/generate_dashboard.py` compiles `dashboard.html`, it embeds a single JSON document inside `<script id="msa-data" type="application/json">` with the following top-level keys:

| Key | Type | Description |
| :--- | :--- | :--- |
| `title` | `string` | Human-readable dashboard title. |
| `generated_at` | `string` | ISO 8601 UTC generation timestamp. |
| `summary` | `object` | Summary statistics from `metrics.json`. |
| `sequences` | `array` | Ordered array of aligned sequence objects. |
| `columns` | `array` | Per-column entropy, conservation, consensus, and residue distribution objects. |
| `identity_matrices` | `object` | `global`, `coverage`, `overlap`, and optional `external_pim` $N \times N$ matrices plus `pairwise_details`. |
| `features` | `object` | Reference coordinate mappings and per-sequence residue projections. |
| `phylotree` | `object` | Raw Newick string (`newick`) and parsed hierarchical tree structure (`root`, `leaves`, `max_depth`). |
| `color_palettes` | `object` | Biochemical color scheme dictionaries (`ClustalX`, `Zappo`, `Taylor`, `Hydrophobicity`, `Monochromatic`). |
| `citations` | `array` | Mandatory academic citations (Sievers et al. 2011; Madeira et al. 2024; Goujon et al. 2010). |
| `disclaimer` | `string` | Scientific AS-IS disclaimer text. |
