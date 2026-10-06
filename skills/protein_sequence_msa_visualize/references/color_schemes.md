# Biochemical Amino Acid Color Schemes Reference

This document specifies the hex color codes and physicochemical groupings implemented by `scripts/visualize_alignment.py` and rendered in `dashboard.html`.

---

## 1. ClustalX Palette (Default)

The ClustalX color scheme groups amino acids by side-chain physicochemical properties.

| Group | Amino Acids (1-Letter Codes) | Hex Color | Visual Hue |
| :--- | :--- | :--- | :--- |
| **Hydrophobic / Aliphatic** | `A`, `I`, `L`, `M`, `F`, `W`, `V` | `#1F77B4` | Blue |
| **Polar / Uncharged / Glycine** | `N`, `Q`, `S`, `T`, `G` | `#2CA02C` | Green |
| **Acidic (Negatively Charged)** | `D`, `E` | `#E377C2` | Magenta |
| **Basic (Positively Charged)** | `K`, `R` | `#D62728` | Red |
| **Cysteine / Proline** | `C`, `P` | `#FF7F0E` | Orange |
| **Aromatic / Histidine / Tyrosine** | `H`, `Y` | `#17BECF` | Cyan |
| **Gap / Unknown** | `-`, `.`, `X` | `#94A3B8` | Slate Gray |

---

## 2. Zappo Palette

The Zappo scheme highlights physicochemical classes with high-contrast colors for rapid secondary and domain inspection.

| Group | Amino Acids (1-Letter Codes) | Hex Color | Visual Hue |
| :--- | :--- | :--- | :--- |
| **Aliphatic / Hydrophobic** | `I`, `L`, `V`, `A`, `M` | `#FFA07A` | Light Salmon |
| **Aromatic** | `F`, `W`, `Y` | `#FFA07A` | Light Salmon |
| **Positively Charged (Basic)** | `K`, `R`, `H` | `#0000FF` | Blue |
| **Negatively Charged (Acidic)** | `D`, `E` | `#FF0000` | Red |
| **Hydrophilic / Polar** | `S`, `T`, `N`, `Q` | `#00FF00` | Green |
| **Conformationally Special** | `P`, `G` | `#FF00FF` | Magenta |
| **Cysteine (Disulfide Capable)** | `C` | `#FFFF00` | Yellow |
| **Gap / Unknown** | `-`, `.`, `X` | `#94A3B8` | Slate Gray |

---

## 3. Taylor Palette

The Taylor (1997) color scheme assigns a unique, continuous rainbow hue to each of the 20 standard amino acids ordered by steric and hydrophobic similarity.

| Amino Acid | Name | Hex Color |
| :--- | :--- | :--- |
| `A` | Alanine | `#CCFF00` |
| `R` | Arginine | `#0000FF` |
| `N` | Asparagine | `#CC00FF` |
| `D` | Aspartic Acid | `#FF0000` |
| `C` | Cysteine | `#FFFF00` |
| `Q` | Glutamine | `#FF00CC` |
| `E` | Glutamic Acid | `#FF0066` |
| `G` | Glycine | `#FF9900` |
| `H` | Histidine | `#0066FF` |
| `I` | Isoleucine | `#66FF00` |
| `L` | Leucine | `#33FF00` |
| `K` | Lysine | `#6600FF` |
| `M` | Methionine | `#00FF00` |
| `F` | Phenylalanine | `#00FF66` |
| `P` | Proline | `#FFCC00` |
| `S` | Serine | `#FF3300` |
| `T` | Threonine | `#FF6600` |
| `W` | Tryptophan | `#00CCFF` |
| `Y` | Tyrosine | `#00FFCC` |
| `V` | Valine | `#99FF00` |
| `-` | Gap | `#94A3B8` |

---

## 4. Kyte-Doolittle Hydrophobicity Scale

Maps each residue along a continuous Blue (hydrophilic, score -4.5) to Red (hydrophobic, score +4.5) gradient based on the Kyte and Doolittle (1982) hydropathy index.

| Amino Acid | Hydropathy Index | Hex Color | Polarity Regime |
| :--- | :--- | :--- | :--- |
| `I` | `+4.5` | `#D7191C` | Strongly Hydrophobic |
| `V` | `+4.2` | `#DC2B24` | Strongly Hydrophobic |
| `L` | `+3.8` | `#E3432F` | Strongly Hydrophobic |
| `F` | `+2.8` | `#EF7546` | Moderately Hydrophobic |
| `C` | `+2.5` | `#F48850` | Moderately Hydrophobic |
| `M` | `+1.9` | `#F9AA60` | Weakly Hydrophobic |
| `A` | `+1.8` | `#FBB265` | Weakly Hydrophobic |
| `G` | `-0.4` | `#E6F598` | Neutral |
| `T` | `-0.7` | `#D5SR9E` / `#D1EC9C` | Neutral / Weakly Polar |
| `S` | `-0.8` | `#C7E89E` | Neutral / Weakly Polar |
| `W` | `-0.9` | `#BBE3A1` | Neutral / Weakly Polar |
| `Y` | `-1.3` | `#99D5A4` | Weakly Hydrophilic |
| `P` | `-1.6` | `#80C7A5` | Weakly Hydrophilic |
| `H` | `-3.2` | `#4393C3` | Moderately Hydrophilic |
| `E` | `-3.5` | `#3682BA` | Strongly Hydrophilic |
| `Q` | `-3.5` | `#3682BA` | Strongly Hydrophilic |
| `D` | `-3.5` | `#3682BA` | Strongly Hydrophilic |
| `N` | `-3.5` | `#3682BA` | Strongly Hydrophilic |
| `K` | `-3.9` | `#2B6CB0` | Strongly Hydrophilic |
| `R` | `-4.5` | `#1E4E8C` | Strongly Hydrophilic |
| `-` | `0.0` | `#94A3B8` | Gap |
