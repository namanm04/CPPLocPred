# CPPLocPred

**Prediction of Cell-Penetrating Peptides and their Subcellular Localization**

CPPLocPred is a two-stage machine learning tool:
- **Stage 1** — Classifies peptides as CPP or Non-CPP using an ExtraTrees model (AAC features)
- **Stage 2** — Predicts subcellular localization (Cytoplasm, Nucleus, Mitochondria, Endo/Lysosome, Others) using CatBoost models (DDR features)

It also supports **motif scanning** of protein sequences using MERCI for location-specific CPP motifs.

---

## Requirements

```
python >= 3.6
pandas
scikit-learn
catboost
perl (for motif search only)
```

Install dependencies:
```bash
pip install pandas scikit-learn catboost
```

---

## Directory Structure

```
CPPLocPred/
├── CPPLocPred.py                   # main script
├── MERCI_motif_locator.pl          # required for motif search (Job 2)
├── CPP_vs_NonCPP_ET_AAC.pkl
├── Cytoplasm_DDR_CatBoost.pkl
├── Nucleus_DDR_CatBoost.pkl
├── Mitochondria_DDR_CatBoost.pkl
├── Endo_lyso_DDR_CatBoost.pkl
├── Other_DDR_CatBoost.pkl
├── Localization_thresholds.pkl
└── motifs/
    ├── None/
    │   ├── Cytoplasm_motif
    │   ├── Nucleus_motif
    │   ├── Mitochondria_motif
    │   ├── Endo_lysosome_motif
    │   └── Others_motif
    ├── Koolman/          (same 5 files)
    ├── Betts-Russell/    (same 5 files)
    └── Rasmol/           (same 5 files)
```

---

## Usage

```
python CPPLocPred.py -i INPUT -o OUTPUT [-j {1,2}] [options]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Input FASTA file | required |
| `-o` | Output CSV file | required |
| `-j` | Job: `1` = Prediction, `2` = Motif Search | `1` |

---

### Job 1 — CPP Prediction

```bash
python CPPLocPred.py -i input.fasta -o results.csv
python CPPLocPred.py -i input.fasta -o results.csv -j 1 -t 0.5 -m /path/to/models/
```

| Flag | Description | Default |
|------|-------------|---------|
| `-t` | CPP probability threshold | `0.44` |
| `-m` | Directory containing model `.pkl` files | `./` |

**Output file columns:**

| Column | Description |
|--------|-------------|
| ID | Sequence identifier |
| Sequence | Amino acid sequence |
| CPP_Probability | Stage 1 probability (ExtraTrees) |
| CPP_Prediction | `CPP` or `Non-CPP` |
| Cytoplasm_Probability | Stage 2 localization score |
| Nucleus_Probability | |
| Mitochondria_Probability | |
| Endo_lysosome_Probability | |
| Others_Probability | |
| Final_Localization | Predicted location(s), semicolon-separated |

---

### Job 2 — Motif Search

Produces two output files:
- `output.csv` — one row per sequence (summary)
- `output_hits.csv` — one row per motif hit (detail)

```bash
python CPPLocPred.py -i input.fasta -o motifs.csv -j 2 -l Nucleus
python CPPLocPred.py -i input.fasta -o motifs.csv -j 2 -l Mitochondria -c Koolman
python CPPLocPred.py -i input.fasta -o motifs.csv -j 2 -l Nucleus -c Rasmol --motif_dir /path/to/motifs
```

| Flag | Description | Default |
|------|-------------|---------|
| `-l` | Location: `Cytoplasm`, `Nucleus`, `Mitochondria`, `Endo_lysosome`, `Others` | `Cytoplasm` |
| `-c` | Motif class subfolder: `None`, `Koolman`, `Betts-Russell`, `Rasmol` | `None` |
| `--motif_dir` | Directory containing motif class subfolders | `./motifs` |
| `--perl` | Path to perl executable | `perl` |

Motif file resolved as: `<motif_dir>/<class>/<location>_motif`

**Summary output (`output.csv`):**

| Column | Description |
|--------|-------------|
| ID | Sequence identifier |
| Sequence | Full amino acid sequence |
| Length | Sequence length |
| Total_Hits | Number of motif hits found |
| Motif_Patterns | Distinct motif patterns matched (`;` separated) |
| Hit_Positions | Start-end of each hit (`;` separated) |
| Matched_Residues | Matched amino acid strings (`;` separated) |

**Detail output (`output_hits.csv`):**

| Column | Description |
|--------|-------------|
| ID | Sequence identifier |
| Sequence | Full amino acid sequence |
| Start | Hit start position (1-based) |
| End | Hit end position |
| Hit_Length | Length of matched region |
| Motif_Pattern | MERCI motif pattern |
| Matched_Residues | Matched amino acid string |

---

## Input Format

Standard FASTA format:

```
>seq1
NALAALAKKRQIKIW
>seq2
RQIKIWFQNRRMKWKK
>seq3
GRKKRRQRRRPPQ
```

Only standard 20 amino acid single-letter codes (ACDEFGHIKLMNPQRSTVWY) are accepted.
Sequences with non-standard characters are flagged as `Invalid`.

---

## Examples

```bash
# Predict CPPs with default settings
python CPPLocPred.py -i peptides.fasta -o predictions.csv

# Predict with stricter threshold
python CPPLocPred.py -i peptides.fasta -o predictions.csv -j 1 -t 0.6

# Scan for Nucleus-targeting motifs
python CPPLocPred.py -i proteins.fasta -o nucleus_motifs.csv -j 2 -l Nucleus

# Scan for Mitochondria motifs using Koolman class
python CPPLocPred.py -i proteins.fasta -o mito_motifs.csv -j 2 -l Mitochondria -c Koolman

# Scan with custom motif directory
python CPPLocPred.py -i proteins.fasta -o mito_motifs.csv -j 2 -l Mitochondria -c Rasmol --motif_dir /data/motifs
```

---

## Citation

If you use CPPLocPred, please cite:

> Raghava et al. (2025) CPPLocPred: Machine learning-based prediction and subcellular localization of cell-penetrating peptides. IIIT Delhi.

---

## Web Server

[https://webs.iiitd.edu.in/raghava/cpplocpred/](https://webs.iiitd.edu.in/raghava/cpplocpred/)

## Contact

Raghava Group, IIIT Delhi — [https://webs.iiitd.edu.in/raghava/](https://webs.iiitd.edu.in/raghava/)