# NetMHCpan–IEDB Evaluation

An evaluation of NetMHCpan 4.1 binding affinity predictions against experimental peptide–MHC binding data from the Immune Epitope Database (IEDB).

## Overview

Major histocompatibility complex class I (MHC-I) molecules present peptide antigens to CD8+ T cells. Computational tools such as NetMHCpan predict peptide–MHC binding affinity, helping researchers identify potential T-cell epitopes.

This project evaluates how well NetMHCpan 4.1 binding affinity predictions agree with experimental measurements and how effectively predicted affinities distinguish stronger from weaker binders.

The analysis focuses on **HLA-A*02:01 and 9-mer peptides**.

## Objectives

- Compare predicted and experimentally measured peptide–MHC binding affinities.
- Quantify the association between predicted and experimental affinities.
- Evaluate binder classification at experimental IC50 thresholds of 500 nM and 50 nM.
- Compare predicted binding affinities across evidence groups: T-cell-positive, binding-positive, and peptides with both types of evidence.

## Methodology

### 1. Data preparation

The processed dataset contains predicted IC50 values, experimental affinity measurements where available, and evidence labels. Experimental measurements are used for the prediction-versus-experiment comparison and binder classification.

IC50 values are expressed in nanomolar (nM), with logarithmic transformations used for correlation and error analysis.

### 2. Statistical evaluation

The analysis calculates:

- **Spearman's rank correlation coefficient (ρ):** Measures the strength of the association between predicted and experimental affinities. A 95% confidence interval is estimated using bootstrap resampling.
- **Mean absolute error (MAE):** Quantifies the average absolute difference between predicted and experimental log10 IC50 values.
- **Linear regression:** Examines the relationship between predicted and experimental log-transformed affinities.
- **ROC-AUC:** Evaluates binder classification using experimental thresholds of 500 nM and 50 nM.
- **Mann–Whitney U test:** Compares predicted affinities between T-cell-positive and binding-positive peptide groups.

For binder classification, peptides below the experimental IC50 threshold are labelled binders. Lower predicted IC50 values indicate stronger predicted binding.

### 3. Visualisation

The script generates three figures:

1. **Predicted vs. experimental affinity:** A hexbin plot showing the relationship between predicted and experimental log10 IC50 values, with reference and regression lines.
2. **ROC curves:** Compares binder classification performance at 500 nM and 50 nM thresholds.
3. **Cumulative distributions by evidence group:** Compares predicted affinities across the available evidence categories.

## Repository structure

```text
NetMHCpan-IEDB-Benchmark/
├── data/
│   └── processed/
│       └── final_dataset.csv
├── scripts/
│   └── analysis.py
├── results/
│   └──                     # Generated figures
├── poster/
├── requirements.txt
└── README.md
```

## Installation

Python 3 and the packages listed in `requirements.txt` are required.

Clone the repository and navigate to its directory:

```bash
git clone https://github.com/rj3-1/NetMHCpan-IEDB-Benchmark.git
cd NetMHCpan-IEDB-Benchmark
```

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

## Running the analysis

Ensure the processed dataset is available at `data/processed/final_dataset.csv` and that `CSV_PATH` in `scripts/analysis.py` points to this location.

Run the script from the repository root:

```bash
python scripts/analysis.py
```

The script prints the statistical evaluation results and saves the generated figures in the `results/` directory.

## Results

The analysis produces quantitative measures of agreement between predicted and experimental binding affinities, classification performance at two IC50 thresholds, and differences between evidence groups.

Numerical findings and their interpretation should be reported after running the analysis on the final dataset.

## Limitations

- The analysis is limited to the peptides, experimental measurements, and HLA allele represented in the dataset.
- Experimental measurements may be missing for some peptides, limiting the subset available for direct comparison.
- Predicted binding affinity is not equivalent to experimental validation of T-cell immunogenicity.
- Results depend on the composition and quality of the underlying dataset.

## References

- **NetMHCpan:** Tool for predicting peptide–MHC binding.
- **IEDB (Immune Epitope Database):** Database of experimentally characterised immune epitopes and related assays.

Refer to the original NetMHCpan publication and official IEDB resources for tool and data documentation.
