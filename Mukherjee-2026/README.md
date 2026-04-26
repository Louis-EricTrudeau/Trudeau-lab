# 🧪 Image Analysis - Data Compilation Toolkit

This repository contains a set of Python scripts developed to automate and standardize the compilation of quantitative data from immunohistochemistry image analyses performed using ImageJ/Fiji. These analyses include soma count and morphology of Iba1-positive microglia, arborisation features, and protein expression intensity for TH, DAT, and tdTomato in striatal regions of the mouse brain.

The dataset processed by these scripts typically originates from `.csv` exports generated via custom ImageJ macros (also included separately), which extract per-cell or per-ROI measurements from brain slice images. These outputs are parsed, filtered, and aggregated across samples and brain regions to create tidy Excel spreadsheets for downstream statistical analysis.

---

## 📚 Project Overview

The goal of this pipeline is to facilitate high-throughput, reproducible quantification of:
- **Microglial morphology and distribution (Iba1)** in dorsal and ventral striatum
- **Fluorescence intensity and area coverage** of tyrosine hydroxylase (TH), dopamine transporter (DAT), and tdTomato in the same regions

Each dataset is associated with a brain sample (identified as "Brain #") and spatially registered to either the **dorsal** or **ventral striatum**. The analysis is structured to maintain regional specificity, biological replicates, and quantitative feature integrity.

**Important**: You must first run the provided ImageJ (.ijm) macros to analyze `.oif` image files and export `.csv` measurement files. These macros extract relevant features such as area, mean intensity, circularity, etc., from the images.

### Biological Context

This analysis framework was designed in the context of neurodegenerative research, particularly focusing on:
- **Microglial activation patterns**, assessed via Iba1-positive cell count and morphological metrics (e.g., roundness, perimeter, solidity)
- **Neuronal and dopaminergic system integrity**, assessed via intensity and %area of marker proteins such as TH and DAT

By averaging results across ROIs and across animals, the scripts produce output files that can be used for statistical comparisons between experimental groups (e.g., control vs treated, wild-type vs knockout, etc.)

---

## 🧩 Features

- **Automatic data aggregation**: Combines measurements from multiple CSVs into a unified DataFrame
- **Per-brain averaging**: Ensures each brain is treated as a biological replicate
- **Flexible CSV handling**: Tolerant to variations in file structure if naming conventions are maintained
- **Export to Excel**: Outputs `.xlsx` spreadsheets ready for statistical software or plotting libraries

---

## 🧪 Script Overview

### 1. `1. Iba1 count compilation.py`
- **Inputs**: `.csv` files from the `Iba1 soma/` folder
- **Measurements**:
  - Cell count (normalized by area)
  - Mean Area, Mean Intensity, Perimeter, Circularity, Roundness, Solidity
- **Region split**: Dorsal (D) and Ventral (V) striatum
- **Output**: `compilation-Iba1-soma.xlsx`

### 2. `2. Iba1 arborisation compilation.py`
- **Inputs**: `.csv` files from the `Iba1 arborisation/` folder
- **Measurements**:
  - Area and Intensity of microglial arbors
- **Output**: `compilation-Iba1-arborisation.xlsx`

### 3. `TH DAT compile.py`
- **Inputs**: `intensity.csv` files inside each brain folder
- **Measurements**:
  - Mean intensity and %area for:
    - Tyrosine Hydroxylase (TH)
    - Dopamine Transporter (DAT)
    - tdTomato reporter
- **Output**: `compilation-ec-gaussian-threshold.xlsx`

---

## 🛠 Installation & Dependencies

These scripts require Python 3.8+ and the following Python packages:

- `pandas`
- `numpy`

You can install dependencies with:

```bash


pip install pandas numpy
