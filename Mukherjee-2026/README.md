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

---

## 🧲 MRI Acquisition & Reconstruction

### Acquisition Method: RARE_VTR

The MRI data are acquired using the **RARE with Variable Repetition Times (RARE_VTR)** sequence on a preclinical Bruker MRI scanner. Each scan consists of multiple volumes acquired at different repetition times (TRs), enabling voxel-wise quantitative T1 mapping. Two TR vectors are supported:

- **6-volume VTR**: TRs = [8000, 2400, 1480, 940, 650, 418] msec
- **7-volume VTR**: TRs = [8000, 3600, 2400, 1480, 940, 650, 501.1] msec

Each animal contributes two longitudinal scan sessions (pre- and post-injection/treatment), acquired at two timepoints (Day 13 and Day 26 post-injection). Three cohorts (cohort_0, cohort_1, cohort_2) are included in the full dataset.

---

## 🔁 MRI Processing Pipeline

### 1. `mass_recon.py` — Batch Bruker-to-NIfTI Reconstruction

Converts raw Bruker PVDatasets to NIfTI format for all subjects in a directory, using the `bruker2nifti` package. Scan directories are automatically renamed according to their acquisition method.

**Usage:**
```bash
conda activate <your_env>
pip install bruker2nifti
python mass_recon.py -i <path/to/cohort_dir> -o <reconstruction_folder_name>
```

**Key functions:**
- `reconstruct_study()` — Instantiates a `Bruker2Nifti` object and converts all scans within a study directory
- `rename_scan_dirs()` — Renames reconstructed scan folders using the acquisition method string parsed from `acquisition_method.txt`

**Dependencies:** `bruker2nifti`, `shutil`, `os`

---

### 2. `rarevtr_pipeline.py` — Per-Scan Preprocessing, T1 Fitting, and Atlas Registration

Applies the full image processing pipeline to each RARE_VTR NIfTI file. The pipeline is run automatically over all matching files discovered via directory walk, or manually via explicit cohort path lists (`cohort_0_list`, `cohort_1_list`, `cohort_2_list`).

**Processing steps (in order):**

| Step | Function | Method |
|------|----------|--------|
| Motion correction | `preprocessing()` | Rigid-body registration of each TR volume to the first TR (mutual information metric, 3-level multiresolution); `dipy` |
| Denoising | `preprocessing()` | MP-PCA patch-based denoising (`dipy.denoise.localpca.mppca`) |
| Gibbs ringing removal | `preprocessing()` | Slice-wise Gibbs suppression (`dipy.denoise.gibbs.gibbs_removal`) |
| Skull stripping | `skull_stripping()` | SHERM rodent brain segmentation tool (called via `matlab.engine`) |
| Bias field correction | `bias_field_correcting()` | N3 bias field correction per TR volume (`ANTs`) |
| T1 map fitting | `t1map_fitting()` | Calls `rarevtr_t1fit.m` via `matlab.engine` (see below) |
| Atlas registration | `inverse_reg()` | Affine registration of masked T1-weighted image to DSURQE atlas (FSL FLIRT, 12 DOF, mutual information); inverse transform applied to atlas labels |

**Outputs per scan directory:**
- `corrected_data.nii.gz` — Motion-corrected, denoised, Gibbs-corrected 4D NIfTI
- `binary_mask.nii.gz` — Brain mask from SHERM
- `bias_corrected_data.nii.gz` — Bias-corrected 4D NIfTI
- `masked_first_vtr.nii.gz` — Masked first TR volume (used as registration target)
- `t1_map_corrected.nii.gz` — Voxel-wise T1 map (in msec)
- `labels_inv.nii.gz` — DSURQE atlas labels warped to native space

**Dependencies:** `dipy`, `ants`, `nipype` (FSL interface), `matlab.engine`

---

### 3. `rarevtr_t1fit.m` — Voxel-wise T1 Map Fitting (MATLAB)

Fits a mono-exponential signal recovery model to the multi-TR RARE_VTR data at every voxel within the brain mask, using non-linear least squares (`lsqcurvefit`).

**Signal model:**
```
S(TR) = M0 * (1 - exp(-TR / T1))
```

**Key parameters:**
- Starting guess: `[max(signal_intensity), 1800]` msec
- Bounds: `[0, 0]` to `[Inf, 10000]` msec
- Optimization: `lsqcurvefit` with `TolFun = 1e-3`, max 100 iterations (parallelized via `parfor`)

**Outputs:**
- `t1_map_corrected.nii.gz` — Full T1 map saved in native space
- `masked_first_vtr.nii.gz` — Masked first TR image (also used as FLIRT input)
- Return value: vector of T1 values within the brain mask

**Dependencies:** MATLAB Optimization Toolbox, `niftitools` (`load_untouch_nii` / `save_untouch_nii`)

---

### 4. `rarevtr_roi_analysis.m` — ROI-Level T1 Histogram Analysis

Applies the full per-subject pipeline across all cohorts and computes T1 histogram distributions and Earth Mover's Distance (EMD) metrics for each ROI, comparing pre- and post-injection timepoints.

**ROIs** (defined by DSURQE atlas label indices):

| ROI | Atlas Labels |
|-----|-------------|
| Thalamus | 4, 204 |
| Striatum | 7, 17 |
| Primary Somatosensory Cortex (SSC) | 111, 280 |
| Dentate Gyrus | 326–331 |

**Per-subject workflow:**
1. Locates pre- and post-injection scan directories from the cohort path list
2. Calls `segment_t1_map()` to extract T1 histograms per ROI (bin width = 25 msec; range = 1500–3000 msec)
3. Bootstraps histogram densities (1000 iterations) for robust estimation
4. Computes **EMD** between pre- and post-injection T1 distributions for each ROI
5. Parses subject metadata (genotype, infection status, timepoint, sex) from Bruker `subject` files
6. Generates per-subject summary figures and saves them to the study directory
7. Aggregates all subject data into `t1_RAREVTR_AllData.mat`

**Subject metadata parsed:**
- `Genotype`: WT or KO (from subject name string)
- `Infected`: true/false (non-infected, control, or PTX-treated → false)
- `Timepoint`: 1 (Day 13), 2 (Day 26), 0 (control), -1 (PTX-treated)
- `Sex`: M or F (from Bruker subject file)
- `Cohort`: cohort index (0, 1, 2)

**Output:** `t1_RAREVTR_AllData.mat` — structure with per-subject EMD values and metadata for all cohorts

**Dependencies:** MATLAB, `niftitools`, custom `emd()` function

---

### 5. `rarevtr_roi_analysis_ptx.m` — ROI Analysis Variant for PTX Controls

Variant of `rarevtr_roi_analysis.m` adapted for PTX-treated control subjects (timepoint = −1). Follows the same pipeline but handles the specific naming conventions and directory structure of the PTX cohort.

---

## 📊 Statistical Analysis

All statistical analyses are performed on the compiled EMD values extracted from `t1_RAREVTR_AllData.mat` (or its CSV export `pd_mouse_data_emd_only.csv`). The four outcome variables are EMD values for each ROI: **Striatum**, **Thalamus**, **Primary SSC**, and **Dentate Gyrus**.

The experimental design is a 2×2 factorial with two between-subject factors: **Genotype** (WT vs. KO) and **Infection status** (non-infected vs. infected), analyzed separately at each longitudinal timepoint.

---

### 6. `stat_analysis_all_cohorts.m` — MATLAB Statistical Analysis

Performs group-level statistical comparisons using the compiled dataset, separately for Timepoint 1 and Timepoint 2.

**Analyses performed:**
- **Violin plots** of EMD distributions per group (WT-NI, WT-IF, KO-NI, KO-IF) using the `violinplot` function
- **One-way groupwise ANOVA** followed by **multiple comparisons** (`multcompare`) across the four groups
- **Two-way interaction ANOVA** (`anovan` with `'model','interaction'`) with Genotype and Infection as factors

**Input:** `t1_RAREVTR_AllData.mat`

**Dependencies:** MATLAB Statistics and Machine Learning Toolbox, `violinplot` (third-party)

---

### 7. `Script.py` — Python Statistical Analysis

Replicates and extends the MATLAB analysis in Python, with additional cross-timepoint comparisons and Bonferroni correction.

**Analyses performed:**

- **Two-way ANOVA** (Type II SS) for each ROI × timepoint combination, with Genotype × Infection interaction term (`statsmodels.formula.api.ols`)
- **Tukey HSD post-hoc tests** for significant main effects and interactions (`statsmodels.stats.multicomp.pairwise_tukeyhsd`)
- **Bonferroni correction** applied across all post-hoc p-values (`statsmodels.stats.multitest.multipletests`)
- **Welch's t-test** comparing Timepoint 1 vs. Timepoint 2 for each ROI (pooled across genotype/infection groups)
- **Genotype-stratified interaction ANOVA**: for KO and WT groups separately, tests Infection × Timepoint interaction for each ROI

**Input:** `pd_mouse_data_emd_only.csv`

**Required columns:** `Genotype`, `Infected`, `Timepoint`, `EMDstriatum`, `EMDthalamus`, `EMDpssc`, `EMDdg`

**Dependencies:** `pandas`, `numpy`, `statsmodels`, `scipy`

Install with:
```bash
pip install pandas numpy statsmodels scipy
```

---

## 🗂 Full Pipeline Summary

```
Raw Bruker Data
      │
      ▼
mass_recon.py           → NIfTI conversion + scan directory renaming
      │
      ▼
rarevtr_pipeline.py     → Motion correction → Denoising → Gibbs removal
                           → Skull stripping (SHERM) → Bias correction (N3)
                           → T1 fitting (rarevtr_t1fit.m)
                           → Atlas registration (FSL FLIRT + inverse warp)
      │
      ▼
rarevtr_roi_analysis.m  → ROI segmentation (DSURQE atlas)
                           → T1 histogram extraction + bootstrapping
                           → EMD computation (pre vs. post)
                           → Subject metadata parsing
      │
      ▼
t1_RAREVTR_AllData.mat / pd_mouse_data_emd_only.csv
      │
      ├──▶ stat_analysis_all_cohorts.m   (MATLAB: violin plots, ANOVA)
      └──▶ Script.py                     (Python: two-way ANOVA, Tukey, Bonferroni)
```
