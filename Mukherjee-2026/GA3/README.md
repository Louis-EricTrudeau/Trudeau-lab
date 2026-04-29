# NIS-Elements GA3 Macro Pipeline
### Dopaminergic Neuron Quantification and Microglial Analysis

**Software:** NIS-Elements AR 6.20.02  
**Format:** General Analysis 3 (`.ga3`)  
**Application:** Quantification of dopaminergic neurons in SNc and VTA, and microglial marker analysis (Iba1 / CD68), for Parkinson's disease research models

---

## Overview

This pipeline consists of seven GA3 macro files that together automate the full analysis workflow from raw multi-channel fluorescence images through AI-based segmentation to final regional cell counts. Steps are numbered to reflect execution order; steps 1 and 5 (acquisition and intermediate review, respectively) are performed outside this pipeline.

The two branches of the pipeline are:

- **Dopaminergic neuron branch (steps 2 to 8):** Exports, preprocesses, runs two rounds of AI segmentation, defines anatomical ROIs, and counts TH+ / DA neurons in the SNc and VTA.
- **Microglial analysis branch (standalone):** Runs EDF projection, AI segmentation, and threshold-based detection for Iba1 and CD68 markers.

---

## File Descriptions

### `2-Export_with_1_channel.ga3`
**Purpose:** Exports a single fluorescence channel from a multi-channel ND2 acquisition file into a standalone TIFF for downstream processing.

**What it does:**
- Opens the source ND2 file from the acquisition folder
- Isolates the target channel (typically the TH or DAT channel marking dopaminergic neurons)
- Exports a single-channel TIFF to a designated output directory

**When to run:** After image acquisition (step 1). Must be run before any preprocessing or AI segmentation steps.

**Input:** Multi-channel `.nd2` file  
**Output:** Single-channel `.tif` file

---

### `3-Pre-AI-1_macro_on_decon.ga3`
**Purpose:** Preprocesses deconvolved images in preparation for the first round of AI segmentation. Uses a **Bright Spot Detection** node to enhance and isolate neuronal soma signals.

**What it does:**
- Opens the deconvolved TIFF (output from a deconvolution step)
- Applies a composed bright spot detection filter (`SegComposedBrightSpotDetectionImpl`) to enhance punctate or soma-like structures
- Saves the processed image as a new TIFF ready for AI-1 input

**When to run:** After deconvolution of the exported channel TIFF. Run before step 4.

**Input:** Deconvolved single-channel `.tif`  
**Output:** Pre-processed `.tif` optimized for AI segmentation

**Note:** The bright spot detection node is tuned for the morphology of dopaminergic neuron cell bodies. If applied to a different cell type or marker, the detection parameters may need adjustment inside the GA3 graph editor.

---

### `4-Run_AI-1_on_samples_TIFF.ga3`
**Purpose:** Runs the **first AI segmentation pass** (AI-1) on preprocessed TIFF files using NIS-Elements' deep learning segmentation engine.

**What it does:**
- Loads each preprocessed TIFF from the input directory
- Applies a trained AI model via the `ComposedSegmentObjectsAiImpl` deep learning node to segment candidate objects (neuronal soma)
- Outputs segmentation masks or labeled binary images for each sample

**When to run:** After step 3. Requires a trained NIS-Elements AI model to be loaded and associated with the `ComposedSegmentObjectsAiImpl` node.

**Input:** Pre-processed `.tif` from step 3  
**Output:** AI-1 segmentation mask `.tif`

**Prerequisite:** The AI model file (`.h5` or NIS model package) must be accessible and correctly referenced within the GA3 node. If the model path has changed, re-link it inside the deep learning node settings.

---

### `6-Run_AI-2_on_TIFF.ga3`
**Purpose:** Runs a **second AI segmentation pass** (AI-2) to refine or complement the detections from AI-1, using a separate trained model or detection strategy.

**What it does:**
- Loads TIFF files (typically the output of AI-1 or the original pre-processed images)
- Applies a second `ComposedSegmentObjectsAiImpl` deep learning node with a different model or parameter set
- Produces refined segmentation masks for the counting step

**When to run:** After step 4 (and any intermediate review or correction at step 5). Run before step 7.

**Input:** TIFF files from step 4 output  
**Output:** AI-2 refined segmentation mask `.tif`

**Note:** Steps 4 and 6 use the same underlying segmentation node type but are intended to apply different models or handle different object classes (e.g., AI-1 for soma detection, AI-2 for exclusion of false positives or detection of a second marker). Confirm model assignment in each file separately before running.

---

### `7-Draw_SNc_and_VTA.ga3`
**Purpose:** Guides the user through **manually drawing anatomical ROIs** for the substantia nigra pars compacta (SNc) and ventral tegmental area (VTA) on each brain section.

**What it does:**
- Opens the relevant image (typically the TH-channel TIFF or a merged reference image)
- Prompts the user to draw or load ROI boundaries corresponding to the SNc and VTA
- Saves the ROI definitions (coordinates or binary masks) for use in step 8

**When to run:** After AI segmentation is complete (step 6). Must be run before step 8 since the counting macro depends on the ROI definitions produced here.

**Input:** Reference image TIFF (TH channel or overlay)  
**Output:** ROI file(s) defining SNc and VTA boundaries per section

**Note:** ROI placement should follow a standardized atlas reference (e.g., Paxinos and Franklin Mouse Brain Atlas) to ensure reproducibility across sections and animals. This step requires manual intervention — it cannot be run in fully automated batch mode.

---

### `8-Count_SNc_and_VTA_DA_neurons.ga3`
**Purpose:** **Counts dopaminergic (DA) neurons** within the SNc and VTA ROIs defined in step 7, using the AI segmentation masks from steps 4 and 6.

**What it does:**
- Loads the segmentation masks and the corresponding ROI definitions
- Restricts object detection to within each anatomical ROI
- Counts the number of segmented objects (TH+ neuronal soma) in the SNc and VTA separately
- Exports count data to a results table or `.csv` file

**When to run:** After step 7 (ROIs must exist). This is the final step of the dopaminergic neuron quantification branch.

**Input:** AI segmentation masks from steps 4/6 and ROI files from step 7  
**Output:** Cell count results table (per ROI, per section, per animal)

---

### `Iba1_and_CD68_-_only_EDF.ga3`
**Purpose:** **Standalone microglial analysis pipeline** for Iba1 (ionized calcium-binding adapter molecule 1, pan-microglial marker) and CD68 (lysosomal marker for activated microglia). Uses Extended Depth of Focus (EDF) projection combined with AI and threshold-based segmentation.

**What it does:**
- Applies EDF projection to collapse z-stack images into a single in-focus plane for each channel
- Runs AI-based segmentation (`ComposedSegmentObjectsAiImpl`) to detect microglial cells or processes
- Applies four independent threshold segmentation nodes (`SegComposedThresholdImpl`) for channel-specific intensity-based detection of Iba1 and CD68 signals
- Outputs segmented objects and/or area measurements for both markers

**When to run:** Independently of the dopaminergic neuron branch. Can be run on its own set of sections or the same tissue if both markers were acquired.

**Input:** Multi-channel z-stack `.nd2` or `.tif` files containing Iba1 and CD68 channels  
**Output:** Segmentation results and/or measurement tables for Iba1 and CD68

**Note:** The four threshold nodes suggest separate segmentation passes — likely one per channel and one or more for co-localization or morphological sub-classification (e.g., CD68+ area within Iba1+ cells as a proxy for microglial activation). Threshold values should be validated against a set of representative images before batch processing.

---

## Execution Order

```
Step 1   [Acquisition — external]
   ↓
Step 2   2-Export_with_1_channel.ga3
   ↓
          [Deconvolution — external]
   ↓
Step 3   3-Pre-AI-1_macro_on_decon.ga3
   ↓
Step 4   4-Run_AI-1_on_samples_TIFF.ga3
   ↓
Step 5   [Review / QC — manual, external]
   ↓
Step 6   6-Run_AI-2_on_TIFF.ga3
   ↓
Step 7   7-Draw_SNc_and_VTA.ga3       ← requires manual input
   ↓
Step 8   8-Count_SNc_and_VTA_DA_neurons.ga3


Iba1_and_CD68_-_only_EDF.ga3        ← independent branch, run separately
```

---

## Prerequisites

| Requirement | Details |
|---|---|
| NIS-Elements AR | Version 6.20.02 or later |
| GA3 module | General Analysis 3 license required |
| Deep learning module | Required for steps 3, 4, 6, and the Iba1/CD68 pipeline |
| AI model files | Must be linked inside each `ComposedSegmentObjectsAiImpl` node before running |
| Atlas reference | Paxinos & Franklin (or equivalent) for SNc/VTA ROI drawing in step 7 |

---

## Notes on AI Model Linking

Each file containing a `ComposedSegmentObjectsAiImpl` node (`4-Run_AI-1_on_samples_TIFF.ga3`, `6-Run_AI-2_on_TIFF.ga3`, `Iba1_and_CD68_-_only_EDF.ga3`) requires a pre-trained model to be associated with that node. If macros are transferred to a new workstation:

1. Open the `.ga3` file in NIS-Elements GA3 editor
2. Double-click the AI segmentation node
3. Re-link the model file using the model path field
4. Validate on a test image before batch processing

---

## Output Files

| Step | Output |
|---|---|
| Step 2 | Single-channel `.tif` per sample |
| Step 3 | Pre-processed `.tif` per sample |
| Step 4 | AI-1 segmentation mask `.tif` per sample |
| Step 6 | AI-2 segmentation mask `.tif` per sample |
| Step 7 | ROI definition files (SNc, VTA) per section |
| Step 8 | Cell count table `.csv` (SNc count, VTA count, per section/animal) |
| Iba1/CD68 | Segmentation masks and measurement table for Iba1 and CD68 |

---

## Contact

For questions about this pipeline, contact the lab member responsible for image analysis protocol development.