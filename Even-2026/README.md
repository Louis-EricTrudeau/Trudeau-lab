# 🧪 Image Analysis - Data Compilation Toolkit #

This repository contains NIS-Elements AR `.ga3` analysis recipes, developed for the image-analysis workflows associated with the following study:

***Repeated exposure to viral and bacterial danger signals in Parkin deficient mice induces inflammation but fails to trigger parkinsonism* — Even et al., 2026**

## 📚 Project Overview

The `GA3` directory contains **two separate image-analysis workflows** (organized into separate directories):

1. **"IBA1/GFAP" analysis**

   * Analysis of
   *  a) IBA1 and GFAP immunofluorescence volumic intensity
   *  b) IBA1 and GFAP immunofluorescence area and volume
   *  c) Count of Iba1 cells
   *  d) Count of GFAP cells
   * Contains two separate `.ga3` recipes because they correspond to different biological targets (**mesencephalon** and **striatum**), and so different measurements and image-analysis procedures.

2. **"TH-DAT" analysis**

   * Quantitative analysis of TH/DAT-related signal in the **striatum**.

---

# Software

The recipes were saved from:

* **NIS-Elements AR**
* Application version: **6.10.02**
* Application build: **2031**

The `.ga3` files are native NIS-Elements General Analysis Graph recipes and should therefore be opened with a compatible NIS-Elements installation.

> **Important:** Compatibility may depend on the NIS-Elements version and on the analysis modules installed on the workstation.

---

# 1. IBA1/GFAP analysis

## Overview

The **IBA1/GFAP workflow** is dedicated to quantitative analysis of IBA1 and GFAP immunofluorescence and to count Iba1 and GFAP cells.

## Repository

The IBA1/GFAP recipes are located in:

https://github.com/Louis-EricTrudeau/Trudeau-lab/tree/main/Even-2026/GA3/IBA1_GFAP


It contains two separate `.ga3` recipes corresponding to two anatomical regions:

```text
GA3/
└── IBA1_GFAP/
    ├── Mes_AI_Iba1_GFAP.ga3
    └── Stri_AI_Iba1_GFAP.ga3
```

| Recipe                  | Region        | Analysis  |
| ----------------------- | ------------- | --------- |
| `Mes_AI_Iba1_GFAP.ga3`  | Mesencephalon | IBA1/GFAP |
| `Stri_AI_Iba1_GFAP.ga3` | Striatum      | IBA1/GFAP |


---

## AI segmentation

The recipe use an **AI-based object segmentation** step.

The AI model is not provided as part of the repository because AI segmentation should be performed using a model **trained for the specific dataset and analysis requirements**.

The AI training procedure is described in the following protocol:

[AI training protocol — Protocols.io](https://dx.doi.org/10.17504/protocols.io.x54v9qn5zl3e/v1)

The use of AI segmentation is therefore not tied to a specific pre-existing model. When reproducing the analysis, users should train an appropriate AI model following the protocol and apply it to their own datasets.


## 1.1 `Mes_AI_Iba1_GFAP.ga3`

This recipe is intended for **IBA1/GFAP analysis in the mesencephalon**.
In this recipe, we select the TH volume/area and the analysis is done inside this volume/area

The analysis workflow includes operations such as:

* channel/image selection;
* intensity equalization;
* threshold-based segmentation;
* morphological object separation;
* AI-based object segmentation when required;
* object counting;
* object intensity measurements;
* object-area measurements;
* aggregation of measurements;
* export to a `.csv` file.

### Main measurements

The analysis contains measurements corresponding to:

* IBA1 object count (using AI);
* IBA1 mean intensity per object;
* IBA1 summed object area and volume;
* GFAP mean intensity per object;
* GFAP summed object area and volume;
* DAPI/object-count related measurements;
* TH intensity;
* TH area and volume.

---

## 1.2 `Stri_AI_Iba1_GFAP.ga3`

This recipe is intended for **IBA1/GFAP analysis in the striatum**.

It is a separate `.ga3` recipe from the mesencephalon analysis because it does not include TH selection for the analysis, but follows the same general IBA1/GFAP analysis strategy.

The workflow includes:

1. Input image/channel selection.
2. Intensity equalization.
3. Threshold-based segmentation.
4. Object separation.
5. AI-based segmentation when required.
6. Object measurements.
7. Aggregation of measurements.
8. Export to a `.csv` file.


### Main measurements

The analysis contains measurements for:

* DAPI/object counts;
* IBA1 AI object count;
* IBA1 mean intensity per object;
* IBA1 summed object area;
* GFAP mean intensity per object;
* GFAP summed object area;
* object intensity;
* object area;
* aggregated object measurements.

The final step of the analysis is **export to a `.csv` file**.

---

## 1.3 IBA1/GFAP analysis workflow

The general IBA1/GFAP workflow can be summarized as:

```text
Raw ND2 image
      │
      ▼
Channel selection
      │
      ▼
Intensity equalization / preprocessing
      │
      ▼
Segmentation
      │
      ├── Conventional segmentation
      │
      └── AI object segmentation, when required
      │
      ▼
Object measurements
      │
      ├── Object count
      ├── Mean intensity
      └── Object area
      │
      ▼
Aggregation
      │
      ▼
Export to a .csv file
```


# 2. TH-DAT analysis

## Overview

The **TH-DAT workflow** is a separate analysis dedicated to quantitative measurement of TH/DAT-related signal in the striatum.

## Repository

The TH-DAT recipe is located in:

https://github.com/Louis-EricTrudeau/Trudeau-lab/tree/main/Even-2026/GA3/TH_DAT

It contains one `.ga3` recipe:

```text
GA3/
└── TH_DAT/
    └── TH-DAT_intensity_area_striatum.ga3


```

## 2.1 `TH-DAT_intensity_area_striatum.ga3`

This recipe is designed to quantify **TH/DAT-related signal intensity and segmented area in the striatum**.

### General workflow

The analysis graph contains the following major operations:

1. **Rolling-ball preprocessing**
2. **Threshold-based segmentation**
3. **Object selection/separation**
4. **Object smoothing**
5. **Object/field intensity measurement**
6. **Object/field area measurement**
7. **Column calculations**
8. **Joining and modifying measurement tables**
9. **Calculation of percentage area**
10. **Final column selection/compaction**
11. **Export to a `.csv` file**

### Main measurements

The recipe contains measurements including:

* Fiber_intensity
* Whole_field_intensity
* Fiber_area
* Whole_field_area
* Fiber_area_in_pct

The percentage-area calculation is based on:

```text
Fiber area / Whole-field area × 100
```

This provides an estimate of the proportion of the analyzed field occupied by the segmented signal.

---

## 2.2 TH-DAT analysis workflow

The general TH-DAT workflow can be summarized as:

```text
Raw ND2 image
      │
      ▼
Channel selection
      │
      ▼
Rolling-ball preprocessing
      │
      ▼
Threshold segmentation
      │
      ▼
Object selection / smoothing
      │
      ▼
Intensity and area measurements
      │
      ▼
Percentage-area calculation
      │
      ▼
Export to a .csv file
```

---

# 3. Input images and reproduction

All `.ga3` recipes should be used with the appropriate `.nd2` images and a compatible NIS-Elements AR installation.
Input image paths and other workstation-specific settings may need to be updated when reproducing the analyses on another computer.

When reproducing an analysis:

1. Open the appropriate `.ga3` recipe in NIS-Elements.
2. Replace the input image with the appropriate `.nd2` dataset.
3. Verify that all required channels are present.
4. Verify the image calibration.
5. Verify the segmentation parameters.
6. For IBA1/GFAP analyses using AI segmentation, train an appropriate AI model following the dedicated protocol.
7. Verify the output settings.
8. Run the analysis and inspect the resulting measurements.

The appropriate recipe should always be selected according to the biological analysis:

* **IBA1/GFAP → `GA3/IBA1_GFAP/`**
* **TH-DAT → `GA3/TH_DAT/`**

---

# 4. Comparison of the two analysis workflows

The `GA3` directory therefore contains **two distinct and independent image-analysis workflows**.

|                    | IBA1/GFAP                                                              | TH-DAT                                                                      |
| ------------------ | ---------------------------------------------------------------------- | --------------------------------------------------------------------------- |
| Directory          | `GA3/IBA1_GFAP/`                                                       | `GA3/TH_DAT/`                                                               |
| Biological targets | IBA1 / GFAP                                                            | TH / DAT                                                                    |
| Anatomical regions | Mesencephalon / Striatum                                               | Striatum                                                                    |
| `.ga3` recipes     | 2                                                                      | 1                                                                           |
| Main approach      | Segmentation + object measurements, with AI segmentation when required | Rolling-ball preprocessing + threshold segmentation                         |
| Main measurements  | Object count, mean intensity, object area                              | Signal intensity, field intensity, signal area, field area, percentage area |
| AI                 | Used when required; users should train an appropriate model            | Not used in this workflow                                                   |
| Output             | Quantitative IBA1/GFAP measurements                                    | Quantitative TH/DAT measurements                                            |

The two workflows should therefore be understood as **separate analysis pipelines**, even though both are implemented using NIS-Elements `.ga3` recipes.

---

# Summary

The `GA3` directory contains the image-analysis recipes used for the study:

***Repeated exposure to viral and bacterial danger signals in Parkin deficient mice induces inflammation but fails to trigger parkinsonism* — Even et al., 2026**

The recipes are divided into two independent workflows:

### IBA1/GFAP

Located in:

```text
GA3/IBA1_GFAP/
```

Contains:

* `Mes_AI_Iba1_GFAP.ga3` — mesencephalon IBA1/GFAP analysis
* `Stri_AI_Iba1_GFAP.ga3` — striatum IBA1/GFAP analysis

These recipes quantify IBA1/GFAP-related object measurements and can use AI-based segmentation when required. For AI segmentation, users should train an appropriate model for their own dataset following the dedicated training protocol.

### TH-DAT

Located in:

```text
GA3/TH_DAT/
```

Contains:

* `TH-DAT_intensity_area_striatum.ga3` — striatal TH/DAT intensity and area analysis

This recipe uses a separate image-analysis workflow based on rolling-ball preprocessing, threshold segmentation, intensity measurements and area measurements.

All `.ga3` recipes should be used with the appropriate `.nd2` images and a compatible NIS-Elements AR installation.
Input image paths and other workstation-specific settings may need to be updated when reproducing the analyses on another computer.





