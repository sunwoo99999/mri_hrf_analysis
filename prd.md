# [Final] PRD for fMRI HRF Analytics

## MATLAB-to-Python Migration

### 1. Data Input & Structure (Input Inference)

- **Format:** MATLAB `.mat` files (`SUBxxx_hrf.mat`).
- **Structure (`hrf_ROI`):**
- **data (7x2x13):** BOLD signals (7 ROIs, 2 Hemispheres, 13 Timepoints).
- **cnt (7x2):** Voxel counts per sub-region (used as weights for averaging).
- **all (1x13):** Pre-calculated global mean signal.

- **Time Vector:** seconds (13 points, 2s interval).

---

### 2. Core Research Logic (MATLAB Requirements)

- **Weighted ROI Averaging:**
- Group sub-regions into V1 (rows 1-2), V2 (rows 3-4), and V3 (rows 5-6).
- Apply weights:

- **7-Parameter Optimization:**
- **Model:** (where is baseline offset).
- **Objective Function:** Minimize Sum of Squared Errors (SSE) + Penalty.
- **Constraints:** Apply a penalty if parameters exceed defined bounds ().

---

### 3. Statistical Metric Extraction (0.01s Resolution)

Interpolate the fitted curve at a **0.01s interval** for high-precision measurement.

| Metric                | Description                                   | Constraints |
| --------------------- | --------------------------------------------- | ----------- |
| **Peak & T-peak**     | Maximum amplitude and its time.               | Ignore if   |
| **FWHM**              | Full Width at Half Maximum.                   | -           |
| **Trough & T-trough** | Minimum amplitude of undershoot and its time. | Ignore if   |

---

### 4. Batch Processing & Automation

- **Target Cohorts:** Automated iteration over `sub_Young`, `sub_Old`, and `sub_T2DM` lists.
- **Output Management:**
- Save optimized parameters as `SUBxxx_params.mat` in respective directories.
- Export consolidated statistics (Peak, FWHM, etc.) in a **tab-separated format**.

- **Visualization:**
- Generate multi-line plots for V1, V2, V3, and All.
- Apply `gem_colors` mapping and include legends/axes labels.

---

### 5. Technical Stack

- **Language:** Python 3.x
- **Core Libraries:**
- `scipy.io`: For loading `.mat` files.
- `numpy`: For matrix operations and math.
- `scipy.optimize`: For `fmin` optimization.
- `scipy.signal`: For peak detection.
- `matplotlib`: For plotting and visualization.

---
