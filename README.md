# fMRI HRF Analytics Pipeline

A Python-based pipeline for analyzing hemodynamic response function (HRF) data from fMRI retinotopic mapping studies. Migrated from MATLAB to Python with enhanced features.

## Overview

This tool processes MATLAB `.mat` files containing BOLD signal data, fits 7-parameter gamma function models to HRF curves, and extracts statistical metrics for visual cortex regions (V1, V2, V3, V4) and overall voxel activity.

## Features

- Weighted ROI averaging across hemispheres and sub-regions
- 7-parameter gamma function optimization with penalty-based constraints
- High-precision metric extraction at 0.01s resolution
  - Peak amplitude and time-to-peak
  - Full Width at Half Maximum (FWHM)
  - Undershoot trough and time-to-trough
- Batch processing with progress tracking
- Dual output format (MATLAB .mat and CSV)
- Visualization with R-squared goodness-of-fit indicators

## Installation

### Requirements

- Python 3.8 or higher

### Setup

1. Clone or download this repository

2. Install dependencies:

```bash
pip install -r requirements.txt
```

The required packages are:

- numpy
- scipy
- matplotlib
- pandas
- tqdm

## Usage

### Basic Command

```bash
python hrf_pipeline.py --input-dir <path> --output-dir <path>
```

### Arguments

- `--input-dir`: Directory containing `SUB*_hrf.mat` files (default: current directory)
- `--output-dir`: Directory for output files (default: `./output`)
- `--subjects`: Optional list of subject IDs (e.g., `SUB701 SUB702`)
- `--verbose`: Enable optimizer progress output

### Example

Process all subjects in the current directory:

```bash
python hrf_pipeline.py --input-dir . --output-dir ./output
```

Process specific subjects with verbose output:

```bash
python hrf_pipeline.py --subjects SUB701 SUB720 --verbose
```

## Input Format

The pipeline expects MATLAB `.mat` files with the following structure:

```
hrf_ROI
  ├── data (7x2x13): BOLD signals [ROIs × Hemispheres × Timepoints]
  ├── cnt (7x2): Voxel counts per sub-region
  └── all (1x13): Pre-calculated global mean signal
```

Time vector: 13 points at 2-second intervals (0, 2, 4, ..., 24s)

## Output Files

### 1. Parameters

`output/params/SUBxxx_params.mat`: MATLAB-compatible structure with fitted parameters

`output/params/SUBxxx_params.csv`: CSV table with columns:

- subject, roi, p1, q1, p2, q2, a1, a2, c

### 2. Statistics

`output/stats.tsv`: Tab-separated table with metrics:

- subject, roi, peak, t_peak, fwhm, trough, t_trough

### 3. Plots

`output/plots/SUBxxx_hrf_fit.png`: Multi-line plots showing:

- Raw data points (circles)
- Fitted gamma curves (lines)
- R-squared values in legend
- V1, V2, V3, V4, and All regions

## Model Details

### HRF Gamma Function

```
y(t) = a1 * Gamma(t; p1, q1) - a2 * Gamma(t; p2, q2) + c
```

Where:

- `p1, q1`: Shape and scale for primary response
- `p2, q2`: Shape and scale for undershoot
- `a1, a2`: Amplitude scaling factors
- `c`: Baseline offset

### Optimization Constraints

- Parameter bounds with penalty function
- Minimization of sum of squared errors (SSE)
- Initial values: `[6, 1, 12, 1.2, 30, 15, baseline]`

### ROI Averaging

- V1: Weighted average of rows 1-2 (both hemispheres)
- V2: Weighted average of rows 3-4
- V3: Weighted average of rows 5-6
- V4: Weighted average of row 7
- All: Pre-computed global mean

## Progress Tracking

The pipeline displays real-time progress during batch processing:

```
Processing subjects: 50%|██████████          | 5/10 [00:45<00:45, 9.1s/it]
```

## License

This code is provided for research purposes. Please cite appropriately if used in publications.

## Contact

For questions or issues, please contact the research team.
