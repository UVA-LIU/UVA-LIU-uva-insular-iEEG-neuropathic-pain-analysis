# Insular iEEG Analysis for Chronic Neuropathic Pain

This repository contains MATLAB code for preprocessing, spectral analysis, connectivity analysis, and mixed-effects modeling of intracranial EEG (iEEG) data from patients with chronic neuropathic pain.

This repository contains the exact code used to generate all results reported in the manuscript.

**Liu et al. (2026)**  
Posterior insular oscillations encode pain intensity and unpleasantness in chronic neuropathic pain  
(under revision at *Brain Communications*)

---

## Overview

The analysis pipeline includes:

### Preprocessing
- Load EEGLAB `.set` files  
- Canonicalize monopolar channel labels  
- Construct bipolar montages  
- Apply channel-level quality control  
- Merge daily pain ratings  

### Spectral Analysis
- Compute power spectral density (PSD) using Welch’s method  
- Extract band-limited power (delta, theta, alpha, beta, gamma)  
- Compute peak metrics (peak frequency and power)  
- Compute slow index (theta/alpha)  
- Build epoch-level long-format tables  

### Aperiodic (1/f) Analysis
- Fit and remove aperiodic spectral slope  
- Extract residual alpha and beta peaks  

### Connectivity Analysis
- Compute coherence (COH)  
- Compute phase-locking value (PLV)  
- Generate pair-level long-format connectivity tables  

### Statistical Modeling
- Linear mixed-effects models  
- Trait (baseline) vs state (pain-related) effects  
- Control for longitudinal drift  
- False discovery rate (FDR) correction  

---

## Requirements

- MATLAB (R2023b or later recommended)  
- Signal Processing Toolbox  
- Statistics and Machine Learning Toolbox  
- EEGLAB  

---

## How to Run

### Step 1: Configure inputs

```matlab
rootPath = '/path/to/data';

subjectIDs = {'DBS001','DBS002','DBS003','DBS004','DBS005','DBS006'};

% Example subjectInfo structure
subjectInfo.DBS001 = {'Day_1','Day_2','Day_3'};
subjectInfo.DBS002 = {'Day_1','Day_2','Day_3'};

regionDef = [];  % optional

exclusionFile = fullfile(rootPath,'bipolar_channel_exclusions.csv');

opts = struct();
```

---

### Step 2: Run full pipeline

```matlab
run_all_analysis(rootPath, subjectIDs, subjectInfo, regionDef, exclusionFile, opts);
```

---

### Step 3: Run individual steps (optional)

```matlab
% Preprocessing
load_raw_eeg_from_setfiles(rootPath, subjectInfo, [], regionDef, opts);
build_bipolar_montage(rootPath, subjectIDs, regionDef);
run_bipolar_qc(rootPath, exclusionFile, subjectIDs);

% Spectral analysis
compute_psd_metrics(rootPath, subjectIDs);

% Connectivity analysis
build_connectivity_longtable(rootPath, subjectIDs);

% Statistical modeling
run_psd_lme_models([], rootPath);
run_connectivity_lme_models([], rootPath);
```

---

## Data Availability

De-identified EEG data are available from the corresponding author upon reasonable request. The dataset is currently being prepared for public release via the DABI repository.

---

## Code Availability

All analysis code used in this study is available at:

https://github.com/UVA-LIU/UVA-LIU-uva-insular-iEEG-neuropathic-pain-analysis

---

## Notes

- This repository contains the final analysis pipeline only  
- Intermediate scripts are not included  
- All code is designed to be modular and reproducible  

---

## Contact

Jeff Liu, PhD  
University of Virginia
