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
- Canonicalize monopolar labels  
- Construct bipolar montages  
- Apply channel-level quality control  
- Merge daily pain ratings  

### Spectral Analysis
- Compute power spectral density (PSD) using Welch’s method  
- Extract band-limited power (delta–gamma)  
- Compute peak metrics and slow index  
- Build epoch-level long-format tables  

### Aperiodic (1/f) Analysis
- Fit and remove aperiodic spectral slope  
- Extract residual alpha and beta peaks  

### Connectivity Analysis
- Compute coherence (COH)  
- Compute phase-locking value (PLV)  
- Generate connectivity long-format tables  

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

```matlab
% Preprocessing
load_raw_eeg_from_setfiles(...)
build_bipolar_montage(...)
run_bipolar_qc(...)

% Spectral
compute_psd_metrics(...)

% Connectivity
build_connectivity_longtable(...)

% Statistics
run_psd_lme_models(...)
run_connectivity_lme_models(...)


## Data Availability

De-identified EEG data and associated metadata are available through DABI:

https://dabi.loni.usc.edu/projects/3AXFYGFQQSCA

## Code Availability

All analysis code is available at:

https://github.com/UVA-LIU/UVA-LIU-uva-insular-iEEG-neuropathic-pain-analysis

Notes
This repository contains the final analysis pipeline only
Intermediate scripts are not included
All code is designed to be modular and reproducible
