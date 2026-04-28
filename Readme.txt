# Insular iEEG Analysis for Chronic Neuropathic Pain

This repository contains MATLAB code for preprocessing, spectral analysis, connectivity analysis, and mixed-effects modeling of intracranial EEG (iEEG) data from patients with chronic neuropathic pain.

The code reproduces the analyses reported in:

> Liu et al., "Posterior insular oscillations encode pain intensity and unpleasantness in chronic neuropathic pain", *Brain* (under revision)

---

## Overview

The analysis pipeline includes:

### 1. Preprocessing
- Load EEGLAB `.set` files  
- Canonicalize monopolar channel labels  
- Construct bipolar montages  
- Apply channel-level quality control  
- Merge daily pain ratings  

### 2. Spectral Analysis
- Compute power spectral density (PSD) using Welch’s method  
- Extract band-limited power (delta, theta, alpha, beta, gamma)  
- Compute peak metrics (peak frequency and power)  
- Compute slow index (theta/alpha)  
- Build epoch-level long-format tables  

### 3. Aperiodic (1/f) Analysis
- Fit aperiodic spectral slope in log–log space  
- Subtract background activity  
- Extract residual alpha and beta peaks  
- Validate oscillatory contributions to pain-related effects  

### 4. Connectivity Analysis
- Compute coherence (COH)  
- Compute phase-locking value (PLV)  
- Generate pair-level long-format connectivity tables  

### 5. Statistical Modeling
- Linear mixed-effects models  
- Separate:
  - Trait (stable anatomical organization)  
  - State (pain-related variation)  
- Control for longitudinal drift  
- Apply false discovery rate (FDR) correction  

---

## Repository Structure
