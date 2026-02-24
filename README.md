# Ultrasonic Mid-Air Haptic SEP Analysis

Analysis code and stimulus control software for a master's thesis investigating the effect of spatiotemporal modulation (STM) drawing speed on somatosensory evoked potentials (SEPs) elicited by ultrasonic mid-air haptic stimulation.

## Overview

This repository contains three components:

| File | Language | Description |
|------|----------|-------------|
| `sep_processing.m` | MATLAB | EEG signal processing pipeline: trigger detection, filtering, epoching, artifact rejection, balanced trial selection, WPSS analysis, and visualisation |
| `sep_statistics.m` | MATLAB | Statistical analysis: repeated-measures ANOVA, post-hoc comparisons, trend analysis, normality testing, power analysis, SNR, and split-half reliability |
| `Program.cs` | C# | Stimulus control application for the Ultraleap Stratos device. Generates STM circle patterns at five speeds with audio trigger synchronisation |

## Experimental Setup

- **Stimulation:** Ultrasonic mid-air haptic stimulation (Ultraleap Stratos Explore) delivering STM circle patterns to the left palm
- **Conditions:** Five drawing speeds — 1.0, 3.0, 5.0, 7.0, 10.0 m/s
- **Recording:** 6-channel EEG (FC3, C3, CP3, FC4, C4, CP4) with Cz re-referencing
- **Trigger synchronisation:** Audio trigger captured by microphone channel, corrected for 55 ms acoustic delay
- **Epochs:** 100 balanced trials per speed per participant after artifact rejection (±20 µV threshold)

## Processing Pipeline

```
Raw EEG → Bandpass filter (2–32 Hz, FIR order 1000)
        → Trigger detection (MIC-first with REF rank matching)
        → Microphone delay correction (55 ms)
        → Re-referencing (Cz)
        → Epoching (−200 to +800 ms)
        → Baseline correction (−200 to 0 ms)
        → Artifact rejection (±20 µV)
        → Balanced selection (100 epochs/speed/participant)
        → Component extraction (P1, N1, P2)
        → WPSS analysis (Morse wavelet, 5–15 Hz band average)
```

## Requirements

### MATLAB scripts
- MATLAB R2020b or later
- Signal Processing Toolbox (for `fir1`, `filtfilt`)
- Wavelet Toolbox (for `cwtfilterbank`)
- Statistics and Machine Learning Toolbox (for `fitrm`, `ranova`)

### Stimulus control (C#)
- .NET Framework 4.7.2+
- [Ultraleap Haptics SDK](https://www.ultraleap.com/developers/)

## Usage

### 1. Processing

Place raw data in `data/mat_files/` and trial CSVs in `data/csv_files/`, then run:

```matlab
sep_processing
```

This produces:
- `GROUP_SEP_Data_Cz.mat` — all epochs before balancing
- `BALANCED_SEP_Data_100eps_Cz.mat` — balanced dataset
- `PEAK_ANALYSIS_Results_Cz.mat` — component amplitudes and latencies
- `WPSS_Results_Cz.mat` — wavelet phase synchronisation stability
- Visualisation PNGs (grand averages, heatmaps, WPSS spectrograms)

### 2. Statistics

Requires the balanced dataset from step 1:

```matlab
sep_statistics
```

Outputs:
- `MULTI_LEAN_Analysis_Cz.mat` — full statistical results
- `MULTI_LEAN_Results_Cz.csv` — per-participant component measures
- `MULTI_LEAN_ThesisReport_Cz.txt` — formatted results text
- 5 publication-ready figures

### 3. Stimulus control

Build and run the C# project with an Ultraleap device connected. The program measures the device sample rate, then presents the randomised trial sequence with audio triggers for EEG synchronisation.

## Data Availability

Raw EEG data are not included in this repository due to ethics and privacy considerations.

## Citation

If you use this code, please cite the associated thesis:

> L. Wolf, *Effect of Spatiotemporal Modulation Drawing Speed on Somatosensory Evoked Potentials Elicited by Ultrasonic Mid-Air Haptic Stimulation*, Master's Thesis, HTW Saar, 2026.

## License

MIT
