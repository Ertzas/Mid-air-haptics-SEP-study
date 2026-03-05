# Spatiotemporal Modulation Speed Effects on Somatosensory Evoked Potentials

This repository contains the experimental control and analysis code for a Master's thesis investigating how the drawing speed of spatiotemporal modulation (STM) stimuli in mid-air ultrasonic haptics affects cortical processing, as measured by somatosensory evoked potentials (SEPs) and subjective preference.

**Five drawing speeds were tested:** 1.0, 3.0, 5.0, 7.0, and 10.0 m/s  
**Stimulation device:** Ultraleap Stratos (single focal point, STM via circular trajectory)  
**Stimulation site:** Left hand (palm)  
**EEG focus:** Contralateral electrodes C4, CP4, FC4

---

## Repository Structure

```
.
├── objective_program.cs          # EEG study – stimulus presentation (C#)
├── processing_data.m         # EEG – preprocessing, epoching, balanced selection, visualization
├── statistics.m                  # EEG – rm-ANOVA, post-hoc, trend analysis, SNR, split-half
│
├── intensity_script.cs           # Subjective study – paired comparison, intensity dimension (C#)
├── valence_script.cs             # Subjective study – paired comparison, valence dimension (C#)
├── paired_comparison_analysis.m  # Subjective study – Bradley-Terry model, win rates
│
└── data/
    ├── csv_files/                # EEG trial logs (one CSV per participant)
    ├── mat_files/                # Raw EEG recordings (one MAT per block)
    └── csv_files_subjective/     # Paired comparison response logs
```

---

## Study 1 – EEG / SEP

### Overview

Participants received mid-air ultrasonic haptic stimuli to the left palm while EEG was recorded. Each trial presented one of five speed conditions. SEPs were extracted from contralateral electrodes and analysed for P1, N1, and P2 components.

### Step 1 – Stimulus Presentation

**File:** `objective_program.cs`

Runs the experiment: generates a randomised trial sequence, presents stimuli via the Ultraleap device, records audio triggers through a microphone, and logs timing to CSV.

**Requirements:**
- .NET 6.0 or later
- [Ultraleap Haptics SDK](https://developer.leapmotion.com/releases) (referenced as `Ultraleap.Haptics`)
- A custom TCP trigger server (`TcpServer_And_Ultrahaptics`) — included in the original project build

**Key parameters (edit at top of file):**

| Parameter | Default | Description |
|---|---|---|
| `TrialsPerBlock` | 150 | Trials per block |
| `TotalBlocks` | 5 | Number of blocks |
| `ISI` | 2000 ms | Inter-stimulus interval |
| Speeds | 1, 3, 5, 7, 10 m/s | Randomised within each block |

**Output:**
- `P{ID}_trials.csv` — trial log with columns: `trial, participant, block, speed, isi_duration, stimulus_onset`

**Usage:**
```
dotnet run
# Enter participant ID when prompted (integer, e.g. 3 → P003)
```

---

### Step 2 – EEG Preprocessing & Analysis

**File:** `processing_data.m`

Loads raw EEG recordings and trial logs, detects triggers, epochs, rejects artifacts, applies balanced epoch selection, and produces grand average figures and WPSS (Wavelet Phase Synchronisation Stability) time-frequency maps.

**Required MATLAB Toolboxes:**
- Signal Processing Toolbox (`filtfilt`, `fir1`)
- Wavelet Toolbox (`cwtfilterbank`, `wt`)
- Statistics and Machine Learning Toolbox (`prctile`)

**Expected input files:**

| File | Location | Naming convention |
|---|---|---|
| Trial log | `data/csv_files/` | `P003_trials.csv` |
| Raw EEG (per block) | `data/mat_files/` | `s03_tr1_SS.mat`, `s03_tr2_SS.mat`, … |

**MAT file channel layout (row index):**

| Row | Signal |
|---|---|
| 2 | Reference trigger |
| 4 | Microphone trigger |
| 6 | C3 |
| 7 | C4 |
| 8 | CP3 |
| 9 | CP4 |
| 10 | FC3 |
| 11 | FC4 |
| 12 | Right earlobe |
| 13 | Left earlobe |
| 14 | Cz |

**Key configuration (edit at top of script):**

```matlab
PARTICIPANTS         = {'P003', 'P004', ...};  % participant IDs to include
REFERENCE_METHOD     = 'Cz';                   % 'Cz', 'earlobes', 'average', etc.
N_EPOCHS_PER_SPEED   = 100;                    % epochs kept per speed per participant
ARTIFACT_THRESHOLD   = 20;                     % µV, peak-to-peak rejection threshold
MICROPHONE_DELAY_MS  = 55;                     % acoustic delay correction
```

**Output files:**

| File | Contents |
|---|---|
| `GROUP_SEP_Data_{ref}.mat` | All clean epochs before balancing |
| `BALANCED_SEP_Data_100eps_{ref}.mat` | Balanced epoch set (input for statistics.m) |
| `PEAK_ANALYSIS_Results_{ref}.mat` | P1/N1/P2 amplitudes and latencies from grand average |
| `WPSS_Results_{ref}.mat` | Wavelet phase synchronisation stability per speed |
| `*.png` | Grand average waveforms, heatmaps, WPSS spectrograms |

**Run:**
```matlab
run('processing_data.m')
```

---

### Step 3 – Statistical Analysis

**File:** `statistics.m`

Loads the balanced dataset and runs a complete statistical pipeline on P1, N1, P2 amplitudes and latencies at electrode C4.

**Required MATLAB Toolboxes:**
- Statistics and Machine Learning Toolbox (`fitrm`, `ranova`, `mauchly`, `epsilon`, `swtest` equivalent)

> **Note:** A custom `swtest()` function implementing the Shapiro-Wilk test (Royston 1992) is included at the bottom of the script and requires no additional toolbox.

**Input:** `BALANCED_SEP_Data_100eps_{ref}.mat` (produced by `processing_data.m`)

**Configuration:**
```matlab
REFERENCE_METHOD = 'Cz';   % must match the mat file used
ELECTRODE_FOCUS  = 'C4';
ALPHA            = 0.05;
```

**Analysis pipeline:**
1. Component extraction — P1 (40–90 ms), N1 (100–160 ms), P2 (180–235 ms)
2. One-way repeated measures ANOVA with Greenhouse-Geisser correction
3. Pairwise post-hoc tests with Bonferroni correction and Cohen's *d*
4. Polynomial trend analysis (linear, quadratic, cubic)
5. Normality testing (Shapiro-Wilk)
6. Signal-to-noise ratio (P1 peak / pre-stimulus RMS)
7. Split-half reliability (Spearman-Brown corrected)

**Output files:**

| File | Contents |
|---|---|
| `MULTI_LEAN_Analysis_{ref}.mat` | All statistical results as struct |
| `MULTI_LEAN_Results_{ref}.csv` | Per-participant, per-speed component values |
| `MULTI_LEAN_ThesisReport_{ref}.txt` | Formatted results text ready for thesis |
| `MULTI_GrandAverage_{ref}.png` | Waveforms with component annotations |
| `MULTI_Amplitudes_{ref}.png` | Amplitude × speed with trend lines |
| `MULTI_Latencies_{ref}.png` | Latency × speed |
| `MULTI_IndividualProfiles_{ref}.png` | Per-participant profiles + normalised comparison |
| `MULTI_ANOVA_Summary_{ref}.png` | Effect sizes and power |

**Run:**
```matlab
run('statistics.m')
```

---

## Study 2 – Subjective Paired Comparison

### Overview

A separate, matched participant group completed two paired comparison sessions (intensity and valence). On each trial, two stimuli were presented sequentially and the participant indicated which felt stronger (intensity) or more pleasant (valence). Responses were analysed using the Bradley-Terry model.

### Step 1 – Data Collection

**Files:** `intensity_script.cs` / `valence_script.cs`

Run one file per session. Both scripts share the same structure and differ only in the question asked and the output filename prefix.

**Requirements:** same as `objective_program.cs` (Ultraleap SDK, .NET 6+)

**Design:**
- 10 repetitions × 2 presentation orders × 10 pairs = **200 trials per session**
- Break every 40 trials
- Key mapping: `1` = Stimulus 1, `2` = Stimulus 2, `0` = No clear preference, `Esc` = Abort

**Stimulus parameters:**

| Parameter | Value |
|---|---|
| Duration per stimulus | 200 ms |
| Inter-stimulus interval | 1000 ms |
| Inter-trial interval | 800 ms |
| Circle radius | 32.5 mm |
| Hand distance | 170 mm above device |

**Output:**
- `P{ID}_intensity_{timestamp}.csv`
- `P{ID}_valence_{timestamp}.csv`

CSV columns: `trial, block, participant, speed_first, speed_second, choice, chosen_speed, rt_ms, timestamp`

**Usage:**
```
dotnet run   # intensity_script.cs or valence_script.cs
# Enter participant ID when prompted
```

---

### Step 2 – Paired Comparison Analysis

**File:** `paired_comparison_analysis.m`

Loads all intensity and valence CSV files, fits a Bradley-Terry model, computes win rates, performs per-participant analysis, and checks for circular triads and position bias.

**Required MATLAB Toolboxes:** Statistics and Machine Learning Toolbox (`normcdf`, `corr`)

**Input:** CSV files placed in `data/csv_files_subjective/`

Naming conventions expected:
- `P001_intensity_20240101_120000.csv`
- `P001_valence_20240101_130000.csv`

**Configuration:**
```matlab
dataDir = fullfile(pwd, 'data', 'csv_files_subjective');
speeds  = [1.0, 3.0, 5.0, 7.0, 10.0];
```

**Output:** Figures (Bradley-Terry scale values, win rates, per-participant profiles, intensity vs. valence comparison) displayed in MATLAB. No files are saved automatically — export figures manually as needed.

**Run:**
```matlab
run('paired_comparison_analysis.m')
```

---

## Recommended Execution Order

```
Study 1 (EEG)
  1.  objective_program.cs          →  P00X_trials.csv
  2.  processing_data.m             →  BALANCED_SEP_Data_*.mat, WPSS_Results_*.mat
  3.  statistics.m                  →  MULTI_LEAN_*.mat, figures, report

Study 2 (Subjective)
  1.  intensity_script.cs           →  P00X_intensity_*.csv
  2.  valence_script.cs             →  P00X_valence_*.csv
  3.  paired_comparison_analysis.m  →  figures
```

---

## Dependencies Summary

| Component | Requirement |
|---|---|
| MATLAB | R2024b or later recommended |
| Signal Processing Toolbox | `filtfilt`, `fir1` |
| Wavelet Toolbox | `cwtfilterbank`, `wt` |
| Statistics and Machine Learning Toolbox | `fitrm`, `ranova`, `mauchly`, `epsilon`, `normcdf`, `corr` |
| .NET runtime | 6.0 or later |
| Ultraleap Haptics SDK | Tested with Ultraleap Stratos Explore |

---

## Notes

- The 55 ms microphone delay correction in `processing_data.m` shifts trigger timestamps earlier to compensate for acoustic propagation latency. EEG signals are not modified.
- Balanced epoch selection (`N_EPOCHS_PER_SPEED = 100`) randomly draws an equal number of epochs per speed condition per participant. Participants with fewer than 100 clean epochs in any condition are excluded from the analysis.
- The Bradley-Terry model is fit via iterative maximum likelihood (500 iterations). Log-transformed, mean-centred scale values are reported.
- All decimal values in CSV output use `CultureInfo.InvariantCulture` to ensure `.` as the decimal separator regardless of system locale.

---

## Citation

If you use this code, please cite the associated thesis:

> L. Wolf, *Investigating the Impact of Mid-Air Haptic Feedback Parameters on Human Perception via Somatosensory Evoked Potentials (SEPs)*, Master's Thesis, htw saar, 2026.
