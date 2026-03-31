"""
isometric_work.py

Description: Process all animals in a folder, integrate all positive normalized 
force over displacement into isometric work (force-time integral), and return a 
figure plus the raw results.
"""

__author__ = "Chaoran Yang"
__version__ = "2.7"
__email__ = "cy197@duke.edu"
__date__ = "2026-03-30"

import os
import re
from typing import List, Tuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import zipfile
from datetime import datetime
import random
from pathlib import Path as PlPath

def natural_sort_key(s: str):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r"([0-9]+)", s)]

def _first_upward_crossing_time(signal: np.ndarray, level: float,
                                  sample_freq_hz: float,
                                  after_sample: int = 0) -> float:
    """
    Interpolated time (s) of the first upward crossing of `level` in `signal`,
    starting the search at `after_sample`.
    Falls back to the time of the first sample at-or-above `level`, then to the
    last sample if `level` is never reached.
    """
    n = len(signal)
    for i in range(max(0, after_sample), n - 1):
        if signal[i] < level and signal[i + 1] >= level:
            frac = (level - signal[i]) / (signal[i + 1] - signal[i])
            return (i + frac) / sample_freq_hz
    for i in range(max(0, after_sample), n):
        if signal[i] >= level:
            return i / sample_freq_hz
    return (n - 1) / sample_freq_hz

def _last_upward_crossing_before(signal: np.ndarray, level: float,
                                   sample_freq_hz: float,
                                   before_sample: int) -> float:
    """
    Interpolated time (s) of the LAST upward crossing of `level` strictly
    before `before_sample`. Scanning backward guarantees we pick the crossing
    that is closest in time to the actual rise, so the integration region
    begins right at the edge of the baseline with (force − baseline) = 0.
    Falls back to time 0 if no upward crossing is found.
    """
    limit = min(before_sample, len(signal) - 1)
    for i in range(limit, 0, -1):
        if signal[i - 1] < level and signal[i] >= level:
            frac = (level - signal[i - 1]) / (signal[i] - signal[i - 1])
            return (i - 1 + frac) / sample_freq_hz
    return 0.0

def detect_limits_by_intersection(signal: np.ndarray, sample_freq_hz: float,
                                    n_steps: int = 2000,
                                    file_path: str = "",
                                    ) -> Tuple[float, float]:
    """
    Identify the integration START by sweeping a horizontal line y = n and
    counting intersections with the force-time signal.

    END detection is NOT performed here — it is derived from the length signal
    (see _length_shortening_onset in parse_dmc_file).

    How the crossing-count profile looks
    -------------------------------------
    Sweeping n from 0 upward:

      n inside baseline noise band  ->  many crossings
      n_start (first n with 1 crossing)  ->  signal below n during quiet baseline,
                                              crosses n exactly once on the way up

    Integration start
    -----------------
      1.  n_start  → first n with exactly 1 crossing (on the pre-smoothed signal).
      2.  Indicator time  = first upward crossing of n_start.
      3.  baseline_mean   = mean( signal[ : indicator_sample ] ).
      4.  start_time      = last upward crossing of baseline_mean before the
                            indicator (interpolated) → force = baseline here,
                            so the baseline-corrected integrand starts at 0.

    Parameters
    ----------
    signal          : 1-D force array, truncated to the contraction window.
    sample_freq_hz  : sampling rate (Hz).
    n_steps         : number of horizontal levels to test (default 2000).
    file_path       : included in diagnostic log messages.

    Returns
    -------
    (start_time, baseline_mean)
    """
    signal_max = float(np.max(signal))
    if signal_max <= 0:
        raise RuntimeError(f"{file_path}: signal has no positive values.")

    # Pre-smooth before counting to suppress noise-spike false crossings.
    # Edge-padding avoids boundary roll-off from mode='same'.
    # Raw signal is used for all interpolation steps.
    smooth_k  = max(1, int(0.020 * sample_freq_hz))   # 20 ms
    padded    = np.pad(signal, smooth_k // 2, mode='edge')
    signal_sm = np.convolve(padded, np.ones(smooth_k) / smooth_k, mode='valid')[:len(signal)]

    # Avoid testing exactly 0 or exactly max (endpoint ambiguity)
    n_values = np.linspace(signal_max * 0.001, signal_max * 0.999, n_steps)

    # Vectorised crossing count for all n levels simultaneously.
    above           = signal_sm[:, np.newaxis] >= n_values[np.newaxis, :]  # (N, M)
    diffs           = np.diff(above.astype(np.int8), axis=0)               # (N-1, M)
    crossing_counts = np.sum(np.abs(diffs), axis=0)                        # (M,)

    # ── START ──────────────────────────────────────────────────────────── #
    one_mask = crossing_counts == 1
    if not np.any(one_mask):
        raise RuntimeError(
            f"{file_path}: no horizontal level produced exactly 1 crossing. "
            f"Verify that the recording contains a clear baseline-to-plateau rise."
        )
    first_one_local  = int(np.argmax(one_mask))
    n_start          = float(n_values[first_one_local])

    # Smoothed signal for indicator (avoids noise-spike false crossings);
    # raw signal for the final interpolated start_time.
    t_ind_start      = _first_upward_crossing_time(signal_sm, n_start, sample_freq_hz)
    ind_start_sample = max(1, min(int(round(t_ind_start * sample_freq_hz)),
                                   len(signal) - 1))
    baseline_mean    = float(np.mean(signal[:ind_start_sample]))
    start_time       = _last_upward_crossing_before(signal, baseline_mean,
                                                     sample_freq_hz, ind_start_sample)

    # ── Diagnostics ────────────────────────────────────────────────────── #
    print(f"[detect_limits_by_intersection] {file_path}")
    print(f"  n_start        = {n_start:.4f}  (baseline noise ceiling)")
    print(f"  baseline_mean  = {baseline_mean:.4f}")
    print(f"  start_time     = {start_time:.6f} s")

    return float(start_time), float(baseline_mean)

def _length_shortening_onset(length_mm: np.ndarray, sample_freq_hz: float,
                              onset_sample: int, file_path: str = "") -> float:
    """
    Find the interpolated time at which the length signal first dips below
    its pre-contraction baseline — the moment shortening begins and the
    isometric phase ends.

    Algorithm
    ---------
    1. Baseline mean and std computed from length_mm[ : onset_sample ].
    2. Threshold = baseline_mean − 3 × baseline_std.
       This tolerates small noise fluctuations while reliably catching the
       onset of genuine shortening (which is typically several tenths of a mm).
    3. Lightly smooth the length signal (20 ms) to suppress sub-noise spikes.
    4. Scan forward from onset_sample for the first sample that falls below
       the threshold, then linearly interpolate the exact crossing time.

    Parameters
    ----------
    length_mm      : full-recording length array (mm)
    sample_freq_hz : sampling rate (Hz)
    onset_sample   : integer index of force onset — search starts here
    file_path      : for diagnostic messages

    Returns
    -------
    Interpolated time (s) of shortening onset.
    Falls back to the time of the global length minimum if no crossing is found.
    """
    N = len(length_mm)
    if onset_sample < 1:
        onset_sample = 1

    baseline_mean = float(np.mean(length_mm[:onset_sample]))
    baseline_std  = float(np.std(length_mm[:onset_sample]))
    # Floor the std so a near-zero noise level doesn't produce a hair-trigger
    std_floor     = max(baseline_std, 0.001)
    threshold     = baseline_mean - 6.0 * std_floor

    smooth_k  = max(1, int(0.020 * sample_freq_hz))   # 20 ms
    padded    = np.pad(length_mm, smooth_k // 2, mode='edge')
    length_sm = np.convolve(padded, np.ones(smooth_k) / smooth_k, mode='valid')[:N]

    for i in range(onset_sample, N - 1):
        if length_sm[i] >= threshold and length_sm[i + 1] < threshold:
            frac = (threshold - length_sm[i]) / (length_sm[i + 1] - length_sm[i])
            t    = (i + frac) / sample_freq_hz
            print(f"[_length_shortening_onset] {file_path}")
            print(f"  length_baseline = {baseline_mean:.4f} mm  "
                  f"threshold = {threshold:.4f} mm  (−6σ = −{6*std_floor:.4f})")
            print(f"  shortening onset = {t:.6f} s  (sample {i})")
            return float(t)

    # Fallback: no crossing found — use global minimum after onset
    fallback_idx = onset_sample + int(np.argmin(length_mm[onset_sample:]))
    t_fallback   = fallback_idx / sample_freq_hz
    print(f"[_length_shortening_onset] {file_path}: WARNING — no threshold crossing "
          f"found; falling back to length minimum at {t_fallback:.4f} s")
    return float(t_fallback)

def _plateau_slope_flag(force_seg: np.ndarray, sample_freq_hz: float) -> bool:
    """
    Return True if force_seg shows a significant upward drift, indicating
    the machine plateau was still rising (Image 2 style).

    The drift (last-quarter mean minus first-quarter mean) is compared to the
    detrended residual std (noise floor).  Flagged when drift > 0 AND
    drift / noise_std > 2.0.

    This is purely an annotation for the CSV — it does not affect end_time.
    """
    if len(force_seg) < 8:
        return False
    q          = max(1, len(force_seg) // 4)
    drift      = float(np.mean(force_seg[-q:])) - float(np.mean(force_seg[:q]))
    t_seg      = np.arange(len(force_seg), dtype=float) / sample_freq_hz
    slope, intercept = np.polyfit(t_seg, force_seg, 1)
    residuals  = force_seg - (slope * t_seg + intercept)
    noise_std  = float(np.std(residuals))
    drift_snr  = drift / noise_std if noise_std > 0 else 0.0
    return bool(drift > 0 and drift_snr > 2.0)

def _rough_onset_idx(signal: np.ndarray, sample_freq_hz: float,
                     bootstrap_s: float = 0.08, threshold_std: float = 6.0,
                     smooth_window_s: float = 0.02,
                     min_std_fraction: float = 0.005,
                     sustained_s: float = 0.005) -> int:
    bootstrap_n   = int(bootstrap_s * sample_freq_hz)
    baseline_mean = np.mean(signal[:bootstrap_n])
    baseline_std  = np.std(signal[:bootstrap_n])
    std_floor     = min_std_fraction * abs(baseline_mean) if baseline_mean != 0 else 1e-9
    baseline_std  = max(baseline_std, std_floor)

    kernel   = max(1, int(smooth_window_s * sample_freq_hz))
    smoothed = np.convolve(signal, np.ones(kernel) / kernel, mode='same')

    threshold   = threshold_std * baseline_std
    above       = np.abs(smoothed - baseline_mean) > threshold
    sustained_n = max(1, int(sustained_s * sample_freq_hz))
    count = 0
    for i, flag in enumerate(above):
        if flag:
            count += 1
            if count >= sustained_n:
                return i - sustained_n + 1
        else:
            count = 0
    return 0



def parse_dmc_file(file_path: str) -> dict:
    """
    Parse a DMCv5.x-style data file and return detected integration limits.

    Integration window
    ------------------
    START : last moment before the force rise where force = baseline
            (from detect_limits_by_intersection on the force signal).
    END   : first moment the length signal dips below its pre-contraction
            baseline — i.e. when shortening begins and isometric work ends
            (from _length_shortening_onset on the length signal).

    Returns
    -------
    dict with keys:
      'time'           : time axis (s)
      'length_mm'      : muscle length (mm)
      'force_mN'       : force (mN)
      'raw_df'         : full pandas DataFrame
      'start_time'     : interpolated start time (s) — (force − baseline) = 0
      'end_time'       : interpolated end time (s) — shortening onset
      'baseline_mean'  : pre-contraction force baseline (mN)
      'slope_flag'     : True if the isometric plateau showed significant upward
                         drift (annotation only; does not affect end_time)
      'end_plateau'    : index of minimum length (contraction window limit)
    """
    with open(file_path, "r", encoding="latin-1") as f:
        lines = f.readlines()

    data_marker_idx = None
    sample_freq_hz  = None
    channel_names: List[str] = []
    units:   List[str]   = []
    scales:  List[float] = []
    offsets: List[float] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("Sample Frequency"):
            try:
                sample_freq_hz = float(line.split(":")[1].strip())
            except (IndexError, ValueError):
                raise ValueError(f"Could not parse sample frequency in {file_path}")
            i += 1
            continue
        if line.startswith("Test Data in Volts"):
            data_marker_idx = i
            break
        if line.startswith("Channel"):
            channel_names = lines[i].strip().split("\t")[1:]
            units   = lines[i + 1].strip().split("\t")[1:]
            scales  = [float(x) for x in lines[i + 2].strip().split("\t")[1:]]
            offsets = [float(x) for x in lines[i + 3].strip().split("\t")[1:]]
            tads    = [float(x) for x in lines[i + 4].strip().split("\t")[1:]]
            i += 5
            continue
        i += 1

    if data_marker_idx is None:
        raise ValueError(f"'Test Data in Volts:' not found in {file_path}")
    if sample_freq_hz is None:
        raise ValueError(f"Sample frequency not found in {file_path}")

    calib = {}
    for name, unit, scale, offset, tad in zip(channel_names, units, scales, offsets, tads):
        calib[name] = {"units": unit, "scale": scale, "offset": offset, "tad": tad}

    skiprows = data_marker_idx + 1
    df = pd.read_csv(file_path, delimiter="\t", skiprows=skiprows, engine="python")

    length_chan = None
    force_chan  = None
    for name in calib:
        meta = calib[name]
        if name.startswith("AI") and meta["units"] == "mm"  and length_chan is None:
            length_chan = name
        if name.startswith("AI") and meta["units"] == "Ref" and meta["scale"] > 0 \
                and force_chan is None:
            force_chan = name

    if length_chan is None or force_chan is None:
        raise ValueError(
            f"Could not identify length/force channels in calibration table for {file_path}"
        )

    length_volts = df[length_chan].to_numpy(dtype=float)
    force_volts  = df[force_chan].to_numpy(dtype=float)

    length_mm = ((length_volts - calib[length_chan]["offset"])
                 * calib[length_chan]["scale"] + calib[length_chan]["tad"])
    force_mN  = ((force_volts  - calib[force_chan]["offset"])
                 * calib[force_chan]["scale"]  + calib[force_chan]["tad"])

    time_s = np.arange(len(df), dtype=float) / sample_freq_hz

    # Contraction window: search for length minimum only after rough force onset
    rough_idx     = _rough_onset_idx(force_mN, sample_freq_hz)
    end_idx_force = int(rough_idx + np.argmin(length_mm[rough_idx:]))

    # ── START: from force signal ────────────────────────────────────────── #
    start_time, baseline_mean = detect_limits_by_intersection(
        force_mN[:end_idx_force], sample_freq_hz, file_path=file_path
    )

    # ── END: from length signal ─────────────────────────────────────────── #
    # End of isometric work = start of shortening = first moment length
    # drops below its pre-contraction baseline.
    onset_sample = max(1, int(round(start_time * sample_freq_hz)))
    end_time     = _length_shortening_onset(
        length_mm, sample_freq_hz, onset_sample, file_path=file_path
    )

    # ── SLOPE FLAG: annotate the force plateau shape ─────────────────────── #
    # Extract the force segment between start and end, compute plateau flag.
    # This is annotation only — end_time is already set from length above.
    end_sample  = min(int(round(end_time * sample_freq_hz)), len(force_mN) - 1)
    force_plateau_seg = force_mN[end_sample:]   # force after shortening onset
    slope_flag  = _plateau_slope_flag(force_plateau_seg, sample_freq_hz)

    return {
        "time":          time_s,
        "length_mm":     length_mm,
        "force_mN":      force_mN,
        "raw_df":        df,
        "start_time":    start_time,
        "end_time":      end_time,
        "baseline_mean": baseline_mean,
        "slope_flag":    slope_flag,
        "end_plateau":   end_idx_force,
    }

def isometric_work_from_file(file_path: str) -> Tuple[float, bool]:
    """
    Compute the isometric work (force-time integral) for a single file.

    Integration window
    ------------------
    [start_time, end_time] where:
      start_time : last moment before the rise where force = baseline
                   (baseline-corrected integrand = 0 here by construction)
      end_time   : first moment the length signal drops below its baseline
                   (shortening onset; end of isometric phase)

    The integrand is force − baseline_mean.  The exact boundary values are
    interpolated so the trapezoidal integral captures fractional samples.

    Returns
    -------
    (work_mN_s, slope_flag)
      work_mN_s  : force-time integral (mN · s), baseline corrected
      slope_flag : True if the isometric force plateau showed significant
                   upward drift (annotation for CSV; does not affect the area)
    """
    parsed        = parse_dmc_file(file_path)
    time          = parsed["time"]
    force_mN      = parsed["force_mN"]
    start_time    = parsed["start_time"]
    end_time      = parsed["end_time"]
    baseline_mean = parsed["baseline_mean"]
    slope_flag    = parsed["slope_flag"]

    # Integer samples strictly inside (start_time, end_time)
    start_sample = int(np.searchsorted(time, start_time, side='right'))
    end_sample   = int(np.searchsorted(time, end_time,   side='left')) - 1
    end_sample   = max(start_sample, min(end_sample, len(time) - 1))

    time_seg  = time[start_sample : end_sample + 1]
    force_seg = force_mN[start_sample : end_sample + 1] - baseline_mean

    # Interpolate the force value at end_time for the closing boundary point
    if end_sample < len(time) - 1:
        t0, t1   = time[end_sample], time[end_sample + 1]
        f0, f1   = force_mN[end_sample] - baseline_mean, force_mN[end_sample + 1] - baseline_mean
        frac     = (end_time - t0) / (t1 - t0) if t1 > t0 else 0.0
        force_at_end = f0 + frac * (f1 - f0)
    else:
        force_at_end = float(force_mN[end_sample] - baseline_mean)

    # Prepend exact start (integrand = 0) and append interpolated end value
    time_int  = np.concatenate([[start_time],  time_seg,  [end_time]])
    force_int = np.concatenate([[0.0],         force_seg, [force_at_end]])

    return float(np.trapezoid(force_int, time_int)), bool(slope_flag)

def run_isometric_work(folder_path: str, csa_mm2: float) -> List[Tuple[float, bool]]:
    """
    Process one animal folder and compute normalised isometric work (mN·s/mm²)
    for each contraction file.

    Returns a list of (normalized_work, slope_flag) tuples.
    slope_flag is True when the plateau had a significant positive drift.
    """
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"Folder does not exist: {folder_path}")

    files = [f for f in os.listdir(folder_path)
             if f.lower().endswith(".ddf") and not f.startswith(".")]

    animal_results: List[Tuple[float, bool]] = []
    for f in sorted(files, key=natural_sort_key):
        data_file = os.path.join(folder_path, f)
        print(f"Processing {f} ...")
        work, slope_flag = isometric_work_from_file(data_file)
        animal_results.append((work / csa_mm2, slope_flag))

    return animal_results

def val_isometric_work(file_path: str, i: int):
    """
    Produce a validation figure for a single file showing both force and length.

    Force panel:
      Green dot  : (start_time, baseline_mean) — integration start; force = baseline
      Brown dot  : (end_time, force at end_time) — shortening onset; end of isometric phase
      Blue shade : integrated region

    Length panel:
      Orange dot : (end_time, length_baseline) — the crossing point that defines end_time
    """
    parsed        = parse_dmc_file(file_path)
    time          = parsed["time"]
    force_mN      = parsed["force_mN"]
    length_mm     = parsed["length_mm"]
    start_time    = parsed["start_time"]
    end_time      = parsed["end_time"]
    baseline_mean = parsed["baseline_mean"]
    end_plateau   = parsed["end_plateau"]

    # Interpolate force at end_time for the closing boundary
    start_sample = int(np.searchsorted(time, start_time, side='right'))
    end_sample   = int(np.searchsorted(time, end_time,   side='left')) - 1
    end_sample   = max(start_sample, min(end_sample, len(time) - 1))

    if end_sample < len(time) - 1:
        t0, t1 = time[end_sample], time[end_sample + 1]
        f0, f1 = force_mN[end_sample], force_mN[end_sample + 1]
        frac   = (end_time - t0) / (t1 - t0) if t1 > t0 else 0.0
        force_at_end = f0 + frac * (f1 - f0)
    else:
        force_at_end = float(force_mN[end_sample])

    # Length baseline (pre-onset mean)
    onset_sample    = max(1, int(round(start_time * (len(time) / time[-1]))))
    length_baseline = float(np.mean(length_mm[:onset_sample]))

    time_shade  = np.concatenate([[start_time],
                                   time[start_sample : end_sample + 1],
                                   [end_time]])
    force_shade = np.concatenate([[baseline_mean],
                                   force_mN[start_sample : end_sample + 1],
                                   [force_at_end]])

    short_path = "/".join(PlPath(file_path).parts[-3:])
    window = slice(0, end_plateau + 1)

    fig_l, (ax_f, ax_l) = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
    fig_l.suptitle(f"{short_path} (index {i})")

    # ── Force panel ─────────────────────────────────────────────────────── #
    ax_f.plot(time[window], force_mN[window])
    ax_f.fill_between(time_shade, force_shade, baseline_mean,
                      alpha=0.3, label="Integrated Region")
    ax_f.plot(start_time, baseline_mean,
              'o', color='green', markersize=8,
              label=f"Start  (force = {baseline_mean:.2f} mN)")
    ax_f.plot(end_time, force_at_end,
              'o', color='brown', markersize=8,
              label=f"End  (shortening onset, force = {force_at_end:.2f} mN)")
    ax_f.set_ylabel("Force (mN)")
    ax_f.legend(fontsize=8)
    ax_f.grid(True)

    # ── Length panel ─────────────────────────────────────────────────────── #
    ax_l.plot(time[window], length_mm[window])
    ax_l.axhline(length_baseline, color='grey', linestyle='--', linewidth=0.8,
                 label=f"Length baseline ({length_baseline:.3f} mm)")
    ax_l.plot(end_time, length_baseline,
              'o', color='orange', markersize=8,
              label=f"Shortening onset ({end_time:.4f} s)")
    ax_l.set_xlabel("Time (s)")
    ax_l.set_ylabel("Length (mm)")
    ax_l.legend(fontsize=8)
    ax_l.grid(True)

    fig_l.tight_layout()
    return fig_l



st.title("Isometric Work Analysis")

uploaded_zip = st.file_uploader("Upload a .zip file", type="zip")

folder_structure = """
The folder you upload should be in the format:
    name_of_folder.zip
      |
       -> animal_1
           |
            -> ___.ddf
            -> ...
       -> animal_2
       -> animal_3
       -> ...

- don't change the name of your file when compressing it
"""
st.code(folder_structure, language='text')

if uploaded_zip:
    upload_folder = "uploaded_files"
    if not os.path.exists(upload_folder):
        os.makedirs(upload_folder)

    zip_file_path = os.path.join(upload_folder, uploaded_zip.name)
    with open(zip_file_path, "wb") as f:
        f.write(uploaded_zip.getbuffer())

    st.success(f"File {uploaded_zip.name} uploaded successfully!")

    unzip_folder = os.path.join(upload_folder, uploaded_zip.name.split('.')[0])
    with zipfile.ZipFile(uploaded_zip, 'r') as zip_ref:
        for member in zip_ref.namelist():
            if '__MACOSX' not in member:
                zip_ref.extract(member, unzip_folder)

    unzip_folder = os.path.join(unzip_folder, uploaded_zip.name.split('.')[0])

    st.write(f"Contents of {unzip_folder}:")
    sorted_folder_names = sorted(os.listdir(unzip_folder), key=natural_sort_key)
    st.write(sorted_folder_names)

    csa_mm2_list: List[float] = []

    st.subheader("Cross-Sectional Area (CSA)")
    mass_mode = st.radio("CSA input mode", ["Manual Input", "Upload CSV"], horizontal=True)

    if mass_mode == "Manual Input":
        for i, foldername in enumerate(sorted_folder_names):
            value = st.number_input(
                f"{foldername} CSA (mm²):", min_value=0.0, value=1.0,
                step=0.00001, key=f"csa_{i}"
            )
            csa_mm2_list.append(value)

        csa_df = pd.DataFrame({"folder": sorted_folder_names, "csa_mm2": csa_mm2_list})
        st.download_button(
            label="Save CSA data as CSV",
            data=csa_df.to_csv(index=False).encode("utf-8"),
            file_name="csa_data.csv",
            mime="text/csv",
        )

    else:
        csa_csv = st.file_uploader("Upload CSA CSV", type="csv")
        if csa_csv is not None:
            uploaded_csa_df = pd.read_csv(csa_csv)
            if {"folder", "csa_mm2"}.issubset(uploaded_csa_df.columns):
                uploaded_csa_df = uploaded_csa_df.set_index("folder")
                missing = [f for f in sorted_folder_names if f not in uploaded_csa_df.index]
                if missing:
                    st.error(f"These folders have no CSA entry in the uploaded CSV: {missing}")
                else:
                    for foldername in sorted_folder_names:
                        csa_mm2_list.append(float(uploaded_csa_df.loc[foldername, "csa_mm2"]))
                    st.success("CSA data loaded from CSV.")
                    st.dataframe(uploaded_csa_df.loc[sorted_folder_names])
            else:
                st.error('CSV must have columns "folder" and "csa_mm2".')

for key in ["analysis_done", "csv_output", "animal_folders",
            "csv_name", "csv_x_labels"]:
    if key not in st.session_state:
        st.session_state[key] = False if key == "analysis_done" else None

def list_subfolders(parent_dir: str):
    return sorted(
        [os.path.join(parent_dir, d) for d in os.listdir(parent_dir)
         if os.path.isdir(os.path.join(parent_dir, d))],
        key=natural_sort_key
    )

def list_files(folder: str):
    return sorted(
        [os.path.join(folder, f) for f in os.listdir(folder)
         if os.path.isfile(os.path.join(folder, f))],
        key=natural_sort_key
    )

if st.button("Run Analysis"):
    st.write("Calculating...")

    csv_output   = []
    csv_x_labels = []

    for i, foldername in enumerate(sorted(os.listdir(unzip_folder), key=natural_sort_key)):
        run_path = os.path.join(unzip_folder, foldername)
        results  = run_isometric_work(run_path, csa_mm2=csa_mm2_list[i])
        csv_output.append(results)

        csv_x_labels.append([])
        for filename in sorted(os.listdir(run_path), key=natural_sort_key):
            if filename.lower().endswith(".ddf"):
                csv_x_labels[i].append(filename)

    st.session_state.csv_output     = csv_output
    st.session_state.csv_x_labels   = csv_x_labels
    st.session_state.animal_folders = list_subfolders(unzip_folder)

    now = datetime.now()
    st.session_state.csv_name = f"isometric_work_{now.strftime('%Y_%m_%d_%H_%M_%S')}.csv"
    st.session_state.analysis_done = True

if st.session_state.analysis_done:
    csv_output     = st.session_state.csv_output
    csv_x_labels   = st.session_state.csv_x_labels
    animal_folders = st.session_state.animal_folders or []

    st.write("Graphing...")
    fig, ax = plt.subplots(figsize=(11, 8))
    for idx, result in enumerate(csv_output):
        values  = [v for v, _ in result]
        x_coord = np.arange(len(values))
        ax.plot(x_coord, np.array(values), label=f"Folder {idx + 1}")
    ax.set_xlabel("Contraction Index")
    ax.set_ylabel("Normalized Isometric Work (mN·s/mm^2)")
    ax.set_title("Isometric Work")
    ax.legend()
    ax.grid(True)
    st.pyplot(fig)

    rows = []
    for i in range(len(csv_x_labels)):
        labels = csv_x_labels[i]
        for fname, (val, slope_flag) in zip(labels, csv_output[i]):
            rows.append([fname, val, "*" if slope_flag else ""])
        rows.append(["", "", ""])

    final_df  = pd.DataFrame(rows, columns=[
        "filename",
        "isometric work (mN·s/mm^2)",
        "plateau drift",
    ])
    csv_bytes = final_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download CSV",
        data=csv_bytes,
        file_name=st.session_state.csv_name,
        mime="text/csv",
    )

    st.subheader("Validate Results")

    if len(animal_folders) == 0:
        st.warning("No animal subfolders found.")
    else:
        animal_choice = st.selectbox(
            "Select animal folder",
            options=animal_folders,
            format_func=lambda p: os.path.basename(p),
        )

        open_animal_folder = [f for f in list_files(animal_choice)
                              if f.lower().endswith(".ddf")]

        if len(open_animal_folder) == 0:
            st.warning("No .ddf files found in the selected animal folder.")
        else:
            mode = st.radio(
                "Validation mode",
                ["Random Sample", "Specific Inspection"],
                horizontal=True,
            )

            if mode == "Random Sample":
                pct       = st.slider("Sample fraction", 0.0, 1.0, 0.05, 0.01)
                n_samples = max(1, int(pct * len(open_animal_folder)))

                if st.button("Run Random Sample"):
                    random_indices = sorted(random.sample(
                        range(len(open_animal_folder)),
                        k=min(n_samples, len(open_animal_folder)),
                    ))
                    for idx in random_indices:
                        val_path = open_animal_folder[idx]
                        fig_l = val_isometric_work(val_path, idx)
                        st.pyplot(fig_l)

            else:
                start = st.number_input(
                    "Start index (inclusive)", min_value=0,
                    max_value=len(open_animal_folder) - 1, value=0, step=1,
                )
                end = st.number_input(
                    "End index (inclusive)", min_value=0,
                    max_value=len(open_animal_folder) - 1,
                    value=len(open_animal_folder) - 1, step=1,
                )

                if start > end:
                    st.error("Start index must be ≤ end index.")
                else:
                    if st.button("Run Specific Inspection"):
                        for idx in range(int(start), int(end) + 1):
                            val_path = open_animal_folder[idx]
                            fig_l = val_isometric_work(val_path, idx)
                            st.pyplot(fig_l)