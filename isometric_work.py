"""
isometric_work.py

Description: Process all animals in a folder, integrate all positive normalized 
force over displacement into isometric work (force-time integral), and return a 
figure plus the raw results.
"""

__author__ = "Chaoran Yang"
__version__ = "2.6"
__email__ = "cy197@duke.edu"
__date__ = "2026-03-29"

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

# ---------------------------------------------------------------------------
# Internal helpers for interpolated crossing times
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Core detection: horizontal-line intersection sweep
# ---------------------------------------------------------------------------

def detect_limits_by_intersection(signal: np.ndarray, sample_freq_hz: float,
                                    n_steps: int = 2000,
                                    file_path: str = "",
                                    ) -> Tuple[float, float, float, float]:
    """
    Identify integration limits by sweeping a horizontal line y = n and
    counting how many times it intersects the force-time signal.

    How the crossing-count profile looks
    -------------------------------------
    Sweeping n from 0 upward:

      n < baseline noise floor
          signal is entirely above n  →  0 crossings

      n inside baseline noise band
          baseline oscillations generate many crossings

      n_start  (first n with exactly 1 crossing)
          signal is fully below n during the quiet baseline, crosses n
          exactly ONCE going up during the rise, and stays above n through
          the plateau  →  1 crossing

      [n_start … n_end)  — the "1-crossing band"
          every level in this band is crossed exactly once (the rise)

      n_end  (first n > 1 crossing after the 1-crossing band)
          plateau noise causes the signal to oscillate across n  →  > 1 crossings

      n > max(signal)
          0 crossings

    Integration limits
    ------------------
    START
      1.  n_start  → first n with exactly 1 crossing.
      2.  Indicator time  = first upward crossing of n_start (just inside the rise).
      3.  baseline_mean   = mean( signal[ : indicator_sample ] ).
      4.  start_time      = last upward crossing of baseline_mean before the
                            indicator (interpolated)  →  force = baseline here,
                            so the baseline-corrected integrand starts at 0.

    END
      1.  n_end    → first n with > 1 crossings after the 1-crossing band.
      2.  Indicator time  = first upward crossing of n_end (inside the rise).
      3.  plateau_mean    = mean( signal[ indicator_sample : ] ).
      4.  end_time        = first upward crossing of plateau_mean after the onset
                            (interpolated)  →  force first reaches plateau level.
    Fallback when n_end is not found (very clean plateau, no oscillations):
      plateau_mean estimated from the last 15 % of the signal.

    Parameters
    ----------
    signal          : 1-D force array, already truncated to the contraction
                      window (up to end_idx_force).
    sample_freq_hz  : sampling rate (Hz).
    n_steps         : number of horizontal levels to test (default 2000).
    file_path       : included in diagnostic log messages.

    Returns
    -------
    (start_time, baseline_mean, end_time, plateau_mean)
    Times are in seconds; force values match the units of `signal`.
    """
    signal_max = float(np.max(signal))
    if signal_max <= 0:
        raise RuntimeError(f"{file_path}: signal has no positive values.")

    # Avoid testing exactly 0 or exactly max (endpoint ambiguity)
    n_values = np.linspace(signal_max * 0.001, signal_max * 0.999, n_steps)

    # Vectorised crossing count for all n levels simultaneously.
    # above[i, j] = (signal[i] >= n_values[j])  shape (N, M)
    # |diff| along axis 0 counts every above↔below transition per column.
    above           = signal[:, np.newaxis] >= n_values[np.newaxis, :]  # (N, M) bool
    diffs           = np.diff(above.astype(np.int8), axis=0)             # (N-1, M)
    crossing_counts = np.sum(np.abs(diffs), axis=0)                      # (M,) int

    # ── START ──────────────────────────────────────────────────────────── #
    one_mask = crossing_counts == 1
    if not np.any(one_mask):
        raise RuntimeError(
            f"{file_path}: no horizontal level produced exactly 1 crossing. "
            f"Verify that the recording contains a clear baseline-to-plateau rise."
        )
    first_one_local  = int(np.argmax(one_mask))
    n_start          = float(n_values[first_one_local])

    t_ind_start      = _first_upward_crossing_time(signal, n_start, sample_freq_hz)
    ind_start_sample = max(1, min(int(round(t_ind_start * sample_freq_hz)),
                                   len(signal) - 1))
    baseline_mean    = float(np.mean(signal[:ind_start_sample]))
    start_time       = _last_upward_crossing_before(signal, baseline_mean,
                                                     sample_freq_hz, ind_start_sample)

    # ── END ────────────────────────────────────────────────────────────── #
    n_end = None
    for k in range(first_one_local, len(crossing_counts)):
        if crossing_counts[k] > 1:
            n_end = float(n_values[k])
            break

    if n_end is not None:
        t_ind_end      = _first_upward_crossing_time(signal, n_end, sample_freq_hz,
                                                       after_sample=ind_start_sample)
        ind_end_sample = max(ind_start_sample + 1,
                             min(int(round(t_ind_end * sample_freq_hz)),
                                  len(signal) - 1))
        plateau_mean   = float(np.mean(signal[ind_end_sample:]))
    else:
        # Very clean plateau — fall back to tail of signal
        plateau_n    = max(1, int(0.15 * len(signal)))
        plateau_mean = float(np.mean(signal[-plateau_n:]))

    end_time = _first_upward_crossing_time(signal, plateau_mean, sample_freq_hz,
                                             after_sample=ind_start_sample)

    # ── Diagnostics ────────────────────────────────────────────────────── #
    n_end_str = f"{n_end:.4f}" if n_end is not None else "not found (fallback used)"
    print(f"[detect_limits_by_intersection] {file_path}")
    print(f"  n_start       = {n_start:.4f}  (baseline noise ceiling)")
    print(f"  n_end         = {n_end_str}  (plateau noise floor)")
    print(f"  baseline_mean = {baseline_mean:.4f}")
    print(f"  plateau_mean  = {plateau_mean:.4f}")
    print(f"  start_time    = {start_time:.6f} s")
    print(f"  end_time      = {end_time:.6f} s")

    return float(start_time), float(baseline_mean), float(end_time), float(plateau_mean)


# ---------------------------------------------------------------------------
# Rough onset helper — used ONLY to anchor the contraction window
# (so the length-signal argmin search starts after the onset).
# The precise integration limits come entirely from detect_limits_by_intersection.
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# DMC file parser
# ---------------------------------------------------------------------------

def parse_dmc_file(file_path: str) -> dict:
    """
    Parse a DMCv5.x-style data file and return detected integration limits.

    Returns
    -------
    dict with keys:
      'time'           : time axis (s)
      'length_mm'      : muscle length (mm)
      'force_mN'       : force (mN)
      'raw_df'         : full pandas DataFrame
      'start_time'     : interpolated start time (s) — (force−baseline) = 0
      'end_time'       : interpolated end time (s)   — force = plateau_mean
      'baseline_mean'  : pre-contraction baseline force (mN)
      'plateau_mean'   : plateau force (mN)
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

    # Contraction window: search for length minimum only after rough onset
    rough_idx     = _rough_onset_idx(force_mN, sample_freq_hz)
    end_idx_force = int(rough_idx + np.argmin(length_mm[rough_idx:]))

    # Precise limits via intersection sweep (operates on the contraction window)
    start_time, baseline_mean, end_time, plateau_mean = detect_limits_by_intersection(
        force_mN[:end_idx_force], sample_freq_hz, file_path=file_path
    )

    return {
        "time":          time_s,
        "length_mm":     length_mm,
        "force_mN":      force_mN,
        "raw_df":        df,
        "start_time":    start_time,
        "end_time":      end_time,
        "baseline_mean": baseline_mean,
        "plateau_mean":  plateau_mean,
        "end_plateau":   end_idx_force,
    }


# ---------------------------------------------------------------------------
# Integration
# ---------------------------------------------------------------------------

def isometric_work_from_file(file_path: str) -> float:
    """
    Compute the isometric work (force-time integral) for a single file.

    The integration window runs from start_time to end_time.
    At start_time  the baseline-corrected force = 0   (by construction).
    At end_time    the baseline-corrected force = plateau_mean − baseline_mean.

    Interpolated boundary points are prepended/appended so np.trapezoid
    captures the exact fractional limits rather than the nearest samples.

    Units: mN · s
    """
    parsed        = parse_dmc_file(file_path)
    time          = parsed["time"]
    force_mN      = parsed["force_mN"]
    start_time    = parsed["start_time"]
    end_time      = parsed["end_time"]
    baseline_mean = parsed["baseline_mean"]
    plateau_mean  = parsed["plateau_mean"]

    # Integer samples strictly inside (start_time, end_time)
    start_sample = int(np.searchsorted(time, start_time, side='right'))
    end_sample   = int(np.searchsorted(time, end_time,   side='left')) - 1
    end_sample   = max(start_sample, min(end_sample, len(time) - 1))

    time_seg  = time[start_sample : end_sample + 1]
    force_seg = force_mN[start_sample : end_sample + 1] - baseline_mean

    # Prepend exact start (integrand = 0) and append exact end
    time_int  = np.concatenate([[start_time],  time_seg,  [end_time]])
    force_int = np.concatenate([[0.0],         force_seg,
                                 [plateau_mean - baseline_mean]])

    return float(np.trapezoid(force_int, time_int))


# ---------------------------------------------------------------------------
# Per-folder processing
# ---------------------------------------------------------------------------

def run_isometric_work(folder_path: str, csa_mm2: float) -> List[float]:
    """
    Process one animal folder and compute normalised isometric work (mN·s/mm²)
    for each contraction file.
    """
    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"Folder does not exist: {folder_path}")

    files = [f for f in os.listdir(folder_path)
             if f.lower().endswith(".ddf") and not f.startswith(".")]

    animal_results: List[float] = []
    for f in sorted(files, key=natural_sort_key):
        data_file = os.path.join(folder_path, f)
        print(f"Processing {f} ...")
        work            = isometric_work_from_file(data_file)
        normalized_work = work / csa_mm2
        animal_results.append(normalized_work)

    return animal_results


# ---------------------------------------------------------------------------
# Validation figure
# ---------------------------------------------------------------------------

def val_isometric_work(file_path: str, i: int):
    """
    Produce a force validation figure for a single file.

    Green dot  : (start_time, baseline_mean)  — integration starts here; force = baseline
    Red dot    : (end_time,   plateau_mean)   — integration ends here; force = plateau
    Blue shade : integrated region
    """
    parsed        = parse_dmc_file(file_path)
    time          = parsed["time"]
    force_mN      = parsed["force_mN"]
    start_time    = parsed["start_time"]
    end_time      = parsed["end_time"]
    baseline_mean = parsed["baseline_mean"]
    plateau_mean  = parsed["plateau_mean"]
    end_plateau   = parsed["end_plateau"]

    start_sample = int(np.searchsorted(time, start_time, side='right'))
    end_sample   = int(np.searchsorted(time, end_time,   side='left')) - 1
    end_sample   = max(start_sample, min(end_sample, len(time) - 1))

    time_shade  = np.concatenate([[start_time],
                                   time[start_sample : end_sample + 1],
                                   [end_time]])
    force_shade = np.concatenate([[baseline_mean],
                                   force_mN[start_sample : end_sample + 1],
                                   [plateau_mean]])

    short_path = "/".join(PlPath(file_path).parts[-3:])

    fig_l, ax_l = plt.subplots(figsize=(11, 8))
    ax_l.plot(time[:end_plateau + 1], force_mN[:end_plateau + 1])
    ax_l.fill_between(time_shade, force_shade, baseline_mean,
                      alpha=0.3, label="Integrated Region")
    ax_l.plot(start_time, baseline_mean,
              'o', color='green', markersize=8,
              label=f"Start  (baseline = {baseline_mean:.2f} mN)")
    ax_l.plot(end_time, plateau_mean,
              'o', color='brown', markersize=8,
              label=f"End  (plateau = {plateau_mean:.2f} mN)")
    ax_l.set_xlabel("Time (s)")
    ax_l.set_ylabel("Force (mN)")
    ax_l.set_title(f"{short_path} (index {i})")
    ax_l.legend()
    ax_l.grid(True)

    return fig_l


# ---------------------------------------------------------------------------
# Streamlit UI
# ---------------------------------------------------------------------------

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
        x_coord = np.arange(len(result))
        ax.plot(x_coord, np.array(result), label=f"Folder {idx + 1}")
    ax.set_xlabel("Contraction Index")
    ax.set_ylabel("Normalized Isometric Work (mN·s/mm^2)")
    ax.set_title("Isometric Work")
    ax.legend()
    ax.grid(True)
    st.pyplot(fig)

    rows = []
    for i in range(len(csv_x_labels)):
        labels = csv_x_labels[i]
        values = csv_output[i]
        for fname, val in zip(labels, values):
            rows.append([fname, val])
        rows.append(["", ""])

    final_df  = pd.DataFrame(rows, columns=["filename", "isometric work (mN·s/mm^2)"])
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
