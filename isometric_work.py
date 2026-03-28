"""
isometric_work.py

Description: Process all animals in a folder, integrate all positive normalized 
force over displacement into isometric work (force-time integral), and return a 
figure plus the raw results.
"""

__author__ = "Chaoran Yang"
__version__ = "2.5"
__email__ = "cy197@duke.edu"
__date__ = "2026-03-28"

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
import warnings

def natural_sort_key(s: str):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r"([0-9]+)", s)]

def detect_onset(signal: np.ndarray, sample_freq_hz: float,
                 bootstrap_s: float = 0.08, threshold_std: float = 6.0,
                 smooth_window_s: float = 0.02,
                 min_std_fraction: float = 0.005, sustained_s: float = 0.005,
                 file_path: str = "") -> Tuple[int, float]:
    """
    Detect the first sample where `signal` departs from its initial baseline,
    and return the interpolated time of the last baseline crossing before the rise.

    The onset index confirms the rise is real (sustained threshold criterion).
    The crossing time is found by walking backward from the onset index along the
    smoothed signal until it falls to or below baseline_mean, then linearly
    interpolating the exact crossing moment. This guarantees the integration
    region always starts at (force - baseline) = 0.

    Parameters
    ----------
    signal          : 1D array (force or length)
    sample_freq_hz  : sampling rate
    bootstrap_s     : seconds at the very start used to estimate baseline
    threshold_std   : number of baseline SDs required to declare onset
    smooth_window_s : smoothing kernel width (seconds) to suppress noise
    min_std_fraction: noise floor — baseline_std is clamped to at least this
                      fraction of |baseline_mean|. Prevents the threshold from
                      collapsing to near-zero on ultra-quiet baselines.
    sustained_s     : seconds the signal must stay above threshold
                      continuously before an onset is declared.
    file_path       : included in log messages for traceability

    Returns
    -------
    (onset_idx, crossing_time_s)
      onset_idx       : integer index of first confirmed above-threshold sample
      crossing_time_s : interpolated time (s) where the smoothed signal last
                        crossed baseline_mean — use as the integration start
    """
    bootstrap_n = int(bootstrap_s * sample_freq_hz)
    bootstrap = signal[:bootstrap_n]
    baseline_mean = np.mean(bootstrap)
    baseline_std  = np.std(bootstrap)

    std_floor    = min_std_fraction * abs(baseline_mean) if baseline_mean != 0 else 1e-9
    baseline_std = max(baseline_std, std_floor)

    kernel   = max(1, int(smooth_window_s * sample_freq_hz))
    smoothed = np.convolve(signal, np.ones(kernel) / kernel, mode='same')

    deviation = np.abs(smoothed - baseline_mean)
    threshold = threshold_std * baseline_std
    above     = deviation > threshold

    sustained_n = max(1, int(sustained_s * sample_freq_hz))
    onset_idx   = None
    count       = 0
    for i, flag in enumerate(above):
        if flag:
            count += 1
            if count >= sustained_n:
                onset_idx = i - sustained_n + 1   # back to where run started
                break
        else:
            count = 0

    if onset_idx is None:
        onset_idx = 0

    # Walk backward from onset_idx along the smoothed signal to find the last
    # sample that is at or below baseline_mean — this is the end of the flat
    # baseline just before the rise begins.
    j = onset_idx
    while j > 0 and smoothed[j] > baseline_mean:
        j -= 1

    # Linearly interpolate the exact crossing moment between sample j and j+1.
    s0 = smoothed[j]
    s1 = smoothed[min(j + 1, len(smoothed) - 1)]
    if s1 > s0:
        frac = float(np.clip((baseline_mean - s0) / (s1 - s0), 0.0, 1.0))
    else:
        frac = 0.0
    crossing_time = (j + frac) / sample_freq_hz

    print(f"[detect_onset] {file_path}")
    print(f"  baseline_mean    = {baseline_mean:.4f}")
    print(f"  baseline_std     = {baseline_std:.6f}")
    print(f"  threshold        = {threshold_std} × std = {threshold_std * baseline_std:.6f}")
    print(f"  onset_idx        = {onset_idx}")
    print(f"  crossing_time_s  = {crossing_time:.6f}  (interpolated baseline crossing)")

    return int(onset_idx), float(crossing_time)


def detect_rise_end(signal: np.ndarray, sample_freq_hz: float,
                    onset_idx: int,
                    smooth_window_s: float = 0.015,
                    deriv_drop_fraction: float = 0.10,
                    sustained_s: float = 0.005,
                    plateau_fraction: float = 0.97,
                    file_path: str = "") -> int:
    """
    Detect where the force rise phase ends — the transition from rapid increase
    to plateau or steady state — without assuming the plateau is flat.

    Two-stage strategy
    ------------------
    Stage 1 — Derivative drop (handles overshoot and clean plateau cases):
        Compute the smoothed first derivative of force.  Find the steepest
        point of the rise (peak df/dt).  After that peak, find the first
        moment the derivative falls to ≤ deriv_drop_fraction of its maximum
        and stays there for sustained_s seconds.
        • Overshoot case: df/dt crosses zero at the force peak → detected
          almost immediately.
        • Clean/gradual rise: df/dt asymptotes to near-zero → detected when
          force rate of change becomes negligible.

    Stage 2 — Plateau-fraction fallback (handles still-rising or very
        slowly settling signals where the derivative never clearly drops):
        Estimate the plateau from the tail of the signal and find where
        force first reaches plateau_fraction of that estimate.

    Parameters
    ----------
    signal               : 1D force array (already truncated to the contraction
                           window, e.g. up to end_idx_force)
    sample_freq_hz       : sampling rate
    onset_idx            : integer index of rise onset (into `signal`)
    smooth_window_s      : smoothing kernel width for derivative computation
                           (wider = less noise sensitivity; 15 ms default)
    deriv_drop_fraction  : df/dt must drop to this fraction of its peak before
                           the rise is declared ended (default 0.10 = 10 %)
    sustained_s          : derivative must stay below threshold for this long to
                           avoid noise spikes triggering early detection
    plateau_fraction     : fallback — fraction of estimated plateau force used
                           to define end of rise (default 0.97 = 97 %)
    file_path            : included in log messages for traceability

    Returns
    -------
    end_idx (int) : index into `signal` where the rise phase ends
    """
    N = len(signal)

    kernel   = max(1, int(smooth_window_s * sample_freq_hz))
    smoothed = np.convolve(signal, np.ones(kernel) / kernel, mode='same')

    seg = smoothed[onset_idx:]
    if len(seg) < 3:
        return N - 1

    dseg = np.diff(seg)
    if len(dseg) == 0:
        return N - 1

    # Steepest point of the rise (peak df/dt in post-onset region)
    peak_deriv_local = int(np.argmax(dseg))
    peak_deriv_val   = dseg[peak_deriv_local]

    if peak_deriv_val <= 0:
        # Signal never rose — return end of array as fallback
        return N - 1

    # --- Stage 1: sustained derivative drop ---
    threshold   = deriv_drop_fraction * peak_deriv_val
    sustained_n = max(1, int(sustained_s * sample_freq_hz))

    post    = dseg[peak_deriv_local:]
    count   = 0
    end_local = None
    for k, d in enumerate(post):
        if d <= threshold:
            count += 1
            if count >= sustained_n:
                # Step back to where the sustained run began
                end_local = peak_deriv_local + (k - sustained_n + 1)
                break
        else:
            count = 0

    if end_local is not None:
        end_idx = onset_idx + end_local
        print(f"[detect_rise_end] {file_path}")
        print(f"  method           = derivative drop  "
              f"(threshold = {deriv_drop_fraction*100:.0f}% × peak df/dt)")
        print(f"  peak_deriv_val   = {peak_deriv_val:.6f}")
        print(f"  end_idx          = {end_idx}")
        return int(end_idx)

    # --- Stage 2: plateau-fraction fallback ---
    plateau_n   = max(1, int(0.15 * N))
    plateau_est = float(np.median(smoothed[-plateau_n:]))
    baseline_est = (float(np.mean(signal[:onset_idx]))
                    if onset_idx > 0 else float(smoothed[0]))
    threshold_val = baseline_est + plateau_fraction * (plateau_est - baseline_est)

    above = np.where(seg >= threshold_val)[0]
    if len(above) > 0:
        end_idx = onset_idx + int(above[0])
        print(f"[detect_rise_end] {file_path}")
        print(f"  method           = plateau-fraction fallback  "
              f"({plateau_fraction*100:.0f}% of estimated plateau)")
        print(f"  plateau_est      = {plateau_est:.4f}")
        print(f"  threshold_val    = {threshold_val:.4f}")
        print(f"  end_idx          = {end_idx}")
        return int(end_idx)

    # Last-resort: end of array
    print(f"[detect_rise_end] {file_path}")
    print(f"  method           = last-resort fallback (end of array)")
    return N - 1


def parse_dmc_file(file_path: str) -> dict[str, np.ndarray]:
    """
    Parse a DMCv5.x-style data file.

    Extracts:
      - sample frequency
      - calibration info per channel (scale, offset, units)
      - data table (Sample, AI0, AI1, ...)

    Returns
    -------
    dict with keys:
      'time'            : time axis in seconds
      'length_mm'       : muscle length in mm (using calibrated AI channel)
      'force_mN'        : force in mN (using calibrated AI channel)
      'raw_df'          : full pandas DataFrame of the test data
      'start_idx'       : integer onset index detected from force signal
      'crossing_time'   : interpolated time (s) of last baseline crossing
                          — use as the actual integration start (force = baseline there)
      'end_idx'         : end of rise phase detected from force signal
      'end_plateau'     : index of minimum length (end of contraction window)
    """

    with open(file_path, "r", encoding="latin-1") as f:
        lines = f.readlines()

    data_marker_idx = None
    sample_freq_hz = None

    # Parse header-level info, including calibration block
    channel_names: List[str] = []
    units: List[str] = []
    scales: List[float] = []
    offsets: List[float] = []

    i = 0

    while i < len(lines):
        line = lines[i].strip()

        # Sample frequency
        if line.startswith("Sample Frequency"):
            try:
                sample_freq_hz = float(line.split(":")[1].strip())
            except (IndexError, ValueError):
                raise ValueError(f"Could not parse sample frequency in {file_path}")
            i += 1
            continue

        # Locate data marker
        if line.startswith("Test Data in Volts"):
            data_marker_idx = i
            break

        # Calibration block
        if line.startswith("Channel"):
            channel_names = lines[i].strip().split("\t")[1:]
            units = lines[i + 1].strip().split("\t")[1:]
            scales = [float(x) for x in lines[i + 2].strip().split("\t")[1:]]
            offsets = [float(x) for x in lines[i + 3].strip().split("\t")[1:]]
            tads = [float(x) for x in lines[i + 4].strip().split("\t")[1:]]
            i += 5
            continue

        i += 1

    if data_marker_idx is None:
        raise ValueError(f"'Test Data in Volts:' not found in {file_path}")

    if sample_freq_hz is None:
        raise ValueError(f"Sample frequency not found in {file_path}")

    # Build calibration dictionary
    calib = {}
    for name, unit, scale, offset, tad in zip(channel_names, units, scales, offsets, tads):
        calib[name] = {"units": unit, "scale": scale, "offset": offset, "tad": tad}

    # Read test data using pandas:
    # skip all lines up to and including 'Test Data in Volts:'
    # so the next line ('Sample\tAI0\tAI1...') becomes the header row.
    skiprows = data_marker_idx + 1
    df = pd.read_csv(
        file_path,
        delimiter="\t",
        skiprows=skiprows,
        engine="python"
    )

    # Identify channels:
    #   - Length: first AI channel with units 'mm'
    #   - Force : first AI channel with units 'Ref' and positive scale
    length_chan = None
    force_chan = None

    for name in calib:
        meta = calib[name]
        if name.startswith("AI") and meta["units"] == "mm" and length_chan is None:
            length_chan = name
        if name.startswith("AI") and meta["units"] == "Ref" and meta["scale"] > 0 and force_chan is None:
            force_chan = name

    if length_chan is None or force_chan is None:
        raise ValueError(
            f"Could not identify length/force channels in calibration table for {file_path}"
        )

    # Raw volts
    length_volts = df[length_chan].to_numpy(dtype=float)
    force_volts  = df[force_chan].to_numpy(dtype=float)

    # Convert to physical units: units = (volts - offset) * scale
    length_mm = (length_volts - calib[length_chan]["offset"]) * calib[length_chan]["scale"] + calib[length_chan]["tad"]
    force_ref = (force_volts  - calib[force_chan]["offset"])  * calib[force_chan]["scale"]  + calib[force_chan]["tad"]

    # Treat 'Ref' for force as mN (matches your previous scaling)
    force_mN = force_ref

    # Build time axis
    num_samples    = len(df)
    sample_indices = np.arange(num_samples, dtype=float)
    time_s         = sample_indices / sample_freq_hz

    # Detect onset: returns integer index + interpolated baseline-crossing time
    start_idx, crossing_time = detect_onset(force_mN, sample_freq_hz,
                                            file_path=file_path)

    # Contraction window: from onset to the point of minimum length
    end_idx_force = int(start_idx + np.argmin(length_mm[start_idx:]))

    # Detect end of rise phase within the contraction window
    end_idx = detect_rise_end(force_mN[:end_idx_force], sample_freq_hz,
                              onset_idx=start_idx, file_path=file_path)

    return {
        "time":           time_s,
        "length_mm":      length_mm,
        "force_mN":       force_mN,
        "raw_df":         df,
        "start_idx":      start_idx,
        "crossing_time":  crossing_time,
        "end_idx":        end_idx,
        "end_plateau":    end_idx_force,
    }


def isometric_work_from_file(file_path: str) -> float:
    """
    Compute the isometric work (force-time integral) for a single file.

    The isometric work is the area under the baseline-corrected force curve
    during the rising phase.

    The integration window runs from the interpolated baseline crossing
    (where force = baseline, so the integrand = 0) to the detected end of
    the rise phase.  Prepending the crossing point ensures the integration
    starts at zero regardless of noise or timing resolution.

    Units: mN · s  (force-time integral / impulse)
    """
    parsed         = parse_dmc_file(file_path)
    time           = parsed["time"]
    force_mN       = parsed["force_mN"]
    start_idx      = parsed["start_idx"]
    crossing_time  = parsed["crossing_time"]
    end_idx        = parsed["end_idx"]

    # Baseline from pre-contraction region
    force_baseline = float(np.mean(force_mN[:start_idx]))

    # Baseline-corrected force from onset to end of rise
    force_seg = force_mN[start_idx : end_idx + 1] - force_baseline
    time_seg  = time[start_idx : end_idx + 1]

    # Prepend the interpolated baseline-crossing point.
    # At that point the baseline-corrected force is exactly 0 by construction,
    # so the integration cleanly starts at zero.
    time_int  = np.concatenate([[crossing_time], time_seg])
    force_int = np.concatenate([[0.0],           force_seg])

    total_work = float(np.trapezoid(force_int, time_int))
    return total_work


def run_isometric_work(folder_path: str, csa_mm2: float) -> List[float]:
    """
    Process one animal folder and compute normalized isometric work (mN·s/mm²)
    for each contraction file.

    Parameters
    ----------
    folder_path : str
    csa_mm2     : float  (muscle cross-sectional area in mm²)

    Returns
    -------
    animal_results : list of normalized isometric work values (mN·s/mm²)
    """
    animal_results: List[float] = []

    if not os.path.isdir(folder_path):
        raise FileNotFoundError(f"Folder does not exist: {folder_path}")

    files = [
        f for f in os.listdir(folder_path)
        if f.lower().endswith(".ddf") and not f.startswith(".")
    ]

    for f in sorted(files, key=natural_sort_key):
        data_file = os.path.join(folder_path, f)
        print(f"Processing {f} ...")

        work = isometric_work_from_file(data_file)
        normalized_work = work / csa_mm2   # mN·s/mm^2
        animal_results.append(normalized_work)

    return animal_results


def val_isometric_work(file_path: str, i: int):
    """
    Produce a force validation figure for a single file.

    - Green dot : interpolated baseline crossing (true integration start,
                  where baseline-corrected force = 0)
    - Red dot   : end of rise phase (integration end)
    - Blue shading: integrated region
    """
    parsed        = parse_dmc_file(file_path)
    time          = parsed["time"]
    force_mN      = parsed["force_mN"]
    start_idx     = parsed["start_idx"]
    crossing_time = parsed["crossing_time"]
    end_idx       = parsed["end_idx"]
    end_plateau   = parsed["end_plateau"]

    force_baseline = float(np.mean(force_mN[:start_idx]))

    short_path = "/".join(PlPath(file_path).parts[-3:])

    # Build the shaded integration region with the prepended baseline-crossing point
    time_plot  = np.concatenate([[crossing_time],
                                  time[start_idx : end_idx + 1]])
    force_plot = np.concatenate([[force_baseline],
                                  force_mN[start_idx : end_idx + 1]])

    fig_l, ax_l = plt.subplots(figsize=(11, 8))
    ax_l.plot(time[: end_plateau + 1], force_mN[: end_plateau + 1])
    ax_l.fill_between(time_plot, force_plot, force_baseline,
                      alpha=0.3, label="Integrated Region")

    # Green dot: interpolated baseline crossing (guaranteed at force = baseline)
    ax_l.plot(crossing_time, force_baseline,
              'o', color='green', label="Start Point (baseline crossing)")

    # Red dot: end of rise phase
    ax_l.plot(time[end_idx], force_mN[end_idx],
              'o', color='brown', label="End Point (rise end)")

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

    st.session_state.csv_output   = csv_output
    st.session_state.csv_x_labels = csv_x_labels
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