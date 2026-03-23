"""
isometric_work.py

Description: Process all animals in a folder, integrate all positive normalized 
force over displacement into isometric work (force-time integral), and return a 
figure plus the raw results.
"""

__author__ = "Chaoran Yang"
__version__ = "2.3"
__email__ = "cy197@duke.edu"
__date__ = "2026-03-18"

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
from scipy.optimize import curve_fit, OptimizeWarning
import warnings

def natural_sort_key(s: str):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r"([0-9]+)", s)]

def detect_onset(signal: np.ndarray, sample_freq_hz: float,
                 bootstrap_s: float = 0.08, threshold_std: float = 4.0,
                 smooth_window_s: float = 0.02) -> int:
    """
    Detect the first sample where `signal` departs from its initial baseline.

    Parameters
    ----------
    signal          : 1D array (force or length)
    sample_freq_hz  : sampling rate
    bootstrap_s     : seconds at the very start used to estimate baseline
    threshold_std   : number of baseline SDs required to declare onset
    smooth_window_s : smoothing kernel width (seconds) to suppress noise

    Returns
    -------
    onset index (int)
    """
    bootstrap_n = int(bootstrap_s * sample_freq_hz)
    bootstrap = signal[:bootstrap_n]
    print("(cy 03/18) bootstrap =",bootstrap)
    baseline_mean = np.mean(bootstrap)
    baseline_std  = np.std(bootstrap)

    kernel = max(1, int(smooth_window_s * sample_freq_hz))
    smoothed = np.convolve(signal, np.ones(kernel) / kernel, mode='same')

    deviation = np.abs(smoothed - baseline_mean)
    crossings = np.where(deviation > threshold_std * baseline_std)[0]

    print(f"[detect_onset] {file_path if 'file_path' in dir() else ''}")
    print(f"  baseline_mean = {baseline_mean:.4f}")
    print(f"  baseline_std  = {baseline_std:.6f}")
    print(f"  threshold     = {threshold_std} × std = {threshold_std * baseline_std:.6f}")
    print(f"  max deviation = {float(np.max(deviation)):.6f}")
    print(f"  bootstrap_n   = {bootstrap_n} samples")

    if len(crossings) == 0:
        raise ValueError(
            "No onset detected — signal never exceeds threshold. "
            "Try lowering threshold_std."
        )
    return int(crossings[0])

def detect_plateau_onset(signal: np.ndarray, sample_freq_hz: float,
                         bootstrap_s: float = 0.08, threshold_std: float = 4.0,
                         smooth_window_s: float = 0.02,
                         file_path: str = "") -> int:
    """
    Detect the start of the plateau by scanning backwards from the end of
    the recording. The plateau onset is defined as the last sample where
    the signal is below the plateau band.

    This function fits the bootstrap region to a quadratic curve to accommodate 
    a concave-downward plateau shape. The lower band is then curve_fit(x) - N*residual_std
    evaluated at every sample, rather than a fixed horizontal threshold.

    Parameters
    ----------
    signal          : 1D array (force)
    sample_freq_hz  : sampling rate
    bootstrap_s     : seconds at the very end used to estimate the plateau curve
    threshold_std   : number of residual SDs below the fitted curve defining
                      the lower edge of the plateau band
    smooth_window_s : smoothing kernel width (seconds) to suppress noise
    file_path       : included in error messages for traceability

    Returns
    -------
    plateau onset index (int)
    """
    N           = len(signal)
    bootstrap_n = int(bootstrap_s * sample_freq_hz)
    boot_start  = N - bootstrap_n

    # Use zero-based x inside the fitter for numerical stability,
    # then evaluate the same relative offsets over the full signal.
    x_boot_rel = np.arange(bootstrap_n, dtype=float)
    y_boot     = signal[boot_start:]
    x_full_rel = np.arange(N, dtype=float) - boot_start  

    kernel   = max(1, int(smooth_window_s * sample_freq_hz))
    smoothed = np.convolve(signal, np.ones(kernel) / kernel, mode='same')

    def exp_decay(x, A, lam, C):
        return A * np.exp(-lam * x) + C

    A0   = float(np.max(y_boot) - np.min(y_boot))   # rough overshoot amplitude
    C0   = float(np.min(y_boot))                     # asymptotic floor
    lam0 = 1.0 / bootstrap_n                         # decays over ~bootstrap window

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            popt, _ = curve_fit(
                exp_decay, x_boot_rel, y_boot,
                p0=[A0, lam0, C0],
                bounds=([-np.inf, 0, -np.inf], [np.inf, np.inf, np.inf]),
                maxfev=10_000
            )
    except (RuntimeError, OptimizeWarning) as e:
        raise RuntimeError(
            f"{file_path}: exponential fit failed in detect_plateau_onset. "
            f"Consider adjusting bootstrap_s or inspecting the tail region. "
            f"Original error: {e}"
        )

    A_fit, lam_fit, C_fit = popt

    fitted_boot = exp_decay(x_boot_rel, *popt)
    fitted_full = exp_decay(x_full_rel, *popt)

    # Band: curve-relative, using residual spread in bootstrap region
    residual_std = np.std(y_boot - fitted_boot)
    lower_band   = fitted_full - threshold_std * residual_std

    # Last sample still below the plateau band = end of the rise
    below = np.where(smoothed < lower_band)[0]

    print(f"[detect_plateau_onset] {file_path}")
    print(f"  curve_type    = exponential decay (A·exp(–λx) + C)")
    print(f"  A             = {A_fit:.6f}")
    print(f"  λ             = {lam_fit:.6f}")
    print(f"  C             = {C_fit:.6f}")
    print(f"  residual_std  = {residual_std:.6f}")
    print(f"  threshold     = {threshold_std} × residual_std = {threshold_std * residual_std:.6f}")
    print(f"  bootstrap_n   = {bootstrap_n} samples")

    if len(below) == 0:
        raise ValueError(
            f"{file_path}: detect_plateau_onset found no samples below the plateau band. "
            f"A={A_fit:.6f}, λ={lam_fit:.6f}, C={C_fit:.6f}, "
            f"residual_std={residual_std:.6f}, "
            f"min(smoothed)={float(np.min(smoothed)):.4f}, "
            f"min(lower_band)={float(np.min(lower_band)):.4f}. "
            f"Try increasing threshold_std or bootstrap_s."
        )

    return int(below[-1])

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
      'time'             : time axis in seconds
      'length_mm'        : muscle length in mm (using calibrated AI channel)
      'force_mN'         : force in mN (using calibrated AI channel)
      'raw_df'           : full pandas DataFrame of the test data
      'start_idx_force'  : onset index detected from force signal
      'start_idx_length' : onset index detected from length signal
      'end_idx_force'    : offset index detected from force signal
      'end_idx_length'   : offset index detected from length signal
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
    force_volts = df[force_chan].to_numpy(dtype=float)

    # Convert to physical units: units = (volts - offset) * scale
    length_mm = (length_volts - calib[length_chan]["offset"]) * calib[length_chan]["scale"] + calib[length_chan]["tad"]
    force_ref = (force_volts - calib[force_chan]["offset"]) * calib[force_chan]["scale"] + calib[force_chan]["tad"]

    # Treat 'Ref' for force as mN (matches your previous scaling)
    force_mN = force_ref

    # Build time axis
    num_samples = len(df)
    sample_indices = np.arange(num_samples, dtype=float)
    time_s = sample_indices / sample_freq_hz

    # Detect onset and offset from each signal independently
    start_idx = detect_onset(force_mN,  sample_freq_hz)
    end_idx_force = int(start_idx + np.argmin(length_mm[start_idx:]))
    print("end_idx_force =",end_idx_force)
    end_idx = detect_plateau_onset(force_mN[:end_idx_force], sample_freq_hz, file_path=file_path)

    return {
        "time": time_s,
        "length_mm": length_mm,
        "force_mN": force_mN,
        "raw_df": df,
        "start_idx": start_idx,
        "end_idx": end_idx,
        "end_plateau": end_idx_force
    }

def isometric_work_from_file(file_path: str) -> float:
    """
    Compute the isometric work (force-time integral) for a single file.

    The isometric work is the area under the baseline-corrected force curve
    during the rising phase — from signal onset to plateau entry.

    Units: mN · s  (force-time integral / impulse)

    Steps:
      - Parse DMC file.
      - Baseline-correct force using the pre-contraction mean.
      - Integrate (force - baseline) over time from onset to plateau entry.
    """
    parsed        = parse_dmc_file(file_path)
    time          = parsed["time"]
    force_mN      = parsed["force_mN"]
    start_idx     = parsed["start_idx"]
    end_idx       = parsed["end_idx"]

    print("(cy) start_idx =", start_idx)
    print("(cy) end_idx   =", end_idx)

    # Baseline from pre-contraction region
    force_baseline = float(np.mean(force_mN[:start_idx]))

    force_seg = force_mN[start_idx:end_idx + 1] - force_baseline
    time_seg  = time[start_idx:end_idx + 1]

    total_work = float(np.trapezoid(force_seg, time_seg))

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
    Produce force and work validation figures for a single file.
    Overlays onset and plateau-entry markers on the force plot.
    Shades the integrated region on the work (force-time) plot.
    """
    parsed = parse_dmc_file(file_path)
    time = parsed["time"]
    force_mN = parsed["force_mN"]
    length_mm = parsed["length_mm"]
    start_idx = parsed["start_idx"]
    end_idx = parsed["end_idx"]
    end_plateau = parsed["end_plateau"]

    force_baseline = float(np.mean(force_mN[:start_idx]))
    force_seg = force_mN[start_idx:end_idx + 1] - force_baseline
    time_seg  = time[start_idx:end_idx + 1]

    short_path = "/".join(PlPath(file_path).parts[-3:])

    fig_l, ax_l = plt.subplots(figsize=(11, 8))
    ax_l.plot(time[:end_plateau + 1], force_mN[:end_plateau + 1])
    ax_l.fill_between(time[start_idx:end_idx], force_mN[start_idx:end_idx], 0,
                      alpha=0.3, label="Integrated Work")
    ax_l.plot(time[start_idx], force_mN[start_idx], 'o', color='green', label="Start Point")
    ax_l.plot(time[end_idx], force_mN[end_idx], 'o', color='brown', label="End Point")
    ax_l.set_xlabel("Time (s)")
    ax_l.set_ylabel("Force (mN))")
    ax_l.set_title(f"{short_path} (index {i})")
    ax_l.grid(True)

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

    csv_output  = []
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
    csv_output   = st.session_state.csv_output
    csv_x_labels = st.session_state.csv_x_labels
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
                pct      = st.slider("Sample fraction", 0.0, 1.0, 0.05, 0.01)
                n_samples = max(1, int(pct * len(open_animal_folder)))

                if st.button("Run Random Sample"):
                    random_indices = sorted(random.sample(
                        range(len(open_animal_folder)),
                        k=min(n_samples, len(open_animal_folder)),
                    ))
                    for idx in random_indices:
                        val_path = open_animal_folder[idx]
                        fig_l, fig_f, fig_w = val_isometric_work(val_path, idx)
                        st.pyplot(fig_l)
                        st.pyplot(fig_f)
                        st.pyplot(fig_w)

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
                            fig_l, fig_f, fig_w = val_isometric_work(val_path, idx)
                            st.pyplot(fig_l)
                            st.pyplot(fig_f)
                            st.pyplot(fig_w)