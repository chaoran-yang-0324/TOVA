"""
power.py

Description: Process all animals in a folder, compute normalized maximum instantaneous 
power for each contraction, and return a figure plus the raw results.
 [y] fixed contraction detector
 [y] add debug function
 [y] write start & end detection into max_inst function, not parse
 [y] make x axis of graph the file names
 [] add the option to save mass data
    add the option to upload mass data (from previous generation)
 [y] fix baseline_start detection 
 [] add units to the output CSV file
"""

__author__ = "Chaoran Yang"
__version__ = "2.2"
__email__ = "cy197@duke.edu"
__date__ = "2026-03-09"

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

def natural_sort_key(s: str):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r"([0-9]+)", s)]

def detect_onset(signal: np.ndarray, sample_freq_hz: float,
                 bootstrap_s: float = 0.02, threshold_std: float = 3.0,
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
    baseline_mean = np.mean(bootstrap)
    baseline_std  = np.std(bootstrap)

    kernel = max(1, int(smooth_window_s * sample_freq_hz))
    smoothed = np.convolve(signal, np.ones(kernel) / kernel, mode='same')

    deviation = np.abs(smoothed - baseline_mean)
    crossings = np.where(deviation > threshold_std * baseline_std)[0]

    if len(crossings) == 0:
        raise ValueError(
            "No onset detected — signal never exceeds threshold. "
            "Try lowering threshold_std."
        )
    return int(crossings[0])

def detect_offset(signal: np.ndarray, sample_freq_hz: float,
                  bootstrap_s: float = 0.02, threshold_std: float = 4.0,
                  smooth_window_s: float = 0.02) -> int:
    """
    Detect the last sample where `signal` is still outside its initial baseline.
    Scans from the end of the array backward to find where the signal
    re-enters the baseline band.

    Parameters
    ----------
    signal          : 1D array (force or length)
    sample_freq_hz  : sampling rate
    bootstrap_s     : seconds at the very start used to estimate baseline
    threshold_std   : number of baseline SDs required to declare offset
    smooth_window_s : smoothing kernel width (seconds) to suppress noise

    Returns
    -------
    offset index (int)
    """
    bootstrap_n = int(bootstrap_s * sample_freq_hz)
    bootstrap = signal[:bootstrap_n]
    baseline_mean = np.mean(bootstrap)
    baseline_std  = np.std(bootstrap)

    kernel = max(1, int(smooth_window_s * sample_freq_hz))
    smoothed = np.convolve(signal, np.ones(kernel) / kernel, mode='same')

    deviation = np.abs(smoothed - baseline_mean)
    crossings = np.where(deviation > threshold_std * baseline_std)[0]

    if len(crossings) == 0:
        raise ValueError(
            "No offset detected — signal never exceeds threshold. "
            "Try lowering threshold_std."
        )
    return int(crossings[-1])

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

    # Detect onset from each signal independently
    start_idx_force  = detect_onset(force_mN,  sample_freq_hz)
    start_idx_length = detect_onset(length_mm, sample_freq_hz)
    end_idx_force    = detect_offset(force_mN,  sample_freq_hz)
    end_idx_length   = detect_offset(length_mm, sample_freq_hz)

    return {
        "time":             time_s,
        "length_mm":        length_mm,
        "force_mN":         force_mN,
        "raw_df":           df,
        "start_idx_force":  start_idx_force,
        "start_idx_length": start_idx_length,
        "end_idx_force":    end_idx_force,
        "end_idx_length":   end_idx_length,
    }

def max_instantaneous_power_from_file(file_path: str) -> float:
    """
    Compute the maximum instantaneous power for a single file.

    Steps:
      - Parse DMC file (header + data).
      - Detect contraction window automatically from force trace.
      - Use a pre-contraction middle segment as baseline.
      - Compute instantaneous power: P(t) = F(t) * v(t), converted to Watts.
    """
    parsed = parse_dmc_file(file_path)
    time = parsed["time"]
    length_mm = parsed["length_mm"]
    force_mN = parsed["force_mN"]
    # start_idx_force = parsed["start_idx_force"]
    start_idx_length = parsed["start_idx_length"]
    # end_idx_force = parsed["end_idx_force"]
    end_idx_length = parsed["end_idx_length"]

    # Baseline correction using the pre-contraction middle segment
    length_baseline = float(np.mean(length_mm[0:start_idx_length-5]))  # small buffer 
    force_baseline = float(np.mean(force_mN[0:start_idx_length-5]))

    length_seg = length_mm[start_idx_length:end_idx_length] - length_baseline
    force_seg = force_mN[start_idx_length:end_idx_length] - force_baseline
    time_seg = time[start_idx_length:end_idx_length]

    # Velocity (mm/s)
    velocity_mm_s = np.gradient(length_seg, time_seg)

    # Power: force (mN) * velocity (mm/s) -> W via 1e-6 factor
    inst_power_W = - 1e-6 * force_seg * velocity_mm_s
    min_idx = np.argmin(inst_power_W)
    inst_power_sliced = inst_power_W[:min_idx + 1]
    max_power_W = float(np.max(inst_power_sliced))

    return max_power_W

def run_max_inst_power(folder_path: str,
                       mass_kg: float) -> List[List[float]]:
    """
    Process one animal folder under `folder_path` and compute normalized
    peak instantaneous power (W/kg) for each contraction.

    Parameters
    ----------
    folder_path : str
        Root directory of the one animal folder.
    mass_kg : float
        Animal mass in kilograms used to normalize peak power.

    Returns
    -------
    outputs : list[list[float]]
        Normalized peak power values (W/kg) for each contraction and animal.
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

        peak_power_W = max_instantaneous_power_from_file(data_file)
        normalized_power = peak_power_W / mass_kg  # W/kg
        animal_results.append(normalized_power)

    return animal_results

def val_max_inst_power(file_path: str,
                       i: float) -> plt:
    """
    Selects specific .ddf files and shows the power vs. time plot
    to validate run_max_inst_power results 

    Steps:
      - Parse DMC file (header + data).
      - Detect contraction window automatically from force trace.
      - Use a pre-contraction middle segment as baseline.
      - Compute instantaneous power: P(t) = F(t) * v(t), converted to Watts.
    """
    parsed = parse_dmc_file(file_path)
    time = parsed["time"]
    length_mm = parsed["length_mm"]
    force_mN = parsed["force_mN"]
    # start_idx_force = parsed["start_idx_force"]
    start_idx_length = parsed["start_idx_length"]
    # end_idx_force = parsed["end_idx_force"]
    end_idx_length = parsed["end_idx_length"]

    # Baseline correction using the pre-contraction middle segment
    length_baseline = float(np.mean(length_mm[0:start_idx_length-5]))  # small buffer 
    force_baseline = float(np.mean(force_mN[0:start_idx_length-5]))

    length_seg = length_mm[start_idx_length:end_idx_length] - length_baseline
    force_seg = force_mN[start_idx_length:end_idx_length] - force_baseline
    time_seg = time[start_idx_length:end_idx_length]

    # Velocity (mm/s)
    velocity_mm_s = np.gradient(length_seg, time_seg)

    # Power: force (mN) * velocity (mm/s) -> W via 1e-6 factor
    inst_power_W = - 1e-6 * force_seg * velocity_mm_s
    min_idx = np.argmin(inst_power_W)
    inst_power_sliced = inst_power_W[:min_idx + 1]
    max_power_W = float(np.max(inst_power_sliced))
    max_idx = np.argmax(inst_power_sliced)

    fig_l, ax_l = plt.subplots(figsize=(11, 8))
    ax_l.plot(time_seg[:min_idx + 1], length_seg[:min_idx + 1])
    ax_l.set_xlabel("Time (s)")
    ax_l.set_ylabel("Length (mm))")
    ax_l.set_title(f"{file_path} (index {i}) Length")
    ax_l.grid(True)

    fig_f, ax_f = plt.subplots(figsize=(11, 8))
    ax_f.plot(time_seg[:min_idx + 1], force_seg[:min_idx + 1])
    ax_f.set_xlabel("Time (s)")
    ax_f.set_ylabel("Force (mN)")
    ax_f.set_title(f"{file_path} (index {i}) Force")
    ax_f.grid(True)

    fig_p, ax_p = plt.subplots(figsize=(11, 8))
    ax_p.plot(time_seg[:min_idx + 1], inst_power_W[:min_idx + 1], label="Power")
    ax_p.plot(time_seg[max_idx], max_power_W, marker="o", label="Max Power")
    ax_p.set_xlabel("Time (s)")
    ax_p.set_ylabel("Power (W)")
    ax_p.legend()
    ax_p.set_title(f"{file_path} (index {i}) Power")
    ax_p.grid(True)

    return fig_l, fig_f, fig_p

st.title("Peak Power Analysis")

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
           |
            -> ...
       -> animal_3
           |
            -> ...
       -> ...

- don't change the name of your file when compressing it
"""

st.code(folder_structure, language='text')

if uploaded_zip:
    # Define the folder where we will save the uploaded file
    upload_folder = "uploaded_files"

    # Create the folder if it doesn't exist
    if not os.path.exists(upload_folder):
        os.makedirs(upload_folder)

    # Save the uploaded zip file
    zip_file_path = os.path.join(upload_folder, uploaded_zip.name)

    with open(zip_file_path, "wb") as f:
        f.write(uploaded_zip.getbuffer())

    st.success(f"File {uploaded_zip.name} uploaded successfully!")

    # Unzip the file to a folder
    unzip_folder = os.path.join(upload_folder, uploaded_zip.name.split('.')[0])

    # with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        # zip_ref.extractall(unzip_folder)

    with zipfile.ZipFile(uploaded_zip, 'r') as zip_ref:
        for member in zip_ref.namelist():
            if '__MACOSX' not in member:
                zip_ref.extract(member, unzip_folder)

    unzip_folder = os.path.join(unzip_folder, uploaded_zip.name.split('.')[0])

    st.write(f"Contents of {unzip_folder}:")
    sorted_folder_names = sorted(os.listdir(unzip_folder), key=natural_sort_key)
    st.write(sorted_folder_names)  # Display contents

    sorted_folder_names = sorted(os.listdir(unzip_folder), key=natural_sort_key)
    mass_g: List[float] = []

    st.subheader("Animal Mass")
    mass_mode = st.radio("Mass input mode", ["Manual Input", "Upload CSV"], horizontal=True)

    if mass_mode == "Manual Input":
        for i, foldername in enumerate(sorted_folder_names):
            value = st.number_input(f"{foldername} Mass (g):", min_value=0.0, value=1.0, step=0.00001, key=f"mass_{i}")
            mass_g.append(value)

        # Build and offer download of current mass values
        mass_df = pd.DataFrame({
            "folder": sorted_folder_names,
            "mass_g": mass_g,
        })
        st.download_button(
            label="Save mass data as CSV",
            data=mass_df.to_csv(index=False).encode("utf-8"),
            file_name="mass_data.csv",
            mime="text/csv",
        )

    else:  # Upload CSV
        mass_csv = st.file_uploader("Upload mass CSV", type="csv")
        if mass_csv is not None:
            uploaded_mass_df = pd.read_csv(mass_csv)
            if set(["folder", "mass_g"]).issubset(uploaded_mass_df.columns):
                uploaded_mass_df = uploaded_mass_df.set_index("folder")
                missing = [f for f in sorted_folder_names if f not in uploaded_mass_df.index]
                if missing:
                    st.error(f"These folders have no mass entry in the uploaded CSV: {missing}")
                else:
                    for foldername in sorted_folder_names:
                        mass_g.append(float(uploaded_mass_df.loc[foldername, "mass_g"]))
                    st.success("Mass data loaded from CSV.")
                    st.dataframe(uploaded_mass_df.loc[sorted_folder_names])
            else:
                st.error('CSV must have columns "folder" and "mass_g".')

if "analysis_done" not in st.session_state:
    st.session_state.analysis_done = False
if "csv_output" not in st.session_state:
    st.session_state.csv_output = None
if "animal_folders" not in st.session_state:
    st.session_state.animal_folders = None
if "csv_name" not in st.session_state:
    st.session_state.csv_name = None
if "csv_x_labels" not in st.session_state:
    st.session_state.csv_x_labels = None

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

    csv_output = []
    csv_x_labels = []

    for i, foldername in enumerate(sorted(os.listdir(unzip_folder), key=natural_sort_key)):
        run_path = os.path.join(unzip_folder, foldername)
        value = run_max_inst_power(run_path, mass_kg=mass_g[i] * 0.001)
        csv_output.append(value)

        csv_x_labels.append([])  # new row for this run
        for filename in sorted(os.listdir(run_path), key=natural_sort_key):
            if filename.lower().endswith(".ddf"):
                csv_x_labels[i].append(filename)

    # store results so they survive reruns
    st.session_state.csv_output = csv_output
    st.session_state.csv_x_labels = csv_x_labels   
    st.session_state.animal_folders = list_subfolders(unzip_folder)

    now = datetime.now()
    st.session_state.csv_name = f"peak_power_{now.strftime('%Y_%m_%d_%H_%M_%S')}.csv"
    st.session_state.analysis_done = True

if st.session_state.analysis_done:
    csv_output = st.session_state.csv_output
    csv_x_labels = st.session_state.csv_x_labels   
    animal_folders = st.session_state.animal_folders or []

    # main plot
    st.write("Graphing...")
    fig, ax = plt.subplots(figsize=(11, 8))
    for idx, result in enumerate(csv_output):
        x_coord = np.arange(len(result))
        ax.plot(x_coord, np.array(result), label=f"Folder {idx + 1}")

    ticks = np.arange(len(result))
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks)
    ax.set_xlabel("Contraction Index")
    ax.set_ylabel("Normalized Power (W/kg)")
    ax.set_title("Peak Power")
    ax.legend()
    ax.grid(True)
    st.pyplot(fig)

    # download
    df = pd.DataFrame(csv_output)

    rows = []
    for i in range(len(csv_x_labels)):
        labels = csv_x_labels[i]
        values = df.iloc[i].tolist()

        for fname, val in zip(labels, values):
            rows.append([fname, val])

        rows.append(["", ""])  # blank line between runs

    final_df = pd.DataFrame(rows, columns=["filename", "peak power"])
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
        # pick which animal folder to inspect
        animal_choice = st.selectbox(
            "Select animal folder",
            options=animal_folders,
            format_func=lambda p: os.path.basename(p),
        )

        open_animal_folder = [f for f in list_files(animal_choice)
                              if f.lower().endswith(".ddf")]

        if len(open_animal_folder) == 0:
            st.warning("No files found in the selected animal folder.")
        else:
            mode = st.radio(
                "Validation mode",
                ["Random Sample", "Specific Inspection"],
                horizontal=True,
            )

            if mode == "Random Sample":
                pct = st.slider("Sample fraction", 0.0, 1.0, 0.05, 0.01)
                n_samples = max(1, int(pct * len(open_animal_folder)))

                if st.button("Run Random Sample"):
                    random_indices = sorted(random.sample(
                        range(len(open_animal_folder)),
                        k=min(n_samples, len(open_animal_folder)),
                    ))
                    for i in random_indices:
                        val_path = open_animal_folder[i]
                        fig_l, fig_f, fig_p = val_max_inst_power(val_path, i)
                        st.pyplot(fig_l)
                        st.pyplot(fig_f)
                        st.pyplot(fig_p)

            else:  # Specific Inspection
                start = st.number_input(
                    "Start index (inclusive)",
                    min_value=0,
                    max_value=len(open_animal_folder) - 1,
                    value=0,
                    step=1,
                )
                end = st.number_input(
                    "End index (inclusive)",
                    min_value=0,
                    max_value=len(open_animal_folder) - 1,
                    value=len(open_animal_folder) - 1,
                    step=1,
                )

                if start > end:
                    st.error("Start index must be less than or equal to end index.")
                else:
                    if st.button("Run Specific Inspection"):
                        for i in range(int(start), int(end) + 1):
                            val_path = open_animal_folder[i]
                            fig_l, fig_f, fig_p = val_max_inst_power(val_path, i)
                            st.pyplot(fig_l)
                            st.pyplot(fig_f)
                            st.pyplot(fig_p)
